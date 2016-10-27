module PLX
using ..Struct, ..Common, PiecewiseIncreasingRanges
import Compat.String

export PLXSpikeChannel, PLXEventChannel, PLXContinuousChannel, PLXFile, readplx

const MagicNumber = 0x58454c50

function struct_string(fn::Function, bytes::Vector{UInt8})
    last_byte = findfirst(bytes, 0)
    fn(last_byte == 0 ? bytes : bytes[1:last_byte-1])
end

@struct PL_FileHeader begin
    MagicNumber::UInt32
    Version::Int32
    Comment::String(128)
    ADFrequency::Int32
    NumDSPChannels::Int32
    NumEventChannels::Int32
    NumSlowChannels::Int32
    NumPointsWave::Int32
    NumPointsPreThr::Int32

    Year::Int32
    Month::Int32
    Day::Int32
    Hour::Int32
    Minute::Int32
    Second::Int32

    FastRead::Int32
    WaveformFreq::Int32
    LastTimestamp::Float64

    Trodalness::UInt8
    DataTrodalness::UInt8
    BitsPerSpikeSample::UInt8
    BitsPerSlowSample::UInt8
    SpikeMaxMagnitudeMV::UInt16
    SlowMaxMagnitudeMV::UInt16
    SpikePreAmpGain::UInt16

    Padding::UInt8(46)

    TSCounts::Int32(5, 130)
    WFCounts::Int32(5, 130)

    EVCounts::Int32(512)
end

@struct PL_ChanHeader begin
    Name::String(32)
    SIGName::String(32)
    Channel::Int32
    WFRate::Int32
    SIG::Int32
    Ref::Int32
    Gain::Int32
    Filter::Int32
    Threshold::Int32
    Method::Int32
    NUnits::Int32
    Template::Int16(64, 5)
    Fit::Int32(5)
    SortWidth::Int32
    Boxes::Int16(4, 2, 5)
    SortBeg::Int32
    Comment::String(128)
    Padding::Int32(11)
end

@struct PL_EventHeader begin
    Name::String(32)
    Channel::Int32
    Comment::String(128)
    Padding::Int32(33)
end

@struct PL_SlowChannelHeader begin
    Name::String(32)
    Channel::Int32
    ADFreq::Int32
    Gain::Int32
    Enabled::Int32
    PreAmpGain::Int32

    SpikeChannel::Int32

    Comment::String(128)
    Padding::Int32(28)
end

type PLXSpikeChannel <: SpikeChannel
    header::PL_ChanHeader
    sample_rate::Int32
    voltage_multiplier::Float64
    bits_per_sample::Int8
    times::Vector{Float64}
    unit_numbers::Vector{Int16}
    data::Matrix{Int16}

    function PLXSpikeChannel(header::PL_ChanHeader, fheader::PL_FileHeader)
        voltage_multiplier = 1e-3*(if fheader.Version >= 105
                fheader.SpikeMaxMagnitudeMV/
                    ((1 << (fheader.BitsPerSpikeSample-1))*header.Gain*fheader.SpikePreAmpGain*1000)
            elseif fheader.Version >= 103
                fheader.SpikeMaxMagnitudeMV/
                    ((1 << (fheader.BitsPerSpikeSample-1))*header.Gain*1000)
            else
                3000/(2048*header.Gain*1000)
            end)
        return new(header, fheader.ADFrequency, voltage_multiplier, fheader.BitsPerSpikeSample, Float64[], Int16[])
    end
end

type PLXEventChannel <: EventChannel
    header::PL_EventHeader
    times::Vector{Float64}
    data::Vector{Int16}

    PLXEventChannel(header::PL_EventHeader) = new(header, Float64[])
end

type PLXContinuousChannel <: ContinuousChannel
    header::PL_SlowChannelHeader
    voltage_multiplier::Float64
    data::Vector{Int16}
    times::PiecewiseIncreasingRange{Float64,StepRange{Int64,Int64},Int32}

    function PLXContinuousChannel(header::PL_SlowChannelHeader, fheader::PL_FileHeader)
        voltage_multiplier = 1e-3*(if fheader.Version >= 103
                fheader.SlowMaxMagnitudeMV/
                    ((1 << (fheader.BitsPerSpikeSample-1))*header.Gain*header.PreAmpGain)
            elseif fheader.Version == 102
                5000/(2048*header.Gain*header.PreAmpGain)
            else
                5000/(2048*1000*header.Gain)
            end)

        return new(header, voltage_multiplier)
    end
end

Common.label(ch::Union{PLXSpikeChannel,PLXEventChannel,PLXContinuousChannel}) =
    ch.header.Name

immutable PLXFile <: File
    header::PL_FileHeader
    spike_channels::Channels{PLXSpikeChannel}
    event_channels::Channels{PLXEventChannel}
    continuous_channels::Channels{PLXContinuousChannel}
end
PLXFile(header::PL_FileHeader) =
    PLXFile(header, Channels{PLXSpikeChannel}(), Channels{PLXEventChannel}(),
            Channels{PLXContinuousChannel}())

function readplx(ios::IOStream; lfps::Bool=true, waveforms::Bool=true)
    const maxChannels = 32*1024

    fheader = read(ios, PL_FileHeader)
    if fheader.MagicNumber != MagicNumber
        error("$file_name does not appear to be a PLX file")
    elseif fheader.NumDSPChannels < 0 || fheader.NumDSPChannels > maxChannels ||
            fheader.NumEventChannels < 0 || fheader.NumEventChannels > maxChannels ||
            fheader.NumSlowChannels < 0 || fheader.NumSlowChannels > maxChannels
        error("PLX file header specifies an invalid number of channels")
    elseif fheader.ADFrequency <= 0
        error("PLX file header specifies ADFrequency <= 0")
    end
    x = PLXFile(fheader)

    for i = 1:x.header.NumDSPChannels
        header = read(ios, PL_ChanHeader)
        x.spike_channels[header.Channel] = PLXSpikeChannel(header, fheader)
    end
    if x.header.NumDSPChannels > 0 && fheader.NumPointsWave <= 0
        error("PLX file header specifies no points in waveform")
    end

    for i=1:x.header.NumEventChannels
        header = read(ios, PL_EventHeader)
        x.event_channels[header.Channel] = PLXEventChannel(header)
    end

    for i=1:x.header.NumSlowChannels
        header = read(ios, PL_SlowChannelHeader)
        lfps && (x.continuous_channels[header.Channel+1] = PLXContinuousChannel(header, fheader))
    end

    data_offset = position(ios)

    # Read through the file once to determine how much memory to allocate
    nspikes = zeros(Int, endof(x.spike_channels))
    nevents = zeros(Int, endof(x.event_channels))
    ntimestamps = zeros(Int, endof(x.continuous_channels))
    nsamples = zeros(Int, endof(x.continuous_channels))

    nblocks = 0
    max_offset = div(filesize(ios)-data_offset, 2)
    contents = Mmap.mmap(ios, Vector{Int16}, max_offset)::Vector{Int16}
    cur_offset = 1
    while cur_offset < max_offset
        block_type = contents[cur_offset]
        ch = contents[cur_offset+4]

        if block_type == 1          # spike
            nspikes[ch] += 1
        elseif block_type == 4      # event
            nevents[ch] += 1
        elseif block_type == 5      # continuous
            if lfps
                ntimestamps[ch+1] += 1
                nsamples[ch+1] += contents[cur_offset+7]
            end
        else
            error("invalid data block type ", block_type)
        end

        cur_offset += contents[cur_offset+7]+8
        nblocks += 1
    end

    # Allocate
    for (i, channel) = enumerate(x.spike_channels)
        resize!(channel.times, nspikes[i])
        resize!(channel.unit_numbers, nspikes[i])
        if waveforms
            channel.data = Array(Int16, fheader.NumPointsWave, nspikes[i])
        end
    end
    for (i, channel) = enumerate(x.event_channels)
        resize!(channel.times, nevents[i])
        if i == 257
            channel.data = Array(Int16, nevents[i])
        end
    end

    if lfps
        sample_dts = Array(Int64, length(ntimestamps))
        timestamp_ranges = Array(Vector{StepRange{Int64,Int64}}, length(ntimestamps))
        for (i, channel) = enumerate(x.continuous_channels)
            sample_dt = div(fheader.ADFrequency, channel.header.ADFreq)
            if sample_dt * channel.header.ADFreq != fheader.ADFrequency
                warn("""
                    channel $i frequency $(channel.header.ADFreq) is non-integer
                    multiple of AD frequency $(x.header.ADFrequency)
                """)
            end
            sample_dts[i] = sample_dt 
            timestamp_ranges[i] = Array(StepRange{Int64,Int64}, ntimestamps[i])
            channel.data = Array(Int16, nsamples[i])
        end
    else
        sample_dts = Array(Int64, 0)
        timestamp_ranges = Array(Vector{StepRange{Int64,Int64}}, 0)
    end

    # Read through the file again
    cur_spike = zeros(Int, size(nspikes))
    cur_event = zeros(Int, size(nevents))
    cur_timestamp = zeros(Int, size(ntimestamps))
    cur_sample = zeros(Int, size(nsamples))
    cur_offset = 1
    while cur_offset < max_offset
        block_type = contents[cur_offset]
        timestamp = (convert(Int64, reinterpret(UInt16, contents[cur_offset+1])) << 32) +
            reinterpret(UInt16, contents[cur_offset+2]) +
            (convert(Int64, reinterpret(UInt16, contents[cur_offset+3])) << 16)
        ch = convert(Int, contents[cur_offset+4])

        if block_type == 1          # spike
            let channel = x.spike_channels[ch]
                num = (cur_spike[ch] += 1)

                channel.unit_numbers[num] = contents[cur_offset+5]
                channel.times[num] = timestamp/x.header.ADFrequency

                block_samples = contents[cur_offset+6]*contents[cur_offset+7]
                if waveforms
                    wf = channel.data
                    if block_samples == 0
                        for i = 1:size(wf, 1)
                            wf[i, num] = 0
                        end
                    else
                        for i = 1:block_samples
                            wf[i, num] = contents[cur_offset+7+i]
                        end
                    end
                end
                cur_offset += block_samples
            end
        elseif block_type == 4      # event
            let channel = x.event_channels[ch]
                num = (cur_event[ch] += 1)

                channel.times[num] = timestamp/x.header.ADFrequency
                if ch == 257
                    channel.data[num] = contents[cur_offset+5]
                end
            end
        elseif block_type == 5      # continuous
            block_samples = contents[cur_offset+7]
            if block_samples > 0
                if lfps
                    t = (cur_timestamp[ch+1] += 1)
                    channel = x.continuous_channels[ch+1]
                    num = cur_sample[ch+1]

                    sample_dt = sample_dts[ch+1]
                    timestamp_ranges[ch+1][t] = timestamp:sample_dt:timestamp+sample_dt*(block_samples-1)

                    start_offset = cur_offset + 7
                    samples = channel.data
                    for i = 1:block_samples
                        samples[num+i] = contents[start_offset+i]
                    end
                    cur_sample[ch+1] += block_samples
                end
                cur_offset += block_samples
            end
        end
        cur_offset += 8
    end

    @assert cur_spike == nspikes
    @assert cur_event == nevents
    @assert cur_timestamp == ntimestamps
    @assert cur_sample == nsamples

    if lfps
        for (i, channel) = enumerate(x.continuous_channels)
            channel.times = PiecewiseIncreasingRange(timestamp_ranges[i], x.header.ADFrequency)
        end
    end

    close(ios)

    return x
end

chmag(ch) = round(UInt16, 2^(bitspersample(ch)-1)*voltagemultiplier(ch)*1000)

function PLXFile(;spike_channels=Channels{PLXSpikeChannel}(),
                  event_channels=Channels{PLXEventChannel}(),
                  continuous_channels=Channels{PLXContinuousChannel}())
    ADFrequency = 0
    NumPointsWave = 0
    LastTimestamp = 0.
    BitsPerSpikeSample = 16
    BitsPerSlowSample = 16
    SpikeMaxMagnitudeMV::UInt16 = 0
    SlowMaxMagnitudeMV::UInt16 = 0
    if !isempty(spike_channels)
        ADFrequency = round(Int32, samplerate(first(spike_channels)))
        for ch in spike_channels
            samplerate(ch) == ADFrequency || throw(ArgumentError("spike channels do not have matching sample rates"))
            if hasdata(ch)
                NumPointsWave = size(data(ch), 1)
            end
            LastTimestamp = max(LastTimestamp, last(times(ch)))
            mag = chmag(ch)
            if SpikeMaxMagnitudeMV == 0
                SpikeMaxMagnitudeMV = mag
            elseif SpikeMaxMagnitudeMV != mag
                SpikeMaxMagnitudeMV = lcm(SpikeMaxMagnitudeMV, mag)
            end
        end
    elseif !isempty(continuous_channels)
        ADFrequency = round(Int32, samplerate(first(continuous_channels)))
        for ch in continuous_channels
            ADFrequency = max(ADFrequency, samplerate(ch))
            LastTimestamp = max(LastTimestamp, last(times(ch)))
            mag = chmag(ch)
            if SlowMaxMagnitudeMV == 0
                SlowMaxMagnitudeMV = mag
            elseif SlowMaxMagnitudeMV != mag
                SlowMaxMagnitudeMV = lcm(SlowMaxMagnitudeMV, mag)
            end
        end
    end
    LastTimestamp *= ADFrequency

    header = PL_FileHeader(MagicNumber, 103, "Created with NeuroIO.jl", ADFrequency,
                           length(spike_channels), length(event_channels), length(continuous_channels),
                           NumPointsWave, div(NumPointsWave, 2), 0, 0, 0, 0, 0, 0, 0, ADFrequency,
                           LastTimestamp, 1, 1, BitsPerSpikeSample, BitsPerSlowSample,
                           SpikeMaxMagnitudeMV, SlowMaxMagnitudeMV, 1, zeros(UInt8, 46),
                           zeros(Int32, 5, 130), zeros(Int32, 5, 130), zeros(Int32, 512))

    new_spike_channels = Channels{PLXSpikeChannel}()
    for (i, ch) in enumerate(spike_channels)
        if isa(ch, PLXSpikeChannel)
            new_spike_channels[i] = ch
            continue
        end
        Gain = round(Int32, SpikeMaxMagnitudeMV/chmag(ch))
        NUnits = maximum(unitnumbers(ch))
        chheader = PL_ChanHeader(label(ch), "", i, 10, i, 0, Gain, 0, 0, 2, NUnits, zeros(Int16, 64, 5),
                                 zeros(Int32, 5), NumPointsWave, zeros(Int16, 4, 2, 5), 0, "", zeros(Int32, 11))
        newch = new_spike_channels[i] = PLXSpikeChannel(chheader, header)
        newch.times = times(ch)
        newch.unit_numbers = unitnumbers(ch)
        hasdata(ch) && (newch.data = data(ch))
    end

    new_event_channels = Channels{PLXEventChannel}()
    have_codes = false
    for (i, ch) in enumerate(event_channels)
        if isa(ch, PLXEventChannel)
            new_event_channels[i] = ch
            continue
        end

        if hasdata(ch) && !isempty(times(ch))
            if !have_codes
                i = 257
                have_codes = true
            else
                warn("multiple event channels have codes; ignoring codes for channel $i")
            end
        end
        chheader = PL_EventHeader(label(ch), i, "", zeros(Int332, 33))
        newch = new_event_channels[i] = PLXEventChannel(chheader)
        newch.times = times(ch)
        hasdata(ch) && i == 257 && (newch.data = data(ch))
    end

    new_continuous_channels = Channels{PLXContinuousChannel}()
    chs = Int[]
    for (i, ch) in enumerate(continuous_channels)
        hasdata(ch) && push!(chs, i)
    end
    dvecs = datavecs(continuous_channels, chs)

    idatavec = 0
    for (i, ch) in enumerate(continuous_channels)
        if isa(ch, PLXContinuousChannel)
            new_continuous_channels[i] = ch
            continue
        end

        Gain = round(Int32, SlowMaxMagnitudeMV/chmag(ch))
        chheader = PL_SlowChannelHeader(label(ch), i, samplerate(ch), Gain, hasdata(ch), 1, 0, "", zeros(Int32, 28))
        newch = new_continuous_channels[i] = PLXContinuousChannel(chheader, header)
        if hasdata(ch)
            newch.times = times(ch)
            newch.data = dvecs[idatavec += 1]
        end
    end

    validate_headers(PLXFile(header, new_spike_channels, new_event_channels, new_continuous_channels))
end

function validate_headers(x::PLXFile)
    x.header.NumDSPChannels = length(x.spike_channels)
    x.header.NumEventChannels = length(x.event_channels)
    x.header.NumSlowChannels = length(x.continuous_channels)

    nwfpoints = 0
    for (i, ch) in enumerate(x.spike_channels)
        length(ch.times) == length(ch.unit_numbers) ||
            throw(ArgumentError("number of spike times and unit numbers do not match for spike channel $i"))
        if isdefined(ch, :data)
            length(ch.times) == size(ch.data, 2) ||
                throw(ArgumentError("number of spike times and waveforms do not match for spike channel $i"))
            if nwfpoints == 0
                nwfpoints = size(ch.data, 1)
            elseif nwfpoints != size(ch.data, 1)
                throw(ArgumentError("all channels must have the same number of points per waveform"))
            end
        end

        ch.header.NUnits = maximum(ch.unit_numbers)
        ch.header.Channel = i
        x.header.NumPointsWave = nwfpoints

        if i <= 129
            counts = x.header.WFCounts
            for unumber in ch.unit_numbers
                unumber > 4 && continue
                counts[unumber+1, i] += 1
            end
        end
    end

    for (i, ch) in enumerate(x.event_channels)
        if isdefined(ch, :data)
            if i == 257
                length(ch.times) == length(ch.data) ||
                    throw(ArgumentError("number of event times and codes do not match"))
            else
                throw(ArgumentError("event codes were present on channel $i; only channel 257 may have associated event codes"))
            end
        end
        ch.header.EVCounts[300+i] = length(ch.times)
        ch.header.Channel = i
    end

    for (i, ch) in enumerate(x.continuous_channels)
        length(ch.times) == length(ch.data) ||
            throw(ArgumentError("number of samples and times do not match for spike $i"))
        ch.header.Channel = i-1
    end
    x
end

function integer_timestamp(time::Float64, ADFrequency::Int32)
    timestamp = round(Int64, time*ADFrequency)
    (Int16(timestamp >>> 32), timestamp % Int16, (timestamp >>> 16) % Int16)
end

const EMPTY_WAVEFORMS = Array(Int16, 0, 0)
const EMPTY_CODES = Int16[]
function Base.write(io::IO, x::PLXFile)
    validate_headers(x)

    write(io, x.header)
    for ch in x.spike_channels
        write(io, ch.header)
    end
    for ch in x.event_channels
        write(io, ch.header)
    end
    for ch in x.continuous_channels
        write(io, ch.header)
    end

    ADFrequency = x.header.ADFrequency
    for (ich, ch) = enumerate(x.spike_channels)
        has_waveforms = isdefined(ch, :data)
        waveforms = has_waveforms ? ch.data : EMPTY_WAVEFORMS
        waveform_size = has_waveforms ? size(waveforms, 1) : 0
        blockhead = Array(Int16, 8+waveform_size)
        blockhead[1] = 1
        blockhead[5] = ich
        blockhead[7] = has_waveforms
        blockhead[8] = waveform_size

        times = ch.times
        unit_numbers = ch.unit_numbers
        for i = 1:length(times)
            (blockhead[2], blockhead[3], blockhead[4]) = integer_timestamp(times[i], ADFrequency)
            blockhead[6] = unit_numbers[i]
            if has_waveforms
                for iwf = 1:waveform_size
                    blockhead[8+iwf] = waveforms[iwf, i]
                end
            end
            write(io, blockhead)
        end
    end

    blockhead = zeros(Int16, 8)
    for (ich, ch) = enumerate(x.event_channels)
        has_codes = isdefined(ch, :data)
        codes = has_codes ? ch.data : EMPTY_CODES
        blockhead[1] = 4
        blockhead[5] = ich
        if !has_codes
            blockhead[6] = 0
        end

        times = ch.times
        for i = 1:length(times)
            (blockhead[2], blockhead[3], blockhead[4]) = integer_timestamp(times[i], ADFrequency)
            has_codes && (blockhead[6] = codes[i])
            write(io, blockhead)
        end
    end

    blockhead[1] = 5
    blockhead[6] = 0
    blockhead[7] = 1
    for (ich, ch) = enumerate(x.continuous_channels)
        times = ch.times
        blockhead[5] = ich-1
        fac = ADFrequency/times.divisor
        n = 0
        for irg = 1:length(times.ranges)
            rg = times.ranges[irg]
            for irgstart = 1:typemax(Int16):length(rg)
                dur = min(typemax(Int16), length(rg)-irgstart+1)
                timestamp = round(Int, rg[irgstart]*fac)
                blockhead[2] = Int16(timestamp >>> 32)
                blockhead[3] = timestamp % Int16
                blockhead[4] = (timestamp >>> 16) % Int16
                blockhead[8] = dur
                write(io, blockhead)
                write(io, sub(ch.data, n+1:n+dur))
                n += dur
            end
        end
    end
end
end