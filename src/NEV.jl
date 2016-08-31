module NEV
using ..Struct, ..Common, PiecewiseIncreasingRanges
import Compat.String
export NEVFile, NSXFile, readnev, readnsx, readnsx_mmap

@struct NEVHeader begin
    file_spec::UInt16
    additional_flags::UInt16
    bytes_in_headers::UInt32
    bytes_in_data_packets::UInt32
    time_resolution_of_time_stamps::UInt32
    time_resolution_of_samples::UInt32
    time_origin::UInt128
    application_to_create_file::String(32)
    comment_field::String(256)
end

@struct FilterDescription begin
    cutoff::UInt32
    order::UInt32
    filter_type::UInt16
end

type NEVSpikeChannel <: SpikeChannel
    label::String

    times::Vector{Float64}
    unit_numbers::Vector{UInt8}
    data::Matrix{Int16}

    physical_connector::Char
    connector_pin::UInt8

    digitization_factor::UInt16
    energy_threshold::UInt16
    high_threshold::Int16
    low_threshold::Int16
    number_of_sorted_units::UInt8
    spike_width::UInt16

    high_pass::FilterDescription
    low_pass::FilterDescription

    bits_per_sample::UInt8
    voltage_multiplier::Float64

    NEVSpikeChannel() = new("", Float64[], UInt8[])
end
samplerate(::NEVSpikeChannel) = 30000

type NEVEventChannel <: EventChannel
    label::String
    times::Vector{Float64}
    data::Vector{UInt16}

    NEVEventChannel() = new("", Float64[], UInt16[])
end

type NEVFile <: File
    header::NEVHeader
    digital_channels::Channels{NEVEventChannel}
    spike_channels::Channels{NEVSpikeChannel}
    array_name::String
    map_file::String

    NEVFile(header) = new(header, Channels{NEVEventChannel}(), Channels{NEVSpikeChannel}())
end

const NEURALEV = reinterpret(UInt64, "NEURALEV".data)[1]
const ARRAYNME = reinterpret(UInt64, "ARRAYNME".data)[1]
const MAPFILE = reinterpret(UInt64, "MAPFILE\0".data)[1]
const NEUEVWAV = reinterpret(UInt64, "NEUEVWAV".data)[1]
const NEUEVLBL = reinterpret(UInt64, "NEUEVLBL".data)[1]
const NEUEVFLT = reinterpret(UInt64, "NEUEVFLT".data)[1]
const DIGLABEL = reinterpret(UInt64, "DIGLABEL".data)[1]
const NEURALCD = reinterpret(UInt64, "NEURALCD".data)[1]
const CC = reinterpret(UInt16, "CC".data)[1]

function readstring(io::IO, n)
    data = read(io, UInt8, n)
    String(data[1:findfirst(data, 0)-1])
end

"""
    readnev(io::IO; waveforms=true)

Read a Blackrock Microsystems NEV (spike/events) file. The `waveforms` keyword
argument determines whether waveforms are read or ignored.
"""
function readnev(io::IO; waveforms::Bool=true)
    read(io, UInt64) == NEURALEV || error("this does not appear to be an NEV file")
    header = read(io, NEVHeader)
    bytes_in_data_packets = header.bytes_in_data_packets
    time_resolution_of_time_stamps = header.time_resolution_of_time_stamps
    nev = NEVFile(header)
    spike_channels = nev.spike_channels
    serial_digital_channel = NEVEventChannel()
    parallel_digital_channel = NEVEventChannel()

    # Extended headers
    n_extended_headers = read(io, UInt32)
    for i = 1:n_extended_headers
        packet_id = read(io, UInt64)
        if packet_id == ARRAYNME
            nev.array_name = readstring(io, 24)
        elseif packet_id == MAPFILE
            nev.map_file = readstring(io, 24)
        elseif packet_id == NEUEVWAV || packet_id == NEUEVLBL || packet_id == NEUEVFLT
            electrode_id::Int = read(io, UInt16)
            channel = isassigned(spike_channels, electrode_id) ?
                      spike_channels[electrode_id] :
                      (spike_channels[electrode_id] = NEVSpikeChannel())

            if packet_id == NEUEVWAV
                channel.physical_connector = 'A'-1+read(io, UInt8)
                channel.connector_pin = read(io, UInt8)
                channel.digitization_factor = read(io, UInt16)
                channel.energy_threshold = read(io, UInt16)
                channel.high_threshold = read(io, Int16)
                channel.low_threshold = read(io, Int16)
                channel.number_of_sorted_units = read(io, UInt8)
                channel.bits_per_sample = read(io, UInt8)*8
                channel.spike_width = read(io, UInt16)
                channel.voltage_multiplier = channel.digitization_factor/10^9
                skip(io, 8)
            elseif packet_id == NEUEVLBL
                channel.label = readstring(io, 16)
                skip(io, 6)
            elseif packet_id == NEUEVFLT
                channel.high_pass = read(io, FilterDescription)
                channel.low_pass = read(io, FilterDescription)
                skip(io, 2)
            end
        elseif packet_id == DIGLABEL
            label::String = readstring(io, 16)
            parallel::Bool = read(io, UInt8)
            if parallel
                parallel_digital_channel.label = label
            else
                serial_digital_channel.label = label
            end
            skip(io, 7)
        else
            skip(io, 24)
        end
    end

    waveform_bytes = bytes_in_data_packets-8
    waveform_samples = div(waveform_bytes, 2)
    waveform_tmp = Array(UInt8, waveform_bytes)
    if waveforms
        electrode_waveforms = [UInt8[] for i = 1:length(spike_channels)]
    else
        electrode_waveforms = Array(Vector{UInt8}, 0)
    end

    # Data packets
    try
        while !eof(io)
            timestamp = read(io, UInt32)
            packet_id = read(io, UInt16)
            if packet_id == 0
                # Digital event
                packet_insertion_reason = read(io, UInt8)
                skip(io, 1)
                data = read(io, UInt16)
                if (packet_insertion_reason & UInt8(1)) != 0
                    if (packet_insertion_reason & (UInt8(1) << 7)) != 0
                        push!(serial_digital_channel.times, timestamp/time_resolution_of_time_stamps)
                        push!(serial_digital_channel.data, data)
                    else
                        push!(parallel_digital_channel.times, timestamp/time_resolution_of_time_stamps)
                        push!(parallel_digital_channel.data, data)
                    end
                end
                skip(io, bytes_in_data_packets-10)
            elseif packet_id < 2048
                # Spike event
                channel = spike_channels[packet_id]
                push!(channel.times, timestamp/time_resolution_of_time_stamps)
                unit_number = read(io, UInt8)
                push!(channel.unit_numbers, unit_number)
                skip(io, 1)

                if waveforms
                    append!(electrode_waveforms[spike_channels.number2idx[packet_id]],
                            read!(io, waveform_tmp))
                else
                    skip(io, waveform_bytes)
                end
            elseif packet_id == 0xFFFFFFFF
                error("packet continuation not supported")
            else
                skip(io, bytes_in_data_packets-6)
            end
        end
    catch err
        if isa(err, EOFError)
            warn("encountered EOF; assuming valid anyway (NSP may have crashed)")
        end
    end

    if waveforms
        i = 0
        for x in spike_channels
            x.data = reinterpret(Int16, electrode_waveforms[i += 1],
                                 (waveform_samples, div(length(electrode_waveforms[i]), 2*waveform_samples)))
        end
    end

    digital_channels = nev.digital_channels
    if serial_digital_channel.label != "" || !isempty(serial_digital_channel.times)
        digital_channels[1] = serial_digital_channel
    end
    if parallel_digital_channel.label != "" || !isempty(parallel_digital_channel.times)
        digital_channels[2] = parallel_digital_channel
    end

    nev
end

@struct NSXHeader begin
    file_spec::UInt16
    bytes_in_headers::UInt32
    label::String(16)
    comment::String(256)
    period::UInt32
    time_resolution_of_time_stamps::UInt32
    time_origin::UInt128
end

type NSXContinuousChannels{C,R} <: AbstractChannels{C}
    channels::Vector{C}
    number2idx::Vector{Int}
    idx2number::Vector{Int}
    idx2segidx::R
    segment_data::Vector{Matrix{Int16}}
    segment_cumsum::Vector{Int}
    times::PiecewiseIncreasingRange{Float64,StepRange{Int64,Int64},Int64}

    NSXContinuousChannels(segmap) = new(C[], Int[], Int[], segmap, Matrix{Int16}[], Int[],
        PiecewiseIncreasingRange{Float64,StepRange{Int64,Int64},Int64}(StepRange{Int64,Int64}[],1))
end

type NSXContinuousChannel <: ContinuousChannel
    cc::NSXContinuousChannels
    number::Int

    label::String

    physical_connector::Char
    connector_pin::UInt8
    
    min_digital_value::Int16
    max_digital_value::Int16
    min_analog_value::Int16
    max_analog_value::Int16
    units::String

    high_pass::FilterDescription
    low_pass::FilterDescription

    bits_per_sample::UInt8
    voltage_multiplier::Float64

    NSXContinuousChannel(cc, number) = new(cc, number)
end

NSXContinuousChannels() = NSXContinuousChannels{NSXContinuousChannel,Void}(nothing)

Common.times(c::NSXContinuousChannels) = c.times

function Common.data(c::NSXContinuousChannels{NSXContinuousChannel,Void})
    if length(c.segment_data) == 1
        return c.segment_data[1]
    else
        return hcat(c.segment_data...)
    end
end

function Common.data(c::NSXContinuousChannels)
    if length(c.segment_data) == 1
        return c.segment_data[1][c.idx2segidx, :]
    else
        return hcat([sub(seg, c.idx2segidx, :) for seg in c.segment_data]...)
    end
end

function Common.data(c::NSXContinuousChannels,
                     channels::Union{Int,AbstractVector{Int}},
                     samples::UnitRange{Int}=1:c.segment_cumsum[end])
    first(samples) > 0 && last(samples) < c.segment_cumsum[end] ||
        throw(ArgumentError("samples must be in range [1,$(length(c.chtimes))]"))
    chidx = c.number2idx[channels]
    for x in chidx
        x == 0 && throw(ArgumentError("channel $x is not a valid channel"))
    end
    if !isa(c.idx2segidx, Void)
        chidx = c.idx2segidx[chidx]
    end
    @assert all(0 .< chidx .<= size(c.segment_data[1], 1))

    out = Array(Int16, length(chidx), length(samples))

    # First segment
    istartsegment = searchsortedfirst(c.segment_cumsum, first(samples))
    segment = c.segment_data[istartsegment]
    segmentoffset = istartsegment == 1 ? 0 : c.segment_cumsum[istartsegment]
    startsample = first(samples) - segmentoffset
    endsample = min(last(samples) - segmentoffset, size(segment, 2))
    @assert 0 < startsample <= size(segment, 2)
    @assert 0 < endsample <= size(segment, 2)
    for j = 1:endsample-startsample+1, i = 1:length(chidx)
        @inbounds out[i, j] = segment[chidx[i], startsample+j-1]
    end
    outoffset = endsample-startsample+1

    # Middle segments
    istopsegment = searchsortedlast(c.segment_cumsum, last(samples))
    for isegment = istartsegment+1:istopsegment-1
        segment = c.segment_data[isegment]
        @assert outoffset + size(segment, 2) <= length(samples)
        for j = 1:size(segment, 2), i = 1:length(chidx)
            @inbounds out[i, outoffset+j] = segment[chidx[i], j]
        end
        outoffset += size(segment, 2)
    end

    # Last segment
    if istopsegment > istartsegment
        segment = c.segment_data[istopsegment]
        @assert length(samples)-outoffset <= size(segment, 2)
        for j = 1:length(samples)-outoffset, i = 1:length(chidx)
            @inbounds out[i, outoffset+j] = segment[chidx[i], j]
        end
    end
    out
end

function Common.datavecs(c::NSXContinuousChannels,
                         channels::Union{Int,AbstractVector{Int}}=validchannels(c))
    chidx = c.number2idx[channels]
    for x in chidx
        x == 0 && throw(ArgumentError("channel $x is not a valid channel"))
    end
    if !isa(c.idx2segidx, Void)
        chidx = c.idx2segidx[chidx]
    end

    datavecs = [Array(Int16, length(c.times)) for i = 1:length(channels)]
    n = 0
    for segment in c.segment_data
        timeblock = div(8192, size(segment, 1))
        for iblockstart = 1:timeblock:size(segment, 2)
            timerg = iblockstart:min(iblockstart+timeblock-1, size(segment, 2))
            for ich = 1:length(chidx)
                datavec = datavecs[ich]
                segch = chidx[ich]
                for itp = timerg
                    datavec[n+itp] = segment[segch, itp]
                end
            end
        end
        n += size(segment, 2)
    end
    @assert n == length(c.times)
    datavecs
end

mapseg(::Void, y) = y
mapseg(x, y) = x[y]
function Base.getindex{C}(c::NSXContinuousChannels{C}, channels::AbstractVector)
    chidx = c.number2idx[channels]
    for x in chidx
        x == 0 && throw(ArgumentError("channel $x is not a valid channel"))
    end
    segmap = mapseg(c.idx2segidx, chidx)
    cc = NSXContinuousChannels{C,typeof(segmap)}(segmap)
    cc.channels = c.channels[chidx]
    cc.number2idx = zeros(Int, maximum(channels))
    cc.idx2number = zeros(Int, maximum(channels))
    cc.segment_data = c.segment_data
    cc.segment_cumsum = c.segment_cumsum
    cc.times = c.times
    for i = 1:length(channels)
        cc.number2idx[channels[i]] = i
        cc.idx2number[i] = channels[i]
    end
    cc
end

Common.data(c::NSXContinuousChannel, samples::UnitRange{Int}=1:length(times(c))) =
    vec(data(c.cc, c.number, samples))

Common.hasdata(c::NSXContinuousChannel) = true

Common.times(c::NSXContinuousChannel) = c.cc.times

type NSXFile
    header::NSXHeader
    continuous_channels::NSXContinuousChannels{NSXContinuousChannel}

    NSXFile(header, continuous_channels) = new(header, continuous_channels)
end

function readnsx_info(io::IO)
    read(io, UInt64) == NEURALCD || error("this does not appear to be an NSX file")
    header = read(io, NSXHeader)
    period::Int = header.period
    channels = NSXContinuousChannels()

    channel_count::Int = read(io, UInt32)
    for i = 1:channel_count
        read(io, UInt16) == CC || error("invalid extended header for continuous data")
        electrode_id::Int = read(io, UInt16)
        channel = channels[electrode_id] = NSXContinuousChannel(channels, electrode_id)

        channel.label = readstring(io, 16)
        channel.physical_connector = 'A'-1+read(io, UInt8)
        channel.connector_pin = read(io, UInt8)
        channel.min_digital_value = read(io, Int16)
        channel.max_digital_value = read(io, Int16)
        channel.min_analog_value = read(io, Int16)
        channel.max_analog_value = read(io, Int16)
        units = channel.units = readstring(io, 16)
        channel.high_pass = read(io, FilterDescription)
        channel.low_pass = read(io, FilterDescription)
        channel.bits_per_sample = ceil(Int, log2(Int(channel.max_digital_value) - channel.min_digital_value))
        voltage_multiplier = channel.max_analog_value/channel.max_digital_value
        if units == "mV"
            voltage_multiplier /= 1000
        elseif units == "uV"
            voltage_multiplier /= 10^6
        else
            warn("unknown units for channel $label: $units")
        end
        channel.voltage_multiplier = voltage_multiplier
    end

    return NSXFile(header, channels)
end

"""
    readnsx(io::IO)

Read a Blackrock Microsystems NSX (continuous data) file.
"""
function readnsx(io::IOStream)
    nsx = readnsx_info(io)
    period::Int = nsx.header.period
    nch = length(nsx.continuous_channels)
    tdiv = Int64(nsx.header.time_resolution_of_time_stamps)
    times = StepRange{Int64,Int64}[]
    data = Matrix{Int16}[]
    ichunk = 1
    while !eof(io)
        byte = read(io, UInt8)
        byte == 0x01 || error("unexpected file contents")
        timestamp = read(io, UInt32)
        while !isempty(times) && timestamp < last(last(times))
            warn(@sprintf(
                 "DISCARDING SEGMENT %d: segment ends at %.3f seconds, but next segment starts at %.3f seconds",
                 ichunk-1, last(last(times))/tdiv, timestamp/tdiv))
            pop!(times)
            pop!(data)
        end
        n_data_points = read(io, UInt32)
        corrupted = n_data_points == 0
        if corrupted
            n_data_points = div(filesize(io) - position(io), nch*sizeof(Int16))
            warn("""
                 There are no data points in chunk $(ichunk). The file is probably corrupted. 
                 Assuming the remaining $(n_data_points/period) seconds belong to this chunk.
                 """)
        end
        push!(data, Mmap.mmap(io, Matrix{Int16}, (nch, Int(n_data_points))))
        push!(times, timestamp:period:timestamp+period*(n_data_points-1))
        corrupted && break
        skip(io, nch*n_data_points*sizeof(Int16))
        ichunk += 1
    end
    nsx.continuous_channels.segment_data = data
    nsx.continuous_channels.segment_cumsum = cumsum(Int[length(x) for x in data])
    nsx.continuous_channels.times = PiecewiseIncreasingRange(times, tdiv)
    nsx
end

"""
   fixnsx(io::IO)

Attempt to fix a Blackrock Microsystems NSX file that was corrupted due to a software crash. 
"""
function fixnsx(io::IO)
    nsx = readnsx_info(io)
    period::Int = nsx.header.period
    nch = nchannels(nsx)
    ichunk = 1
    while !eof(io)
        read(io, UInt8) == 0x01 || error("unexpected file contents")
        timestamp = read(io, UInt32)
        n_data_points = read(io, UInt32)
        println("$(Int(timestamp)) $(Int(n_data_points))")
        if n_data_points == 0
            n_data_points = div(filesize(io) - position(io), nch*sizeof(Int16))
            warn("""
                 There are no data points in chunk $(ichunk). The file is probably corrupted. 
                 Assuming the remaining $(n_data_points/period) seconds belong to this chunk.
                 """)
            truncate(io, position(io)+nch*n_data_points*sizeof(Int16))
            skip(io, -sizeof(UInt32))
            write(io, n_data_points)
            return
        end
        skip(io, nch*n_data_points*sizeof(Int16))
        ichunk += 1
    end
end

end