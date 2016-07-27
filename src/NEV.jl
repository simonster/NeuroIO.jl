module NEV
using ..Struct, PiecewiseIncreasingRanges
export NEVFile, NSXFile, readnev, readnsx, readnsx_mmap, nchannels, channelindexes

@struct NEVHeader begin
    file_spec::UInt16
    additional_flags::UInt16
    bytes_in_headers::UInt32
    bytes_in_data_packets::UInt32
    time_resolution_of_time_stamps::UInt32
    time_resolution_of_samples::UInt32
    time_origin::UInt128
    application_to_create_file::ASCIIString(32)
    comment_field::UTF8String(256)
end

@struct FilterDescription begin
    cutoff::UInt32
    order::UInt32
    filter_type::UInt16
end

type SpikeChannel
    label::UTF8String

    times::Vector{Float64}
    unit_numbers::Vector{UInt8}
    waveforms::Matrix{Int16}

    physical_connector::Char
    connector_pin::UInt8

    digitization_factor::UInt16
    energy_threshold::UInt16
    high_threshold::Int16
    low_threshold::Int16
    number_of_sorted_units::UInt8
    bytes_per_waveform::UInt8
    spike_width::UInt16

    high_pass::FilterDescription
    low_pass::FilterDescription

    SpikeChannel() = new("", Float64[], UInt8[])
end

type DigitalChannel
    label::UTF8String
    times::Vector{Float64}
    data::Vector{UInt16}

    DigitalChannel() = new("", Float64[], UInt16[])
end

type NEVFile
    header::NEVHeader
    digital_channels::Vector{DigitalChannel}
    spike_channels::Vector{SpikeChannel}
    array_name::UTF8String
    map_file::UTF8String

    NEVFile(header) = new(header, DigitalChannel[], SpikeChannel[])
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
    bytestring(data[1:findfirst(data, '\0')-1])
end

function readnev(io::IO; waveforms::Bool=true)
    read(io, UInt64) == NEURALEV || error("this does not appear to be an NEV file")
    header = read(io, NEVHeader)
    bytes_in_data_packets = header.bytes_in_data_packets
    time_resolution_of_time_stamps = header.time_resolution_of_time_stamps
    nev = NEVFile(header)
    spike_channels = nev.spike_channels
    serial_digital_channel = DigitalChannel()
    parallel_digital_channel = DigitalChannel()

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
            electrode_id > length(spike_channels) && (resize!(spike_channels, electrode_id))
            channel = isdefined(spike_channels, electrode_id) ? spike_channels[electrode_id] : (spike_channels[electrode_id] = SpikeChannel())

            if packet_id == NEUEVWAV
                channel.physical_connector = 'A'-1+read(io, UInt8)
                channel.connector_pin = read(io, UInt8)
                channel.digitization_factor = read(io, UInt16)
                channel.energy_threshold = read(io, UInt16)
                channel.high_threshold = read(io, Int16)
                channel.low_threshold = read(io, Int16)
                channel.number_of_sorted_units = read(io, UInt8)
                channel.bytes_per_waveform = read(io, UInt8)
                channel.spike_width = read(io, UInt16)
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
            label::UTF8String = readstring(io, 16)
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

    electrode_waveforms = Array(Vector{UInt8}, waveforms ? length(spike_channels) : 0)
    electrode_waveform_count = zeros(Int, length(electrode_waveforms))
    waveform_bytes = bytes_in_data_packets-8
    waveform_samples = div(waveform_bytes, 2)
    waveform_tmp = Array(UInt8, waveform_bytes)
    if waveforms
        for i = 1:length(spike_channels)
            if isdefined(spike_channels, i)
                electrode_waveforms[i] = UInt8[]
            end
        end
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
                    append!(electrode_waveforms[packet_id], read!(io, waveform_tmp))
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
        for i = 1:length(spike_channels)
            if isdefined(spike_channels, i)
                spike_channels[i].waveforms = reinterpret(Int16, electrode_waveforms[i], (waveform_samples, div(length(electrode_waveforms[i]), 2*waveform_samples)))
            end
        end
    end

    digital_channels = nev.digital_channels
    if serial_digital_channel.label != "" || !isempty(serial_digital_channel.times)
        push!(digital_channels, serial_digital_channel)
    end
    if parallel_digital_channel.label != "" || !isempty(parallel_digital_channel.times)
        push!(digital_channels, parallel_digital_channel)
    end

    nev
end

@struct NSXHeader begin
    file_spec::UInt16
    bytes_in_headers::UInt32
    label::ASCIIString(16)
    comment::ASCIIString(256)
    period::UInt32
    time_resolution_of_time_stamps::UInt32
    time_origin::UInt128
end

type ContinuousChannel
    samples::Vector{Int16}
    times::PiecewiseIncreasingRange{Float64,StepRange{Int64,Int64},Int64}

    label::UTF8String

    physical_connector::Char
    connector_pin::UInt8
    
    min_digital_value::Int16
    max_digital_value::Int16
    min_analog_value::Int16
    max_analog_value::Int16
    units::UTF8String

    high_pass::FilterDescription
    low_pass::FilterDescription

    ContinuousChannel() = new(Int16[])
end

type NSXFile
    header::NSXHeader
    continuous_channels::Vector{ContinuousChannel}
    times::PiecewiseIncreasingRange{Float64,StepRange{Int64,Int64},Int64}
    mmapped_data::Vector{Matrix{Int16}}

    NSXFile(header, continuous_channels) = new(header, continuous_channels)
end

function readnsx_info(io::IO)
    read(io, UInt64) == NEURALCD || error("this does not appear to be an NSX file")
    header = read(io, NSXHeader)
    period::Int = header.period
    channels = ContinuousChannel[]
    times = StepRange{Int64,Int64}[]

    channel_count::Int = read(io, UInt32)
    for i = 1:channel_count
        read(io, UInt16) == CC || error("invalid extended header for continuous data")
        electrode_id::Int = read(io, UInt16)
        electrode_id > length(channels) && (resize!(channels, electrode_id))
        channel = channels[electrode_id] = ContinuousChannel()

        channel.label = readstring(io, 16)
        channel.physical_connector = 'A'-1+read(io, UInt8)
        channel.connector_pin = read(io, UInt8)
        channel.min_digital_value = read(io, Int16)
        channel.max_digital_value = read(io, Int16)
        channel.min_analog_value = read(io, Int16)
        channel.max_analog_value = read(io, Int16)
        channel.units = readstring(io, 16)
        channel.high_pass = read(io, FilterDescription)
        channel.low_pass = read(io, FilterDescription)
    end

    return NSXFile(header, channels)
end

function readnsx(io::IO)
    nsx = readnsx_info(io)
    period::Int = nsx.header.period
    channels = nsx.channels
    times = StepRange{Int64,Int64}[]

    channels_by_index = ContinuousChannel[]
    estimated_samples = isa(io, IOStream) ? div(filesize(io), channel_count*2) : 0
    for i = 1:length(channels)
        if isdefined(channels, i)
            push!(channels_by_index, channels[i])

            if estimated_samples != 0
                sizehint!(channels[i].samples, estimated_samples)
            end
        end
    end
    @assert length(channels_by_index) == channel_count

    n_total_data_points = 0

    # To amortize overhead, read approximately 16 KB at a time
    large_data_tmp = zeros(Int16, channel_count, div(2^13, channel_count))
    small_data_tmp = zeros(Int16, channel_count)

    # Otherwise, read! calls reinterpret internally, which adds overhead
    # TODO PR to change this
    large_data_tmp_r = reinterpret(UInt8, large_data_tmp, (2*channel_count*size(large_data_tmp, 2),))
    small_data_tmp_r = reinterpret(UInt8, small_data_tmp)

    while !eof(io)
        read(io, UInt8) == 0x01 || error("unexpected file contents")
        timestamp = read(io, UInt32)
        n_data_points = read(io, UInt32)

        # Resize channels to new size
        n_old_data_points = n_total_data_points
        n_total_data_points += n_data_points
        for channel in channels_by_index
            resize!(channel.samples, n_total_data_points)
        end

        # First read data in bigger chunks for performance
        lrg = 1:size(large_data_tmp, 2):n_data_points-size(large_data_tmp, 2)
        for k = lrg
            read!(io, large_data_tmp_r)
            @inbounds for j = 1:length(channels_by_index)
                samples = channels_by_index[j].samples
                for i = 1:size(large_data_tmp, 2)
                    samples[n_old_data_points+k+i-1] = large_data_tmp[j, i]
                end
            end
        end

        # Then read remaining data
        srg = last(lrg)+size(large_data_tmp, 2):n_data_points
        for k = srg
            read!(io, small_data_tmp_r)
            @inbounds for j = 1:length(channels_by_index)
                channels_by_index[j].samples[n_old_data_points+k] = small_data_tmp[j]
            end
        end

        # Add times
        push!(times, timestamp:period:timestamp+period*(n_data_points-1))
    end

    times = PiecewiseIncreasingRange(times, Int64(nsx.header.time_resolution_of_time_stamps))
    for ch in channels_by_index
        ch.times = times
    end
    nsx.times = times
    nsx
end

function readnsx_mmap(io::IO)
    nsx = readnsx_info(io)
    period::Int = nsx.header.period
    nch = nchannels(nsx)
    times = StepRange{Int64,Int64}[]
    data = Matrix{Int16}[]
    ichunk = 1
    while !eof(io)
        byte = read(io, UInt8)
        byte == 0x01 || error("unexpected file contents")
        timestamp = read(io, UInt32)
        n_data_points = read(io, UInt32)
        corrupted = n_data_points == 0
        if corrupted
            n_data_points = div(filesize(io) - position(io), nch*sizeof(Int16))
            warn("There are no data points in chunk $(ichunk). The file is probably corrupted. Assuming the remaining $(n_data_points/period) seconds belong to this chunk.")
        end
        push!(data, Mmap.mmap(io, Matrix{Int16}, (nch, n_data_points)))
        push!(times, timestamp:period:timestamp+period*(n_data_points-1))
        corrupted && break
        skip(io, nch*n_data_points*sizeof(Int16))
        ichunk += 1
    end
    nsx.mmapped_data = data
    nsx.times = PiecewiseIncreasingRange(times, Int64(nsx.header.time_resolution_of_time_stamps))
    nsx
end

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
            warn("There are no data points in chunk $(ichunk). The file is probably corrupted. Assuming the remaining $(n_data_points/period) seconds belong to this chunk.")
            truncate(io, position(io)+nch*n_data_points*sizeof(Int16))
            skip(io, -sizeof(UInt32))
            write(io, n_data_points)
            return
        end
        skip(io, nch*n_data_points*sizeof(Int16))
        ichunk += 1
    end
end

nchannels(f::NEVFile) = ndefined(f.spike_channels)
nchannels(f::NSXFile) = ndefined(f.continuous_channels)
function ndefined(arr)
    n = 0
    for i = 1:length(arr)
        if isdefined(arr, i)
            n += 1
        end
    end
    n
end

channelindexes(f::NEVFile) = validindexes(f.spike_channels)
channelindexes(f::NSXFile) = validindexes(f.continuous_channels)
function validindexes(arr)
    inds = Int[]
    for i = 1:length(arr)
        if isdefined(arr, i)
            push!(inds, i)
        end
    end
    inds
end

end