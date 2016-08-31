module Common
using Compat: String
using PiecewiseIncreasingRanges

export File, Channel, ContinuousChannel, SpikeChannel, EventChannel,
       AbstractChannels, Channels, GenericSpikeChannel, validchannels, data, datavecs,
       times, hasdata, samplerate, bitspersample, voltagemultiplier, unitnumbers, label

abstract File
abstract NeuroChannel
abstract ContinuousChannel <: NeuroChannel
abstract SpikeChannel <: NeuroChannel
abstract EventChannel <: NeuroChannel
abstract AbstractChannels{C<:NeuroChannel}

immutable Unit{W,T}
    number::Int
    times::Vector{Float64}
    data::W
    multiplier::T
end

"""
   Channels{C<:NeuroChannel}()

A container for channels of type C that is iterable but indexed by channel number
"""
immutable Channels{C<:NeuroChannel} <: AbstractChannels{C}
    channels::Vector{C}           # all channels that exist
    number2idx::Vector{Int}       # mapping between channel number and index
    idx2number::Vector{Int}       # mapping between index and channel number
    Channels() = new(C[], Int[], Int[])
end

# Iteration returns only channels that exist
Base.start(c::AbstractChannels) = 1
Base.next(c::AbstractChannels, i) = (c.channels[i], i+1)
Base.done(c::AbstractChannels, i) = i == length(c)+1
Base.eltype{C}(c::AbstractChannels{C}) = C
Base.length(c::AbstractChannels) = length(c.channels)

immutable ChannelEnumerator{C<:AbstractChannels}
    x::C
end
Base.start(c::ChannelEnumerator) = 1
Base.next(c::ChannelEnumerator, i) = ((c.x.idx2number[i], c.x.channels[i]), i+1)
Base.done(c::ChannelEnumerator, i) = i == length(c.x)+1
Base.eltype(c::ChannelEnumerator) = Tuple{Int, eltype(c.x)}
Base.length(c::ChannelEnumerator) = length(c.x)
Base.enumerate(c::AbstractChannels) = ChannelEnumerator(c)

function Base.getindex{C}(c::AbstractChannels{C}, x::AbstractVector)
    newch = Channels{C}()
    for i in x
        newch[i] = c[i]
    end
    newch
end

Base.getindex{C}(c::AbstractChannels{C}, i) = c.channels[c.number2idx[i]]
function Base.setindex!(c::AbstractChannels, v, i)
    n = length(c.number2idx)

    if i > n
        resize!(c.number2idx, i)
        c.number2idx[n+1:i] = 0
    end

    if c.number2idx[i] == 0
        push!(c.channels, v)
        push!(c.idx2number, i)
        c.number2idx[i] = length(c.channels)
    else
        c.channels[c.number2idx[i]] = v
    end
    c
end
Base.endof(c::AbstractChannels) = isempty(c.channels) ? 0 : endof(c.number2idx)
Base.haskey(c::AbstractChannels, i::Integer) =
    0 < i <= length(c.number2idx) && c.number2idx[i] != 0

Base.isassigned(c::AbstractChannels, i) =
    i <= length(c.number2idx) && c.number2idx[i] != 0

function Base.show(io::IO, c::AbstractChannels)
    write(io, "$(length(c))-channel $(typeof(c))")
    for (i, ch) in enumerate(c)
        write(io, "\n  $i => $ch")
    end
end

"""
    validchannels(c::AbstractChannels)

Get the channel numbers for all defined channels
"""
validchannels(c::AbstractChannels) = c.idx2number

"""
    NeuroIO.times(c::Union{Unit,NeuroChannel})

Get the times at which spikes, events, or recorded samples were acquired.
"""
times(c::Union{Unit,NeuroChannel}) = c.times

"""
    NeuroIO.times(c::AbstractChannels{C<:ContinuousChannel})

Get the times corresponding to recorded samples on all continuous channels, if
they had a common time basis.
"""
function times{C<:ContinuousChannel}(c::AbstractChannels{C})
    length(c) == 0 && throw(ArgumentError("no channels defined"))
    times = c.channels[1].times
    for i = 2:length(c.channels)
        if c.channels[i].times != times
            throw(ArgumentError("channels do not have a common time basis"))
        end
    end
    return times
end

"""
    NeuroIO.data(c::Union{Unit,NeuroChannel})

Get the spike waveforms (for `SpikeChannel` or `Unit`), event codes (for
`EventChannel`), or continuous data (for `ContinuousChannel`).
"""
data(c::NeuroChannel) = c.data

"""
    NeuroIO.data(c::ContinuousChannel, samples=1:length(NeuroIO.times(c)))

Get the recorded samples for the given channel. If scaled is false, the raw data
is returned. If scaled is true, samples are returned in microvolts.
"""
data(c::NeuroChannel, samples::AbstractVector{Int}) = sub(data(c), samples)

"""
    NeuroIO.data(c::AbstractChannels{C<:ContinuousChannel},
                 channels=validchannels(c),
                 samples=1:length(chtimes(c)))

Get data from continuous channels recorded with a common time basis as an
nchannels × nsamples matrix. `samples` specifies the sample indexes. `channels`
specifies the specific channels to be returned. If set to `nothing`, all
channels are returned.
"""
function data{C<:ContinuousChannel}(c::AbstractChannels{C},
                                    channels::AbstractVector{Int}=validchannels(c),
                                    samples::Range{Int}=1:length(chtimes(c)))
    length(channels) == 0 && throw(ArgumentError("no channels specified"))

    # Ensure common time basis
    times = c[channels[1]].times
    for i = 2:length(channels)
        if c[channels[i]].times != times
            throw(ArgumentError("channels do not have a common time basis"))
        end
    end

    data = Array(eltype(c[1]), length(samples), length(channels))
    for i = 1:length(channels)
        copy!(data, size(data, 1)*(i-1)+1,
              sub(chdata(c[channels[i]]), samples))
    end
    data.'
end

"""
    NeuroIO.datavecs(c::AbstractChannels, channels=nothing)

Get data from multiple channels simultaneously and return a vector of vectors.
This can be more efficient than calling NeuroIO.data on each channel,
depending on the on-disk layout.
"""
datavecs(c::AbstractChannels, channels::AbstractVector{Int}=validchannels(c)) =
    [data(c[x]) for x in channels]

"""
    NeuroIO.hasdata(c::NeuroChannel)

Whether data (waveforms, event codes, or samples) is available for this channel.
"""
hasdata(c::NeuroChannel) = isdefined(c, :data)

"""
    samplerate(c::SpikeChannel)

Get the sample rate for the given spike channel.
"""
samplerate(c::SpikeChannel) = c.sample_rate

"""
    samplerate(c::ContinuousChannel)

Get the sample rate for the given continuous channel.
"""
samplerate(c::ContinuousChannel) = _samplerate(times(c))
_samplerate{T,R<:UnitRange}(r::PiecewiseIncreasingRange{T,R}) = r.divisor
_samplerate(r::PiecewiseIncreasingRange) = r.divisor//step(r.ranges[1])

"""
    bitspersample(c::Union{SpikeChannel,ContinuousChannel})

Get the number of bits encoded for a waveform or continuous sample.
"""
bitspersample(c::Union{SpikeChannel,ContinuousChannel}) = c.bits_per_sample

"""
    voltagemultiplier(c::Union{SpikeChannel,ContinuousChannel})

Get multiplier to convert raw AD values to volts, so that
`NeuroIO.data(c)*voltagemultiplier(c)` gives the recored data voltages.
"""
voltagemultiplier(c::Union{SpikeChannel,ContinuousChannel}) = c.voltage_multiplier

"""
    NeuroIO.label(c::NeuroChannel)

Get the channel label, if one exists.
"""
label(c::NeuroChannel) = c.label

"""
    unitnumbers(c::SpikeChannel)

Unit numbers corresponding to spikes on the given channel. 0 if unsorted.
"""
unitnumbers(c::SpikeChannel) = c.unit_numbers

immutable GenericSpikeChannel{T,S} <: SpikeChannel
    label::String
    sample_rate::S
    voltage_multiplier::Float64
    times::Vector{Float64}
    unit_numbers::Vector{Int16}
    data::Matrix{T}
end
bitspersample{T}(::GenericSpikeChannel{T}) = sizeof(T)*8

Base.show(io::IO, c::SpikeChannel) =
    @printf io "%s %s %d spikes" typeof(c) label(c) length(times(c))
Base.show(io::IO, c::EventChannel) =
    @printf io "%s %s %d spikes" typeof(c) label(c) length(times(c))
Base.show(io::IO, c::ContinuousChannel) =
    @printf io "%s %s %d samples @ %.1f Hz" typeof(c) label(c) length(times(c)) samplerate(c)
end