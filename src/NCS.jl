module NCS
using ..Common, PiecewiseIncreasingRanges
import Compat.String
export readncs

immutable NCSContinuousChannel <: ContinuousChannel
    header::String
    data::Vector{Int16}
    times::PiecewiseIncreasingRange{Float64,StepRange{Int64,Int64},Int64}
    voltage_multiplier::Float64
end

immutable NCSRecordHeader
    qwTimeStamp::UInt64
    dwChannelNumber::UInt32
    dwSampleFreq::UInt32
    dwNumValidSamples::UInt32
end

Base.read(io::IO, ::Type{NCSRecordHeader}) =
    NCSRecordHeader(read(io, UInt64), read(io, UInt32), read(io, UInt32), read(io, UInt32))

function compute_times(rec::NCSRecordHeader, sample_multiplier::Int64, step::Int64)
    first_sample = Int64(rec.qwTimeStamp)*sample_multiplier
    first_sample:step:first_sample+step*(rec.dwNumValidSamples-1)
end

function readncs(filename::AbstractString; kwargs...)
    io = open(filename, "r")
    try
        readncs(io, filesize(io); kwargs...)
    finally
        close(io)
    end
end

function readncs(io::IO, sz::Union{Integer,Void}=nothing; skip_nonmonotonic::Bool=false)
    header = rstrip(bytestring(read(io, UInt8, 16384)), '\0')
    eof(io) && return NCSContinuousChannel(header, Int16[], PiecewiseIncreasingRange(StepRange{Int,Int}[], 1))

    nrecs::Int = isa(sz, Integer) ? div(sz-16384, 1044) : -1
    sample_buffer = Array(Int16, 512)
    sample_buffer_reinterp = reinterpret(UInt8, sample_buffer)
    samples = Array(Int16, isa(sz, Integer) ? nrecs*512 : 512)
    nsamples = 0
    times = Array(StepRange{Int,Int}, isa(sz, Integer) ? nrecs : 1)

    # Read first header so we can save its channel number to make
    # sure we only have one
    rec1 = read(io, NCSRecordHeader)
    read!(io, sample_buffer)
    copy!(samples, nsamples+1, sample_buffer, 1, rec1.dwNumValidSamples)
    nsamples += rec1.dwNumValidSamples

    # Figure out range multipliers

    # This would be correct if the reported sample rate were:
    # @compat divisor = lcm(Int64(rec1.dwSampleFreq), 10^6)
    # sample_multiplier = div(divisor, 10^6)
    # @compat step = div(divisor, Int64(rec1.dwSampleFreq))

    # However, it seems that the reported sample rate is wrong. The
    # actual sample rate is rounded so that the delta between
    # timestamps is constant.
    sample_multiplier = 512
    divisor = 512*10^6
    step = div(divisor, rec1.dwSampleFreq)

    nrecs = 1
    reccount = 1
    times[1] = compute_times(rec1, sample_multiplier, step)

    while !eof(io)
        reccount += 1
        rec = read(io, NCSRecordHeader)
        rec.dwChannelNumber == rec1.dwChannelNumber || error("only one channel supported per file")
        rec.dwSampleFreq == rec1.dwSampleFreq || error("sample rate is non-constant")
        read!(io, sample_buffer_reinterp)
        rectimes = compute_times(rec, sample_multiplier, step)
        if first(rectimes) < last(times[nrecs])
            if skip_nonmonotonic
                warn("""timestamp at record $(reccount) is not monotonically increasing
                        (previous end time = $(last(times[nrecs])/divisor), start time = $(first(rectimes)/divisor); skipping""")
                continue
            else
                error("""timestamp at record $(reccount) is not monotonically increasing.
                         Use readncs(..., skip_nonmotonoic=true) to skip nonmonotonic timestamps""")
            end
        end
        isa(sz, Void) && resize!(samples, nsamples+rec.dwNumValidSamples)
        copy!(samples, nsamples+1, sample_buffer, 1, rec.dwNumValidSamples)
        nsamples += rec.dwNumValidSamples
        nrecs += 1
        if isa(sz, Integer)
            times[nrecs] = rectimes
        else
            push!(times, rectimes)
        end
    end

    resize!(samples, nsamples)
    resize!(times, nrecs)

    voltage_multiplier = NaN
    m = match(r"\n\s*-ADBitVolts\s+([0-9.]+)", header)
    if m !== nothing
        voltage_multiplier = parse(Float64, m.captures[1])
    end

    return NCSContinuousChannel(header, samples, PiecewiseIncreasingRange(times, divisor), voltage_multiplier)
end

end

