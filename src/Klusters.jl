module Klusters
using ..Common, LightXML
export readklusters

"""
readklusters(path::AbstractString) -> GenericSpikeChannel

Read spike data from Klusters.
"""
function readklusters(path::AbstractString)
    dotindex = rsearch(path, '.')
    stem = path[1:rsearch(path, '.', dotindex-1)-1]
    chnumber = parse(Int, path[dotindex+1:end])
    
    xdoc = parse_file(stem*".xml")
    xroot = root(xdoc)
    acquisitionSystem = get_elements_by_tagname(xroot, "acquisitionSystem")[1]
    sample_rate = parse(Float64, content(get_elements_by_tagname(acquisitionSystem, "samplingRate")[1]))
    nbits = parse(Int, content(get_elements_by_tagname(acquisitionSystem, "nBits")[1]))
    voltage_multiplier = parse(Float64, content(get_elements_by_tagname(acquisitionSystem, "voltageRange")[1]))/2^nbits
    nbits == 16 || nbits == 32 || throw(ArgumentError("nBits == $nbits unsupported"))

    spikeDetection = get_elements_by_tagname(xroot, "spikeDetection")[1]
    groups = get_elements_by_tagname(get_elements_by_tagname(spikeDetection, "channelGroups")[1], "group")
    nsamples = 0
    for group in groups
        has_channel = false
        for channel in get_elements_by_tagname(get_elements_by_tagname(group, "channels")[1], "channel")
            if parse(Int, content(channel)) == chnumber
                has_channel = true
                break
            end
        end
        !has_channel && continue
        nsamples = parse(Int, content(get_elements_by_tagname(group, "nSamples")[1]))
    end
    nsamples == 0 && throw(ArgumentError("number of samples for channel $chnumber not specified in $stem.xml"))

    times = open("$stem.fet.$chnumber") do f
        readline(f)
        out = Float64[]
        for line in eachline(f)
            line == "\n" || line == "\r\n" && continue
            push!(out, parse(Int, line[rsearch(line, ' ')+1:end]))
        end
        out
    end
    times /= sample_rate
    if isfile("$stem.clu.$chnumber")
        unit_numbers = open("$stem.clu.$chnumber") do f
            readline(f)
            out = Int16[]
            for line in eachline(f)
                line == "\n" || line == "\r\n" && continue
                push!(out, parse(Int, line))
            end
            out
        end
    else
        unit_numbers = zeros(Int16, length(times))
    end
    data = open(f->read(f, nbits == 16 ? Int16 : Int32, (nsamples, div(filesize(f), nsamples*nbits÷8))), "$stem.spk.$chnumber")

    length(times) == length(unit_numbers) == size(data, 2) ||
        throw(ArgumentError("number of waveforms, times, and unit numbers do not match"))

    GenericSpikeChannel("ch$chnumber", sample_rate, voltage_multiplier, times, unit_numbers, data)
end
end
