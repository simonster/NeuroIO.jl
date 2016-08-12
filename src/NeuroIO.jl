module NeuroIO
using Reexport

include("Struct.jl")
include("Common.jl")
include("NCS.jl")
include("PLX.jl")
include("NEV.jl")
include("Klusters.jl")

@reexport using .NEV
@reexport using .PLX
@reexport using .Klusters
@reexport using .NCS
using .Common
export validchannels, samplerate, bitspersample, voltagemultiplier, unitnumbers

end # module
