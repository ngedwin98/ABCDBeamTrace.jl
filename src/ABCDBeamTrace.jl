module ABCDBeamTrace

export FreeSpace, Interface, ThinLens, Mirror, Tan, Sag, RTM, Beam
export beamtrace, spotsize, location, discretize
export waistradiusfunc
export WithBeam

using RecipesBase
import Colors, Interpolations

include("elements.jl")
include("beamtrace.jl")
include("comparisons.jl")
include("plot-recipes.jl")

end # module
