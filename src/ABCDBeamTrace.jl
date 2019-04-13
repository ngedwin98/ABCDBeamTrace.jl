module ABCDBeamTrace

using RecipesBase

export FreeSpace, Interface, ThinLens, Mirror, Tan, Sag, RTM, Beam
export beamtrace, spotsize, location, waistradiusfunc, discretize

include("elements.jl")
include("beamtrace.jl")
include("comparisons.jl")
include("plot-recipes.jl")

end # module
