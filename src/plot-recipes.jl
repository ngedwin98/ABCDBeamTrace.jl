# plot recipes, see https://docs.juliaplots.org/latest/recipes/ and
# https://github.com/JuliaPlots/RecipesBase.jl

@recipe function f(system::Vector{<:Element})
    return waistradiusfunc(system, beam)
end
