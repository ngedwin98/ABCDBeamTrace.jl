# plot recipes, see https://docs.juliaplots.org/latest/recipes/ and
# https://github.com/JuliaPlots/RecipesBase.jl

mutable struct WithBeam
    system::Vector{<:Element}
    beam::Beam
end

# type recipe, e.g. for `plot(WithBeam(system, beam))`
@recipe f(::Type{WithBeam}, data::WithBeam) = begin
    seriestype := :shape
    linecolor --> color(data.beam.λ)
    fillcolor --> color(data.beam.λ)
    fillalpha --> 0.1
    label --> string(data.beam.λ)
    aspect_ratio --> :equal
    xlabel --> "Distance along Beam Axis"
    ylabel --> "1/e^2 Boundary"
    # choose a very large number of discretizations to have some
    # chance of approximating minimum waist radii decently
    ds = discretize(data.system, 200)
    N = length(ds) + 1
    Tw = typeof(float(spotsize(data.beam)))
    Tz = typeof(float(location(data.beam)))
    ws = Vector{Tw}(undef, N)
    zs = Vector{Tz}(undef, N)
    ws[1] = spotsize(data.beam)
    zs[1] = location(data.beam)
    beam = data.beam
    for i = 1:length(ds)
        beam = transform(ds[i], beam)
        ws[i+1] = spotsize(beam)
        zs[i+1] = location(beam)
    end
    xs, ys = vcat(zs, reverse(zs)), vcat(ws, (-1.0) .* reverse(ws))
    [(xs[i], ys[i]) for i in 1:length(xs)]
end

const _λ_hue_data = [
    380e-9 260 # "far" UV
    440e-9 250 # violet
    485e-9 180 # blue
    523e-9 120 # blue-green
    572e-9 60 # yellow
    620e-9 20
    700e-9 0 # (infra) red
]
const hue_of_λ = Interpolations.LinearInterpolation(
    _λ_hue_data[:,1], _λ_hue_data[:,2]
)

"""

Return a (very approximately) correct color from the package
[Colors](http://docs.juliaplots.org/latest/colors/) for a given
wavelength (specified in SI units, i.e. in the range from
approximately `380e-9` to `700e-9`).

"""
function color(λ)
    # very rough transformation of λ to a color via package Colors
    λmin = 380e-9
    λmax = 700e-9
    if λ < λmin
        # ultraviolet (far spectral violet, dark)
        return Colors.HSL(260.0, 1.0, 0.25)
    elseif λ > λmax
        # infrared (red, dark)
        return Colors.HSL(0, 1.0, 0.25)
    end
    color = Colors.HSL(hue_of_λ(λ), 1.0, 0.5)
end

# user recipe, e.g. for `plot(system, beam)`
@recipe f(system::Vector{<:Element}, beam::Beam) =
    WithBeam(system, beam)
