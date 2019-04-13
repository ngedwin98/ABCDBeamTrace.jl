# plot recipes, see https://docs.juliaplots.org/latest/recipes/ and
# https://github.com/JuliaPlots/RecipesBase.jl

mutable struct WithGaussianBeam
    system::Vector{<:Element}
    beam::GaussianBeam
    unit::Number
end

# type recipe, e.g. for `plot(WithGaussianBeam(system, beam))`
@recipe f(::Type{WithGaussianBeam}, data::WithGaussianBeam) = begin
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
    unit = 1data.unit
    # note: the following unit conversion could be achieved without
    #       making Unitful a dependency, by adding 0.0 and the desired
    #       fraction
    ws1 = Unitful.uconvert(
        Unitful.NoUnits,
        spotradius(data.beam) / unit
    )
    zs1 = Unitful.uconvert(
        Unitful.NoUnits,
        location(data.beam) / unit
    )
    ws = Vector{typeof(ws1)}(undef, N)
    zs = Vector{typeof(zs1)}(undef, N)
    ws[1] = ws1
    zs[1] = zs1
    beam = data.beam
    for i = 1:length(ds)
        beam = transform(ds[i], beam)
        ws[i+1] = Unitful.uconvert(
            Unitful.NoUnits,
            spotradius(beam) / unit
        )
        zs[i+1] = Unitful.uconvert(
            Unitful.NoUnits,
            location(beam) / unit
        )
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
wavelength (specified either as a `Unitful.Quantity` or in SI units,
i.e. if dimensionless, then in the range from approximately `380e-9`
to `700e-9`).

"""
function color(λ::AbstractFloat)
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
color(λ::Unitful.Length) =
    color(Unitful.uconvert(Unitful.NoUnits, float(λ) / Unitful.m))

# user recipe, e.g. for `plot(system, beam)`
@recipe f(system::Vector{<:Element}, beam::GaussianBeam) =
    WithGaussianBeam(
        system,
        beam,
        (
            Unitful.dimension(location(beam)) ==
            Unitful.dimension(Unitful.m)
        ) ? 1Unitful.m : 1
    )
