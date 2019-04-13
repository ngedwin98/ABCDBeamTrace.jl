"""

An element is the abstract supertype of all elements that together
form an optical system, consistently modeled as a `Vector{<:Element}`.

$(SIGNATURES)

"""
abstract type Element end
struct FreeSpace <: Element
    L::Real
end

"""

An optical Interface.

$(SIGNATURES)

The Interface has the parameters `η = n_next / n_prev` where `n` is
the optical density (aka optical index), the angle of incidence `Θ`,
and a radius of curvature `R` which can be infinite. All but the first
arguments are optional; the defaults are zero angle of incidence and
infinite radius of curvature.

"""
struct Interface <: Element
    η::Real; θ::Real; R::Real # R > 0 when light hits concave side
end
Interface(η::Real, θ::Real = 0) = Interface(η, θ, Inf)

"""

A thin lens.

$(SIGNATURES)

The parameters are focal length `f` and an optical angle of incidence
`Θ`.

"""
struct ThinLens <: Element
    f::Real; θ::Real
end
ThinLens(f::Real) = ThinLens(f,0)

"""

A mirror.

$(SIGNATURES)

The arguments are radius of curvature `R` and angle of incidence `θ`.

!!! note "To Do"
    The order of these argument is reversed compared to
    [`Interface`](@ref). Unify the order, or switch to keyword
    arguments.

"""
Mirror(R::Real, θ::Real) = ThinLens(R/2,θ)

"""

Construct a pseudo-element for modeling beam propagation in the
tangential (aka parallel) plane.

$(SIGNATURES)

"""
struct Tan{E} <: Element where {E<:Element}
    e::E
end
Tan(e::FreeSpace) = e

"""

Construct a pseudo-element for modeling beam propagation in the
sagittal plane.

"""
struct Sag{E} <: Element where {E<:Element}
    e::E
end
Sag(e::FreeSpace) = e

"""

RTM(element::Element) returns the Ray Transfer (ABCD) matrix
associated with the given, optical element.

RTM(elements) returns the Ray Transfer (ABCD) matrix associated with
an optical system described by a collection (e.g. a vector or
iteration) of optical elements.

"""
RTM(e::FreeSpace) = [1 e.L ; 0 1]
RTM(e::Interface) = [1 0 ; (e.η-1)/e.R e.η]
RTM(e::ThinLens) = [1 0 ; -1/e.f 1]
RTM(e::Tan{ThinLens}) = (RTM∘ThinLens)(e.e.f*cos(e.e.θ))
RTM(e::Sag{ThinLens}) = (RTM∘ThinLens)(e.e.f/cos(e.e.θ))
# See doi:10.1364/AO.26.000427 for the following matrices
function RTM(e::Tan{Interface})
    θ1, η, R = e.e.θ, e.e.η, e.e.R; θ2 = asin(η*sin(θ1))
    return [cos(θ2)/cos(θ1) 0 ;
        (cos(θ2)-η*cos(θ1))/(R*cos(θ1)*cos(θ2)) η*cos(θ1)/cos(θ2)]
end
function RTM(e::Sag{Interface})
    θ1, η, R = e.e.θ, e.e.η, e.e.R; θ2 = asin(η*sin(θ1))
    return [1 0 ; (cos(θ2)-η*cos(θ1))/R η]
end
RTM(elements) = reduce(*, RTM.(elements))
