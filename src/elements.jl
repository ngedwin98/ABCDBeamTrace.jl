"""

An element is the abstract supertype of all elements that together
form an optical system, consistently modeled as a `Vector{<:Element}`.

$(SIGNATURES)

It is parametrized on types L and N, used for numbers of the
respective dimenions length (L) and dimensionless (N, "number").

"""
abstract type Element{L,N} end

"""

A generic, optical element constructed by giving its ray transfer
matrix (or the ABCD entries, in order).

$(SIGNATURES)

"""
struct ElementABCD{L,N} <: Element{L,N}
    A::N
    B::L
    invC::L
    D::N
    # an inner constructor is needed to avoid that field invC gets
    # assigned directly (with C, not inv(C) which it needs to be
    # assigned)
    function ElementABCD(A::N, B::L, C, D::N) where {L,N}
        if zero(one(N)) != zero(N)
            throw(DomainError("A and D must be dimensionless"))
        end
        if zero(one(B * C)) != zero(B * C)
            throw(DomainError("B and C must have inverse dimensions"))
        end
        L2 = typeof(1.0B) # promote to float even if float(B) is not
                          # defined
        N2 = typeof(1.0A) # promote to float even if ...
        return new{L2,N2}(A, B, inv(C), D)
    end
end
ElementABCD(m::Matrix) = ElementABCD(m[1,1], m[1,2], m[2,1], m[2,2])

"""

An optical element representing propagation over free space.

$(SIGNATURES)

The optical density is assumed to be unity. The propagation length is
the only field in this structure and hence the only argument to the
inherent constructor.

"""
struct FreeSpace{T,N} <: Element{T,N}
    L::T
end
FreeSpace(L) = FreeSpace{typeof(L), typeof(one(L))}(L)

"""

An optical Interface.

$(SIGNATURES)

The Interface has the parameters `η = n1 / n2` where `n1` and `n2` are
the optical densities (aka optical indices) of the previous resp. new
medium, the angle of incidence `aoi`, and a radius of curvature `roc`
which can be infinite. All arguments are optional but the default
arguments have no effect on a beam as the optical densities both
default to zero.

"""
struct Interface{L,N} <: Element{L,N}
    η::N; θ::N; R::L # R > 0 when light hits concave side
end
Interface(; n1=1.0, n2=1.0, η=n2/n1, aoi=0, roc=Inf) =
    Interface(η, aoi, roc)
@deprecate Interface(η,θ) Interface(n1=1,n2=η,aoi=θ)

"""

A thin lens.

$(SIGNATURES)

The parameters are focal length `f` and an optical angle of incidence
`aoi`.

"""
struct ThinLens{L,N} <: Element{L,N}
    f::L; θ::N
end
ThinLens(; f::L, aoi::N = zero(one(f))) where {L,N} =
    ThinLens{L,N}(f, aoi)
@deprecate ThinLens(focallength) ThinLens(f=focallength)

"""

A mirror.

$(SIGNATURES)

The optional keyword arguments are radius of curvature `roc` and angle
of incidence `aoi`.

"""
function Mirror(; roc=Inf, f=0.5*roc, aoi=0)
    if roc ≉ 2f
        throw(ArgumentError("roc and f are incompatible"))
    end
    return ThinLens(f=f, aoi=aoi)
end
@deprecate Mirror(R, θ) Mirror(roc=R, aoi=θ)

"""

Construct a pseudo-element for modeling beam propagation in the
tangential (aka parallel) plane.

$(SIGNATURES)

"""
struct Tan{E,T,N} <: Element{T,N} where {E<:Element{T,N}}
    e::E
end
Tan(e::FreeSpace) = e
Tan(e::Element{T,N}) where {T,N} = Tan{typeof(e),T,N}(e)
Tan(elements::Vector{<:Element}) = Tan.(elements)

"""

Construct a pseudo-element for modeling beam propagation in the
sagittal plane.

$(SIGNATURES)

"""
struct Sag{E,T,N} <: Element{T,N} where {E<:Element{T,N}}
    e::E
end
Sag(e::FreeSpace) = e
Sag(e::Element{T,N}) where {T,N} = Sag{typeof(e),T,N}(e)
Sag(elements::Vector{<:Element}) = Sag.(elements)

"""

Constructing a matrix results in the ray transfer matrix (also known
as ABCD matrix) representing the given, optical [`Element`](@ref) or
an entire system (a vector of `Element`).

$(SIGNATURES)

The matrix is represented as a julia `Matrix` with element type
`number` unless a specialization applies when using plain numbers
(without units).

"""
Base.Matrix(e::ElementABCD) = [e.A e.B; inv(e.invC) e.D]
Base.Matrix(e::FreeSpace) = [1 e.L ; zero(inv(e.L)) 1]
Base.Matrix(e::Interface) = [1 zero(e.R) ; (e.η-1)/e.R e.η]
Base.Matrix(e::ThinLens) = [1 zero(e.f) ; -inv(e.f) 1]
Base.Matrix(e::Tan{ThinLens{T,N},T,N}) where {T,N} =
    Matrix(ThinLens(; f = e.e.f * cos(e.e.θ), aoi = e.e.θ))
Base.Matrix(e::Sag{ThinLens{T,N},T,N}) where {T,N}  =
    Matrix(ThinLens(; f = e.e.f / cos(e.e.θ), aoi = e.e.θ))
function Base.Matrix(e::Tan{Interface{T,N},T,N}) where {T,N}
    # See doi:10.1364/AO.26.000427 for the following matrix
    θ1 = e.e.θ
    η = e.e.η
    R = e.e.R
    θ2 = asin(η * sin(θ1))
    return [
        cos(θ2)/cos(θ1) zero(L) ;
        (cos(θ2)-η*cos(θ1))/(R*cos(θ1)*cos(θ2)) η*cos(θ1)/cos(θ2)
    ]
end
function Base.Matrix(e::Sag{Interface{T,N},T,N}) where {T,N}
    # See doi:10.1364/AO.26.000427 for the following matrix
    θ1 = e.e.θ
    η = e.e.η
    R = e.e.R
    θ2 = asin(η * sin(θ1))
    return [one(N) zero(L) ; (cos(θ2)-η*cos(θ1))/R η]
end
Base.Matrix(elements::Vector{<:Element}) = prod(Matrix.(elements))
@deprecate RTM(element) Matrix(element)

"""

Return `η`, the ratio of optical densities (aka optical indices).

$(SIGNATURES)

The next optical density is in the numerator and the previous optical
density is in the denominator.

"""
η(e::Interface) = e.η
η(e::Union{Tan,Sag}) = η(e.e)
η(e::Element) = 1

"""

Return the effective beam propagation length of an element.

$(SIGNATURES)

The distance is measured along the beam's direction of propagation
(beam axis).

!!! note "To Do"
     The effective propagation length may deviate from the physical
     distance along the beam axis if the optical density deviates from
     unity. Currently, this is not taken into account. Do take it into
     account.

"""
dz(e::FreeSpace) = e.L
dz(e::Union{Tan,Sag}) = dz(e.e)
dz(e::Element{L,N}) where {L,N} = zero(L)

"""

Discretize a system by splitting each [`Element`](@ref) that occupies
space.

$(SIGNATURES)

Each element that occupies space is split into `N` appropriately
shortened versions of itself. A vector of elements is returned.

"""
discretize(e::FreeSpace, N::Int) = fill(FreeSpace(e.L/N), N)
discretize(e::Element, N::Int) = e
discretize(els::Vector{<:Element}, N::Int) = vcat(discretize.(els,N)...)
