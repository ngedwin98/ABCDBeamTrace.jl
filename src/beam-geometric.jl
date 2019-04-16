"""

A beam obeying geometric optics.

$(SIGNATURES)

Its properties are described by a position `z` along the optical axis,
the current index of refraction `n`, a radial position (or extent)
`x`, and a beam slope `k`. The constructor takes some of these as
optional keyword arguments, namely `n`, `x` and a `slope` (which
internally becomes `k`).

!!! note
    Do not confuse the `slope` (or, internally, `k`) with the wave
    vector ``\\vec{k}``.

"""
struct GeometricBeam{L,CL,N} <: AbstractBeam{L,CL,N}
    z::L # position along beam axis
    n::N # inde of refraction aka optical density
    x::L # radial extent
    k::N # (angle/sin/tan of) slope of beam
end
function GeometricBeam(
    ;
    x::L = 0.0,
    n::N = 1.0,
    slope::N = zero(n)
) where {L,N}
    if zero(one(n)) != zero(n)
        throw(DomainError("n must be dimensionless"))
    end
    if zero(one(slope)) != zero(slope)
        throw(DomainError("slope must be dimensionless"))
    end
    return GeometricBeam{L,complex(L),N}(
        zero(x),
        n,
        x,
        slope
    )
end

GeometricBeam{L,CL,N}(beam::GeometricBeam) where {L,CL,N} =
    GeometricBeam{L,CL,N}(
        L(beam.z),
        N(beam.n),
        L(beam.x),
        N(beam.k)
    )

"""

Propagate a [`GeometricBeam`](@ref) through an [`Element`](@ref).

$(SIGNATURES)

"""
transform(
    m::Matrix,
    Γ::GeometricBeam{L,CL,N};
    dz = zero(m[1,1]),
    η = one(m[1,1])
) where {L,CL,N} = GeometricBeam{
    promote_type(L, typeof(m[1,2])),
    complex(promote_type(L, typeof(m[1,2]))),
    promote_type(L, typeof(m[1,1]), typeof(m[2,2]))
}(
    Γ.z + dz,
    Γ.n / η,
    (m * [Γ.x, Γ.k])...
)

"""

Propagate an [`AbstractBeam`](@ref) through a system.

$(SIGNATURES)

The first argument, the optical system, is given as a vector of
[`Element`](@ref). The return value is a vector consisting of the
unpropagated beam, and the beams propagated through the first one, the
first two, etc., and, finally, all elements of which the system
consists.

Note that the beam is the second argument.

"""
function beamtrace(::Vector{<:Element}, ::AbstractBeam) end
function beamtrace(elems::Vector{<:Element}, Γ0::GeometricBeam)
    Γs = Vector{GeometricBeam}(undef, length(elems)+1)
    Γs[1] = Γ0
    for (idx, elem) in enumerate(elems)
        Γs[idx+1] = transform(elem, Γs[idx])
    end
    return Γs
end
