"""

Construct a Beam object for modeling the propagation of a Gaussian
(i.e. typical laser) beam.

$(SIGNATURES)

This type is parametrized on three underlying types, `L` (for lengths
such as the wavelength or the location along the beam axis), `CL` (for
complex-valued lengths such as the beam parameter), and `N` (for
dimensionless numbers such as the optical density and a normalization
constant `k`).

"""
struct GaussianBeam{L,CL,N} <: AbstractBeam{L,CL,N}
    b::GeometricBeam{L,CL,N} # contains all the uncommented fields
                             # below
    λ::L # vacuum wavelength
    #z::L # position along beam axis
    #n::N # inde of refraction aka optical density
    #x::L # radial extent
    #k::N # (angle/sin/tan of) slope of beam
    q::CL # complex beam parameter
end
function GaussianBeam(
    ;
    q = nothing, # complex beam parameter
    λ::Number, # vacuum wavelength
    w0 = nothing, # waist radius
    z0::Number = zero(λ), # waist position along beam axis
    n::Number = 1 # ior: index of refraction (aka optical density)
)
    # determine dimensional compatibility of types using the trick
    # that zero(one(x)) == zero(x) if and only if x is dimensionless
    if zero(one(n)) != zero(n)
        throw(DomainError("n must be dimensionless"))
    end
    if isnothing(w0) == isnothing(q)
        throw(DomainError("either q or w0 must be specified"))
    end
    if (z0 != zero(λ)) && !isnothing(q)
        throw(DomainError("q and z0 cannot both be specified"))
    end
    if isnothing(w0)
        rr = rayleighrange(q)
        w0 = sqrt(rr * λ / (π * n))
        z0 = waistdistance(q)
    end
    if zero(one(w0 / λ)) != zero(w0 / λ)
        throw(DomainError("ratios of lengths must be dimensionless"))
    end
    if zero(one(z0 / λ)) != zero(z0 / λ)
        throw(DomainError("ratios of lengths must be dimensionless"))
    end
    # promote arguments, prepending p for "promoted" to variable
    # names; multiply by 1.0 to promote to float (without throwing an
    # error for symbolic variables that do not have a float method
    # defined)
    pλ, pw0, pz0 = promote(1.0λ, 1.0w0, 1.0z0)
    pn = 1.0n
    # determine values
    pz = zero(pλ) # position along beam axis (zero)
    pzR = π * pn * pw0^2 / pλ # Rayleigh length (aka Rayleigh range)
    # determine types
    L = typeof(pλ)
    CL = complex(L)
    N = typeof(one(pλ))
    return GaussianBeam{L,CL,N}(
        GeometricBeam{L,CL,N}(
            pz,
            n,
            zero(pλ),
            zero(pn)
        ),
        pλ,
        pz - pz0 + 1im * pzR
    )
end
@deprecate Beam(λ_beam, w0_beam) GaussianBeam(λ=λ_beam, w0=w0_beam)
@deprecate Beam(λ_beam, w0_beam, n0) GaussianBeam(
    λ=λ_beam,w0=w0_beam, n=n0
)

"""

Propagate a [`GaussianBeam`](@ref) or complex beam parameter through
an [`Element`](@ref), a System (represented by a vector of elements),
or through a ray transfer matrix.

$(SIGNATURES)

Note that the [`GaussianBeam`](@ref) `Γ` (or the complex beam
parameter `q`) is the second argument. If a ray transfer matrix is
given as first argument, observe the optional keyword arguments to
pass a propagation distance `dz` and a ratio `η` of new to old
refractive index.

"""
function transform(bywhat::Any, ::AbstractBeam) end
function transform(system::Vector{<:Element}, Γ::AbstractBeam)
    for i = 1:length(system)
        Γ = transform(system[i], Γ)
    end
    return Γ
end
transform(e::Element, Γ::AbstractBeam) =
    transform(Matrix(e), Γ; dz = dz(e), η = η(e))
transform(m::Matrix, q) = /((m * [q, 1])...) # beam parameter q
function transform(
    m::Matrix,
    Γ::GaussianBeam;
    dz = zero(m[1,1]),
    η = one(m[1,1])
)
    geometric_beam = transform(m, Γ.b; dz = dz, η = η)
    q = transform(m, Γ.q)
    L = typeof(real(q))
    CL = typeof(q)
    N = typeof(ior(geometric_beam))
    return GaussianBeam{L,CL,N}(
        GeometricBeam{L,CL,N}(geometric_beam),
        Γ.λ,
        q
    )
end

"""

Propagate a [`GaussianBeam`](@ref) through a system.

$(SIGNATURES)

The first argument, the optical system, is given as a vector of
[`Element`](@ref). The return value is a vector consisting of the
unpropagated beam, and the beams propagated through the first one, the
first two, etc., and, finally, all elements of which the system
consists.

Note that the [`GaussianBeam`](@ref) is the second argument.

"""
function beamtrace(elems::Vector{<:Element}, Γ0::GaussianBeam)
    Γs = Vector{GaussianBeam}(undef, length(elems)+1)
    Γs[1] = Γ0
    for (ind, elem) in enumerate(elems)
        Γs[ind+1] = transform(elem, Γs[ind])
    end
    return Γs
end

"""

Return the complex beam parameter (usually denoted `q`) of a
[`GaussianBeam`](@ref), a ray transfer matrix, an [`Element`](@ref) or
a system (a `Vector` of element type [`Element`](@ref)).

$(SIGNATURES)

If an optical system is given, it is assumed to describe one resonator
round-trip; the resulting beam parameter of the fundamental resonant
mode (a Gaussiam beam) is returned or `nothing` if none exists.

This function is also defined (as a stub) for an
[`AbstractBeam`](@ref) to allow other packets to provide
implementations.

!!! note
    This method, along with most others that may seem like they make
    sense only for a [`GaussianBeam`](@ref) is also defined for a
    [`GeometricBeam`](@ref), returning the limit of simultaneously
    vanishing wavelength, waist radius, and Rayleigh range. This
    leads to reasonable return values for these methods.

"""
function beamparameter(::AbstractBeam) end
beamparameter(beam::GaussianBeam) = beam.q
beamparameter(beam::GeometricBeam) = complex(location(beam))
beamparameter(element::Element) = beamparameter(Matrix(element))
beamparameter(system::Vector{<:Element}) =
    beamparameter(Matrix(system))
function beamparameter(m::Matrix)
    A = m[1,1]
    B = m[1,2]
    C = m[2,1]
    D = m[2,2]
    if A isa Number && D isa Number && 4 - (A + D)^2 < 0
        # no stable eigenmode
        return nothing
    end
    # Eq. 57 of Kogelnik and Li, doi:10.1364/AO.5.001550)
    re_inv_q = (D-A) / (2B)
    im_inv_q = sqrt(4 - (A + D)^2) / (2B)
    # decide sign in re_inv_q ± im_inv_q by which leads to a real
    # beam width (see remark after Eq. 46): im_inv_q > 0
    inv_q = complex(re_inv_q, -abs(im_inv_q))
    return inv(inv_q)
end

"""

Return the beam parameter product (usually denoted ``M^2``) of an
[`AbstractBeam`](@ref).
x
$(SIGNATURES)

!!! info "To Do"
    The intention is to use this in generic implementations of
    e.g. [`spotradius`](@ref) to make the code more reusable
    for any package that provides a custom, concrete type derived
    from [`AbstractBeam`](@ref). This still needs to be done.

"""
function beamparameterproduct(::AbstractBeam) end
beamparameterproduct(::GaussianBeam{L,CL,N}) where {L,CL,N} = one(N)
beamparameterproduct(::GeometricBeam{L,CL,N}) where {L,CL,N} = zero(N)

"""

Return the wavefront's radius of curvature of an
[`AbstractBeam`](@ref) or a complex beam parameter `q`.

$(SIGNATURES)

"""
wavefrontroc(beam::AbstractBeam) = wavefrontroc(beamparameter(beam))
wavefrontroc(q) = inv(real(inv(q)))

"""

Return the Rayleigh range (usually denoted ``z_R``) of an
[`AbstractBeam`](@ref) or complex [`beamparameter`](@ref).

$(SIGNATURES)

"""
rayleighrange(beam::AbstractBeam) = rayleighrange(beamparameter(beam))
rayleighrange(q) = imag(q)

"""

Return the waist location (usually denoted `z_0`) of an
[`AbstractBeam`](@ref).

$(SIGNATURES)

!!! note
    This method does not take care of boundaries set by
    elements: It reports a waist location assuming that
    there is an infinite amount of space for the beam
    to propagate, in both forward and backward direction.

!!! warning
    Unlike the other methods that work for a [`GaussianBeam`](@ref),
    this method is not also available for a complex
    [`beamparameter`](@ref) as an argument. The design rationale is
    that the complex beam parameters used in ray transfer matrix
    analysis of Gaussian beams typically do not encode a location but
    merely a distance from the last optical element: Use
    [`waistdistance`](@ref) instead if you only have a complex beam
    parameter rather than an [`AbstractBeam`](@ref).

"""
waistlocation(beam::AbstractBeam) =
    location(beam) + waistdistance(beamparameter(beam))

"""

Return the distance from the current location to the free beam waist
(usually denoted `z_0`) of an [`AbstractBeam`](@ref) or complex
[`beamparameter`](@ref).

$(SIGNATURES)

!!! note
    This method does not take care of boundaries set by
    elements: It reports a waist location assuming that
    there is an infinite amount of space for the beam
    to propagate, in both forward and backward direction.

"""
waistdistance(beam::AbstractBeam) = waistdistance(beamparameter(beam))
waistdistance(q) = -real(q)

"""

Return the ``1/e^2`` radius of the beam waist of a
[`GaussianBeam`](@ref).

$(SIGNATURES)

Note that the beam waist may occur outside of physical
constraints. Its radius is calculated as if it would really occur.

"""
function waistradius(::AbstractBeam) end
function waistradius(beam::GaussianBeam; none = :error)
    rr = rayleighrange(beam)
    # should none be returned (instead of letting sqrt throw a
    # DomainError)?
    if rr isa Number && rr < zero(rr) && none !== :error
        return none
    end
    return sqrt(rr * beam.λ / (π * ior(beam)))
end

"""

Return the ``1/e^2`` radius of a beam at its current location.

$(SIGNATURES)

At the radius returned by this function, the intensity drops to
``1/e^2``. If the keyword argument `none` is specified, return
that argument rather than throw a `DomainError` if the beam is
invalid, i.e. has no (real and finite) spotradius.

"""
function spotradius(Γ::GaussianBeam; none = :error)
    wsquared = -Γ.λ / (π * ior(Γ) * imag(inv(beamparameter(Γ))))
    if isnothing(none) || none !== :error
        # return a custom value if there is no spotradius, rather than
        # the DomainError that would normally result (but only if the
        # result is numeric, i.e. wsquared is a Number, as isless will
        # likely not be implemented for a symbolic result)
        if wsquared isa Number
            if wsquared < zero(wsquared) || isinf(wsquared)
                return none
            end
        end
    end
    return sqrt(wsquared)
end
@deprecate spotsize(beam) spotradius(beam)

"""

Return a function that calculates the beam's waist radius as a
function of the propagation distance.

Expects an optical system (e.g. a vector or collection of
[`Element`](@ref)) as first argument, and a [`GaussianBeam`](@ref) as
second argument. The optional keyword argument `outside` can be used
to define a return value for beam positions outside those covered by
the system; the default is to throw a DomainError.

"""
function spotradiusfunc(elements, beam::GaussianBeam; outside=nothing)
    beams = beamtrace(elements, beam)
    return function(z)
        if (
            z < location(beams[1]) ||
            z > location(beams[length(beams)])
        )
            if isnothing(outside)
                throw(DomainError(string(
                    "system does not cover ",
                    "the requested beam position"
                )))
            end
            return outside
        end
        # find last element before requested beam position z; note
        # that this could be accelerated by a more intelligent
        # algorithm since the beams[i].z should be sorted (and this
        # algorithm already depends on this expected ordering).
        i = 0
        while location(beams[i+1]) < z
            i += 1
        end
        i = max(i, 1)
        # add the effect of a FreeSpace element to reach the requested
        # beam position z
        beam = transform(FreeSpace(z - location(beams[i])), beams[i])
        return spotradius(beam)
    end
end
