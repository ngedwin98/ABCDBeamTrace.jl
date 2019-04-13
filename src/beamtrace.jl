struct Beam
    λ::Real; z::Real; n::Real
    x::Real; k::Real; q::Complex
end
Beam(λ::Real, w0::Real, n0::Real=1) = Beam(λ, 0, n0, 0, 0, 1im*π*n0*w0^2/λ)

η(e::Interface) = e.η
η(e::Union{Tan,Sag}) = dη(e.e)
η(e::Element) = 1

dz(e::FreeSpace) = e.L
dz(e::Union{Tan,Sag}) = dz(e.e)
dz(e::Element) = 0

function transform(e::Element, Γ::Beam)
    M = RTM(e)
    return Beam(Γ.λ, Γ.z+dz(e), Γ.n/η(e), M*[Γ.x,Γ.k]..., /(M*[Γ.q,1]...))
end

function beamtrace(elems::Vector{<:Element}, Γ0::Beam)
    Γs = Vector{Beam}(undef, length(elems)+1); Γs[1] = Γ0
    for (ind, elem) in enumerate(elems)
        Γs[ind+1] = transform(elem, Γs[ind])
    end
    return Γs
end

spotsize(Γ::Beam) = /(-Γ.λ, π*Γ.n*imag(1/Γ.q)) |> sqrt
location(Γ::Beam) = Γ.z

"""

Return a function that calculates the beam's waist radius as a
function of the propagation distance.

Expects an optical system (e.g. a vector or collection of
[`Element`](@ref)) as first argument, and a [`Beam`](@ref) as second
argument. The optional keyword argument `outside` can be used to
define a return value for beam positions outside those covered by
the system; the default is to throw a DomainError.

"""
function waistradiusfunc(elements, beam::Beam; outside=nothing)
    beams = Vector{Beam}(undef, length(elements)+1)
    beams[1] = beam
    for (i, el) in enumerate(elements)
        beam = transform(el, beam)
        beams[i+1] = beam
    end
    return function(z)
        if z < beams[1].z || z > beams[length(beams)].z
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
        while beams[i+1].z < z
            i += 1
        end
        i = max(i, 1)
        # add the effect of a FreeSpace element to reach the requested
        # beam position z
        beam = transform(FreeSpace(z - beams[i].z), beams[i])
        return spotsize(beam)
    end
end

discretize(e::FreeSpace, N::Int) = fill(FreeSpace(e.L/N), N)
discretize(e::Element, N::Int) = e
discretize(els::Vector{<:Element}, N::Int) = vcat(discretize.(els,N)...)
