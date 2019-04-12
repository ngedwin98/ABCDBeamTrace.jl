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

discretize(e::FreeSpace, N::Int) = fill(FreeSpace(L/N), N)
discretize(e::Element, N::Int) = e
discretize(els::Vector{<:Element}, N::Int) = vcat(discretize.(e,N)...)
