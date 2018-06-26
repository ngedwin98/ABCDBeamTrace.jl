struct Beam
    λ::Real; z::Real; n::Real
    x::Real; k::Real; q::Complex
end
Beam(λ::Real, w0::Real, n0::Real=1) = Beam(λ, 0, n0, 0, 0, 1im*π*n0*w0^2/λ)

function Mobius(e::Element)
    M = RTM(e); dz,η = dΓ(e)
    return function(Γ::Beam)
        return Beam(Γ.λ, Γ.z+dz, Γ.n/η, M*[Γ.x,Γ.k]..., /(M*[Γ.q,1]...))
    end
end

function beamtrace(elems::Vector{<:Element}, Γ0::Beam)
    Γ = Vector{Beam}(length(elems)+1); Γ[1] = Γ0
    for (ind, e) in enumerate(elems)
        Γ[ind+1] = Mobius(e)(Γ[ind])
    end
    return Γ
end

discretize(e::Element, N::Int) = e
discretize(e::FreeSpace, N::Int) = FreeSpace(e.L, N)
discretize(els::Vector{<:Element}, N::Int) = vcat([discretize.(e,N) for e in els]...)

spotsize(Γ::Beam) = /(-Γ.λ, π*Γ.n*imag(1/Γ.q)) |> sqrt
location(Γ::Beam) = Γ.z
