abstract type Element end
struct FreeSpace <: Element
    L::Real
end
struct Interface <: Element
    η::Real; θ::Real
end
struct ThinLens <: Element
    f::Real; θ::Real
end
Interface(η::Real) = Interface(η,0)
ThinLens(f::Real) = ThinLens(f,0)
Mirror(R::Real,θ::Real) = ThinLens(R/2,θ)
FreeSpace(L::Real, N::Int) = fill(FreeSpace(L/N), N)
FreeSpace(L::Real, dz::Real) = FreeSpace(L, (Int∘round)(L/dz))

struct Tan{E} <: Element where {E<:Element}
    e::E
end
struct Sag{E} <: Element where {E<:Element}
    e::E
end
Tan(e::FreeSpace) = e
Sag(e::FreeSpace) = e
Sag(e::Interface) = e

RTM(e::FreeSpace) = [1 e.L ; 0 1]
RTM(e::Interface) = [1 0 ; 0 e.η]
RTM(e::ThinLens) = [1 0 ; -1/e.f 1]
RTM(e::Tan{ThinLens}) = (RTM∘ThinLens)(e.e.f*cos(e.e.θ))
RTM(e::Sag{ThinLens}) = (RTM∘ThinLens)(e.e.f/cos(e.e.θ))
function RTM(e::Tan{Interface})
    θ1, η = e.e.θ, e.e.η; θ2 = asin(η*sin(θ1))
    return [cos(θ2)/cos(θ1) 0 ; 0 η*cos(θ1)/cos(θ2)]
end

dΓ(e::Element) = (0,1)
dΓ(e::Union{Tan,Sag}) = dΓ(e.e)
dΓ(e::FreeSpace) = (e.L,1)
dΓ(e::Interface) = (0,e.η)
