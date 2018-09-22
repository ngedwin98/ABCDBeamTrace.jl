abstract type Element end
struct FreeSpace <: Element
    L::Real
end
struct Interface <: Element
    η::Real; θ::Real; R::Real # R > 0 when light hits concave side
end
struct ThinLens <: Element
    f::Real; θ::Real
end
ThinLens(f::Real) = ThinLens(f,0)
Interface(η::Real) = Interface(η,0,Inf)

Mirror(R::Real, θ::Real) = ThinLens(R/2,θ)

struct Tan{E} <: Element where {E<:Element}
    e::E
end
struct Sag{E} <: Element where {E<:Element}
    e::E
end
Tan(e::FreeSpace) = e
Sag(e::FreeSpace) = e

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
