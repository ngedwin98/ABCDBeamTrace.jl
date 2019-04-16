
"""

Return the location of a beam, measured along the beam axis.

$(SIGNATURES)

"""
function location(::AbstractBeam) end
location(Γ::GeometricBeam) = Γ.z
location(Γ::GaussianBeam) = location(Γ.b)

"""

Return the index of refraction (aka optical density).

$(SIGNATURES)

"""
function ior(::AbstractBeam) end
ior(Γ::GeometricBeam) = Γ.n
ior(Γ::GaussianBeam) = ior(Γ.b)

"""

Return the radial position of a beam.

$(SIGNATURES)

"""
function radialpos(::AbstractBeam) end
radialpos(Γ::GeometricBeam) = Γ.x
radialpos(Γ::GaussianBeam) = radialpos(Γ.b)

"""

Return the slope (paraxial angle, sine or tangens of that angle) of a
beam.

$(SIGNATURES)

"""
function slope(::AbstractBeam) end
slope(Γ::GeometricBeam) = Γ.k
slope(Γ::GaussianBeam) = slope(Γ.b)
