"""

An abstract type as supertype of GaussianBeam and GeometricBeam. It is
parametrized on types L, CL, and N for quantities with dimensions
length, complex length, and dimensionless, respectively.

"""
abstract type AbstractBeam{L,CL,N} end
