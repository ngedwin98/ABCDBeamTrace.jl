# Please see the new maintained fork, which follows similar principles: [JuliaPhysics/ABCDMatrixOptics.jl](https://github.com/JuliaPhysics/ABCDMatrixOptics.jl) 

# ABCDBeamTrace.jl
A Julia package for performing calculations with the [ray transfer matrix (or ABCD) formalism](https://en.wikipedia.org/wiki/Ray_transfer_matrix_analysis), for both 1D ray tracing and [Gaussian beam](https://en.wikipedia.org/wiki/Gaussian_beam) propagation in the [paraxial approximation](https://en.wikipedia.org/wiki/Paraxial_approximation).

The following introduction to the package assumes familiarity with the ABCD formalism and its utility in optical analysis and design.  In addition to the above links, the following are classic and useful introductory references:
1. H. Kogelnik and T. Li, "Laser Beams and Resonators", *Applied Optics* **5**, 1550-1567 (1966)
2. A. E. Siegman, *Lasers* (University Science Books, Sausalito, 1986)

## Installation
This package is not yet registered.  It can be installed in Julia with the following:
```julia
Pkg.clone("git://github.com/ngedwin98/ABCDBeamTrace.jl.git")
```
The code currently requires Julia `v0.6`.

## Optical elements
This section discusses how we represent the ABCD model of optical elements and how to access their standard representation as numerical ray transfer matrices.

### Fundamental elements
The basic type representing an optical element is `Element`.  As an abstract type, its implemented concrete types are (all physical distances taken to be in some consistent, arbitrary unit):
* `FreeSpace(L::Real)`: Represents a section of free space of distance represented by `L`.
* `ThinLens(f::Real)`: Represents a focusing element that behaves like a thin lens of focal length represented by `f`.  Optionally, a tilt angle represented by `θ` in radians can be specified using `ThinLens(f::Real,θ::Real)`; otherwise it defaults to `θ=0`.
* `Interface(η::Real)`: Represents an dielectric interface where `η` is the ratio of the initial to final refractive indices.  Optionally, a tilt angle represented by `θ` in radians and a radius of curvature (positive when concave side faces incoming beam) represented by `R` can be specified using `Interface(η::Real,θ::Real,R::Real)`.  Note that if either of these two are specified, they must both be specified, and otherwise they default to `θ=0` and `R=Inf`.

### Additional elements
There are some additional optical elements which are physically distinct but which, in the paraxial approximation, can be represented in terms of the above elements.  Currently, the following convenience methods allow such "elements" to be constructed:
* `Mirror(R::Real,θ::Real) = ThinLens(R/2,θ)` represents a curved mirror with radius of curvature (positive when concave side faces incoming beam) represented by `R` and a tilt angle represented by `θ` in radians.

### Optical systems
An optical system in the ABCD formalism is a cascade of optical elements, represented in this package as `Vector{<:Element}`.  For example, to construct a [Keplerian 2x beam expander](https://www.edmundoptics.com/resources/application-notes/lasers/beam-expanders/) with objective focal length represented by `f::Real`, we can use
```julia
expander_2x = [ThinLens(f), FreeSpace(3f), ThinLens(2f)]
```
This means that we can also compose two optical systems using vector concatenations.  For example, given any length represented by `L::Real`, the following effectively creates a 1-to-1 beam expander with four lenses:
```julia
system = [expander_2x; FreeSpace(L); reverse(expander_2x)]
```
Note the use of `vcat` in composing systems together.

### Ray transfer (ABCD) matrices
The usual ABCD formalism for calculations involve the multiplication of 2x2 matrices, as discussed in the references above. To facilitate this formalism, each element has a well-defined matrix, which is accessed by calling `RTM(e::Element)`.  For example, `RTM(ThinLens(100)) == [1 0 ; -1/100 1]` evaluates as true.  Dot syntax for broadcasting can be used to create a corresponding vector of matrices, and the corresponding elements can be multiplied up, as in
```julia
# Total ABCD matrix for above 4-lens system
system_RTM = reduce(*, RTM.(system))
```

### Sagittal and tangential elements
In the ABCD formalism, optical components that can have a non-zero tilt relative to the optical axis (i.e., with a nonzero field `θ`) bifurcate into two distinct elements, acting differently on the sagittal and tangential components of the ray or Gaussian beam.  To represent this dual behavior, composite `Element` types called `Tan` and `Sag` wrap the fundamental `Element` types and dispatch differently on `RTM` in order to produce the respective ray transfer matrices for the tangential and sagittal components.

For example, suppose the above beam expander were misaligned with a tilt angle represented by `θ::Real`.  We would represent the optical system as
```julia
system = [ThinLens(f,θ), FreeSpace(3f), ThinLens(2f,θ)] # Component-agnostic description of system
system_tan = Tan.(system) # Elements as seen by tangential component
system_sag = Sag.(system) # Elements as seen by sagittal component
system_tan_RTM = reduce(*, RTM.(system_tan)) # Total matrix in the tangential component
system_sag_RTM = reduce(*, RTM.(system_sag)) # Total matrix in the sagittal component
```

**Note**: If no calls are made to either `Tan` and `Sag` (as in the first line above), the `θ` field (if present) confers *no effect*; in this case, `RTM` returns the ray transfer matrix as if we set `θ=0`.  Also, note that `FreeSpace` does not wrap into these composite objects: `Tan(e::FreeSpace) = Sag(e::FreeSpace) = e`.  (This latter behavior could be changed, should there be a reason for it to be otherwise...)

## Beam propagation
The numerical ray transfer matrices in the ABCD formalism facilitate a transfer-function-based approach to optical analysis.  However, it is also possible to take a *state-based* viewpoint, by representing the state of the ray or beam as it passes through the elements.  Here, the parameters that describe the state of the beam are modified by each element, causing the beam to evolve as it propagates through the elements.  This approach is most useful for visualizing the optical propagation of the beam, or for optimizing beam properties in the design.

### Beam parameters
Consistent with the ABCD formalism, we provide a description the state of a Gaussian beam (and the ray parameters) at a particular point along the optical axis.  It is represented with a concrete type `Beam` containing the following fields:
* `q`: The complex Gaussian beam parameter
* `x`: The transverse component of the ray perpendicular to the optical axis
* `k`: The slope of the ray relative to the optical axis
* `z`: The location at which the state is being described, along the optical axis
* `n`: The refractive index at this location
* `λ`: The physical wavelength of the (monochromatic) beam

All optical elements modify `x`, `k`, and `q`, as prescribed by the ABCD formalism.  The parameter `z` is only modified by `FreeSpace` (via addition by the field `L`), while the parameter `n` is only modified by `Interface` (via division by the field `η`).  The wavelength `λ` does not play a direct role in the (scale-free) beam propagation; it is only used for spot size calculations.

A useful constructor method is `Beam(λ::Real,w0::Real,n0::Real=1)`, which represents a beam state with position represented by `z=0` in refractive index represented by `n0` and at a focus with 1/e² waist represented by `w0` (so that the Gaussian beam parameter is imaginary), provided its wavelength `λ`; furthermore, the ray position and slope are represented by `x=0` and `k=0`, respectively.

### Beam transformation by an element
We think of the optical elements as effecting a transformation on the beam state.  If we have an optical element represented by `e::Element` and a beam state `Γ::Beam`, there is a non-exported function `transform(e::Element,Γ::Beam)`, which returns a new instance of `Beam` with fields representing the beam state after propagation through that element.  This function is unexported because it should rarely see utility given the beamtracing functionality described in the next subsection.

**Note**: It is very natural to use Julia's functor mechanism to represent the transformation, via some definition like `(e::Element)(Γ::Beam)` which behaves identically to the currently-implemented `transform(e,Γ)`.  Unfortunately, this is not possible at the moment due to https://github.com/JuliaLang/julia/issues/14919.  While an obvious workaround is to define the transformation for each concrete subtype, this is inelegant enough to be not worth it.

### Beam and ray tracing
Perhaps the most useful part of this whole package is the ability to trace the beam state by providing an initial beam state and propagating it forwards using the optical elements of a given system, recording the state after each element.  For an optical system represented by `elems::Vector{<:Element}` and an initial beam state represented by `Γ0::Beam`, `beamtrace(elems::Vector{<:Element},Γ0::Beam)` returns an instance of `Vector{Beam}` (of length given by `length(elems)+1`) where the first item is `Γ0`, and each subsequent item is the result of applying `transform` to the previous item using the corresponding item of `elems`.

## Examples
More detailed examples are upcoming.  In particular, this package has been developed with the following tasks in mind:
* Cavity design
* Mode-matching

Jupyter notebooks implementing these tasks exist and are currently being cleaned so that they can be pushed to this repository to serve as examples.

Other interesting examples using this package are welcome!

## To Do
* Add examples of beam and ray tracing function calls
* Add detailed Jupyter notebooks showing examples
* Systematic tests of correctness and consistency
