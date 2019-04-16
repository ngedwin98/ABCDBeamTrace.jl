# ABCDBeamTrace.jl
```@meta
CurrentModule = ABCDBeamTrace
```

A Julia package for performing calculations with the [ray transfer
matrix (or ABCD)
formalism](https://en.wikipedia.org/wiki/Ray_transfer_matrix_analysis)
(for which wikipedia has a [second, complimentary
page](https://en.wikipedia.org/wiki/Transfer-matrix_method_(optics))),
for both 1D ray tracing and [Gaussian
beam](https://en.wikipedia.org/wiki/Gaussian_beam) propagation in the
[paraxial
approximation](https://en.wikipedia.org/wiki/Paraxial_approximation). A
concise if incomplete introduction can also be found at [the RP
Photonics
Encyclopedia](https://www.rp-photonics.com/abcd_matrix.html).

The following introduction to the package assumes familiarity with the
ABCD formalism and its utility in optical analysis and design.  In
addition to the above links, the following are classic and useful
introductory references:
1. H. Kogelnik and T. Li, "Laser Beams and Resonators", *Applied Optics* **5**, 1550-1567 (1966) (currently available free of charge online at [doi:10.1364/AO.5.001550](https://doi.org/10.1364/AO.5.001550))
2. A. E. Siegman, *Lasers* (University Science Books, Sausalito, 1986) (see e.g. the [google books preview](https://books.google.de/books?id=1BZVwUZLTkAC&printsec=frontcover))

# Usage

!!! note
    This section is outdated:
    The described methods work but some of them produce
    deprecation warnings (that mention the future syntax
    with keyword arguments replacing positional arguments
    for more clarity and less potential for error).

## Installation
This package is not yet registered.  It can be installed in Julia with the following:
```julia
Pkg.clone("git://github.com/ngedwin98/ABCDBeamTrace.jl.git")
```
The code currently requires Julia `v1.0` or higher. It has been tested with `v1.1`.

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
```jldoctest usage
using ABCDBeamTrace
f=125e-3
expander_2x = [ThinLens(f=f), FreeSpace(3f), ThinLens(f=2f)]

# output

3-element Array{ABCDBeamTrace.Element{Float64,Float64},1}:
 ThinLens{Float64,Float64}(0.125, 0.0)
 FreeSpace{Float64,Float64}(0.375)
 ThinLens{Float64,Float64}(0.25, 0.0)
```

This means that we can also compose two optical systems using vector
concatenations.  For example, the following continuation of the
example above effectively creates a 1-to-1 beam expander with four
lenses:

```jldoctest usage
L = 1000e-3
system = [expander_2x; FreeSpace(L); reverse(expander_2x)]

# output

7-element Array{ABCDBeamTrace.Element{Float64,Float64},1}:
 ThinLens{Float64,Float64}(0.125, 0.0)
 FreeSpace{Float64,Float64}(0.375)
 ThinLens{Float64,Float64}(0.25, 0.0)
 FreeSpace{Float64,Float64}(1.0)
 ThinLens{Float64,Float64}(0.25, 0.0)
 FreeSpace{Float64,Float64}(0.375)
 ThinLens{Float64,Float64}(0.125, 0.0)
```

Note the implicit use of `vcat` (via the semicolon notation) in
composing systems together.

### Ray transfer (ABCD) matrices
The usual ABCD formalism for calculations involve the multiplication
of 2x2 matrices, as discussed in the references above. To facilitate
this formalism, each element has a well-defined matrix, which is
accessed by calling `Matrix(e::Element)`.  For example,
`Matrix(ThinLens(100)) == [1 0 ; -1/100 1]` evaluates as true.  Dot
syntax for broadcasting can be used to create a corresponding vector
of matrices, and the corresponding elements can be multiplied up, as
in

```jldoctest usage
# Total ABCD matrix for above 4-lens system
system_RTM = Matrix(system)

# output

2×2 Array{Float64,2}:
 1.0  -0.125
 0.0   1.0
```

### Sagittal and tangential elements
In the ABCD formalism, optical components that can have a non-zero
tilt relative to the optical axis (i.e., with a nonzero field `θ` for
the angle of incidence) bifurcate into two distinct elements, acting
differently on the sagittal and tangential components of the ray or
Gaussian beam.  To represent this dual behavior, composite
[`Element`](@ref) types called [`Tan`](@ref) and [`Sag`](@ref) wrap
the fundamental [`Element`](@ref) types and dispatch differently on
[`Base.Matrix`](@ref) in order to produce the respective ray transfer
matrices for the tangential and sagittal components.

For example, suppose the above beam expander were misaligned with a
tilt angle represented by `θ::Real`.  We would represent the optical
system as
```jldoctest usage
f = 125e-3 # focal length of 125mm in SI base units (m)
θ = 1.0 * π/180 # one degree of misalignment in radians
# Plane-agnostic description of system
system = [ThinLens(f=f,aoi=θ), FreeSpace(3f), ThinLens(f=2f,aoi=θ)]
system_tan = Tan(system) # Elements as seen by tangential component
system_sag = Sag(system) # Elements as seen by sagittal component
system_tan_RTM = Matrix(system_tan) # Total matrix in the tangential component
system_sag_RTM = Matrix(system_sag) # Total matrix in the sagittal component
system_tan_RTM, system_sag_RTM

# output

([-0.500228 0.375; 0.00182821 -2.00046], [-0.499772 0.375; -0.00182738 -1.99954])
```

!!! note
    If no calls are made to either `Tan` and `Sag` (as in the first
    line above), the `θ` field (if present) confers *no effect*; in
    this case, `Base.Matrix` returns the ray transfer matrix as if we
    set `θ=0`.  Also, note that `FreeSpace` does not wrap into these
    composite objects: `Tan(e::FreeSpace) = Sag(e::FreeSpace) = e`.
    (This latter behavior could be changed, should there be a reason
    for it to be otherwise...)

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

A useful constructor method is [`GaussianBeam`](@ref) with keyword
arguments `λ`, `w0`, `n`, `z0`. It represents a beam state with
position represented by `z=-z0` in refractive index represented by `n`
and at a focus with 1/e² waist represented by `w0` (so that the
Gaussian beam parameter is imaginary for `z0=0` which is the default),
provided its wavelength `λ`; furthermore, the ray position and slope
are represented by `x=0` and `k=0`, respectively.

### Beam transformation by an element
We think of the optical elements as effecting a transformation on the beam state.  If we have an optical element represented by `e::Element` and a beam state `Γ::GaussianBeam`, there is a non-exported function `transform(e::Element,Γ::GaussianBeam)`, which returns a new instance of `GaussianBeam` with fields representing the beam state after propagation through that element.  This function is unexported because it should rarely see utility given the beamtracing functionality described in the next subsection.

**Note**: It is very natural to use Julia's functor mechanism to represent the transformation, via some definition like `(e::Element)(Γ::GaussianBeam)` which behaves identically to the currently-implemented `transform(e,Γ)`.  Unfortunately, this is not possible at the moment due to https://github.com/JuliaLang/julia/issues/14919.  While an obvious workaround is to define the transformation for each concrete subtype, this is inelegant enough to be not worth it.

### Beam and ray tracing
Perhaps the most useful part of this whole package is the ability to trace the beam state by providing an initial beam state and propagating it forwards using the optical elements of a given system, recording the state after each element.  For an optical system represented by `elems::Vector{<:Element}` and an initial beam state represented by `Γ0::GaussianBeam`, `beamtrace(elems::Vector{<:Element},Γ0::GaussianBeam)` returns an instance of `Vector{GaussianBeam}` (of length given by `length(elems)+1`) where the first item is `Γ0`, and each subsequent item is the result of applying `transform` to the previous item using the corresponding item of `elems`.

# Examples

## Unitful Integration

Package [Unitful](https://github.com/ajkeller34/Unitful.jl) can be
used to specify units e.g. for wavelengths, distances along the beam
axis, beam waist radii, etc. A few examples follow:

!!! note
    Currently, a branch of
    [this Unitful fork](https://github.com/Quantum-Factory/Unitful.jl)
    must be used because the registered package Unitful misses
    some functionality regarding unitful, complex quantities such as
    the beam parameter.

```jldoctest
using ABCDBeamTrace, Unitful
using Unitful: nm, µm, mm, cm, m
f = 125mm
L = 50cm
expander_2x = [ThinLens(f=f), FreeSpace(3f), ThinLens(f=2f)]
system = [expander_2x; FreeSpace(L); reverse(expander_2x)]
inputbeam = GaussianBeam(λ = 532nm, w0 = 500µm)
outputbeam = transform(system, inputbeam)
uconvert(µm, spotradius(outputbeam))

# output

507.11840991014856 μm
```

## Plotting

To plot a system (a vector of [`Element`](@ref)) with a [`GaussianBeam`](@ref)
`beam`, pass both as arguments to `Plots.plot` of the plot
meta-package [Plots](https://github.com/JuliaPlots/Plots.jl). For
exampe, to plot a laser beam's ``1/e^2`` full diameter, use the
following:

```@example
using ABCDBeamTrace, Plots
f = 10e-3 # focal length of 10 mm in SI base units (m)
L = 20e-3 # distance of 20 mm (SI) betweeen sub systems
expander_2x = [ThinLens(f=f), FreeSpace(3f), ThinLens(f=2f)]
system = [expander_2x; FreeSpace(L); reverse(expander_2x)]
beam = GaussianBeam(λ = 532e-9, w0 = 1e-3)
plot(system, beam, size = (800, 150))
savefig("plots-01.svg"); nothing # hide
```

![](plots-01.svg)

Note that specifying an `aspect_ratio` (such as `:none`) is required
if you want to exaggerate the beam width with respect to the distance
along the optical axis and that it is helpful to also specify a custom
`size` of the plot, especially if you do not want to exaggerate the
beam width.

It is possible to overlay several beams; for example, consider this
code generating a plot of three beams of different wavelengths near
the focus of a lens and exaggerates the beam width:

```@example
using ABCDBeamTrace, Plots
using Unitful: nm, µm, mm, m
f = 50mm
L = 1000mm
system = [
    FreeSpace(L), ThinLens(f=f), FreeSpace(0.98f), FreeSpace(0.04f)
]
plot(
    xlims = ( # custom range for x axis
        0.0 + (L+0.98f) / m,
        0.0 + (L+1.02f) / m
    ),
    ylims = (-20.0e-6, 20.0e-6), # custom range for y axis
    size = (800, 500)
)
plot!(system, GaussianBeam(λ = 405nm, w0 = 1mm), label="405 nm")
plot!(system, GaussianBeam(λ = 532nm, w0 = 1mm), label="532 nm")
plot!(system, GaussianBeam(λ = 638nm, w0 = 1mm), label="638 nm")
plot!(
    xlabel = "Distance along Beam Axis [m]",
    ylabel = "Beam 1/e^2 Extent [m]",
    aspect_ratio = :none
)
savefig("plots-02.svg"); nothing # hide
```
![](plots-02.svg)

!!! note
    Some arguments (`xlabel`, `ylabel`, `aspect_ratio`) must be given
    in or after the last invocation of `plot` or `plot!` in order to
    overwrite their default values.

## Mode-Matching with SymPy

If you are interested in using this package together with
[SymPy.jl](https://github.com/JuliaPy/SymPy.jl), check the version of
[python's SymPy](https://www.sympy.org) used behind the scenes:

```jldoctest
using SymPy
SymPy.sympy.__version__

# output

"1.4"
```

If your version is lower, you are using an outdated version: Time to
upgrade! Note that the versions shipped with linux distributions tend
to be rather dated, so it is recommended to use other means of
installation.

Otherwise, the usage of `SymPy` is straight-forward. The following
example uses it for mode-matching a fiber output to a cavity with
known beam waists using two lenses of given focal length and the
constraint of given positions of the waists (i.e. a given total
optical path length):

```@example sympy-modematching
using ABCDBeamTrace, SymPy
@vars z1 positive=true # SymPy variable for position of lens with f1
@vars z2z1 positive=true # SymPy var. for distance between lenses
z2 = z1 + z2z1 # position of lens with f2
L = 0.5 # total length
f1 = 0.050 # focal length for lens at z1
f2 = 0.100 # focal length for lens at z2
system = [
    FreeSpace(z1),
    ThinLens(f = f1),
    FreeSpace(z2 - z1),
    ThinLens(f = f2),
    FreeSpace(L-z2)
]
bin = GaussianBeam(λ = 532e-9, w0 = 5e-6) # laser beam from fiber
bout = transform(system, bin) # laser beam to cavity, calculated
bwant = GaussianBeam(λ = 532e-9, w0 = 20e-6) # to cavity, required
qout = beamparameter(bout)
qwant = beamparameter(bwant)
# numerically solve for SymPy variables z1 and z2z1
z1sol, z2z1sol = sympy.nsolve(
    [real(qout - qwant), imag(qout - qwant)], # expr for root finding
    (z1, z2z1), # SymPy variables
    (L/3, 2L/3) # initial values
)
# output lens positions in units of mm, rounded
(
    round(Int, 1000*convert(Float64, z1(z1 => z1sol, z2z1 => z2z1sol))),
    round(Int, 1000*convert(Float64, z2(z1 => z1sol, z2z1 => z2z1sol)))
)
```

Finally, continuing the code above, the resulting beam profile can be
plotted. The main difficulty lies in performing substitution and type
conversion, handled here by simply constructing the `system` again,
this time with substituted values and using the new name
`finalsystem`:

```@example sympy-modematching
using Plots
# system, but with z1sol, z2z1sol substituted for z1, z2z1
finalsystem = [
    FreeSpace(Float64(z1sol)),
    ThinLens(f = f1),
    FreeSpace(Float64(z2z1sol)),
    ThinLens(f = f2),
    FreeSpace(Float64(L - z2z1sol - z1sol))
]
# plot the entire system
plot(finalsystem, bin)
plot!(
    aspect_ratio = :none, # permit exaggeration of beam width
    size = (800, 300)
)
savefig("plots-03.svg"); nothing # hide
```

![](plots-03.svg)

!!! info "To Do"
    The Gaussian nature of the beam would become obvious
    if inset plots for the region around the two waists
    were provided.

## Symmetric Cavity Design with Symata

[Symata.jl](https://github.com/jlapeyre/Symata.jl) is an entire
language for symbolic mathematics that is somewhat close to the
commercial computer algebra language
[Mathematica](https://www.wolfram.com/mathematica/). It also comes
with bindings for julia which will be used in the following example.
The version of Symata (and its main depedencies) used in the following
is:

```jldoctest
using Symata
@symExpr VersionInfo()

# output

Symata version     0.4.5
Julia version      1.1.0
Python version     3.6.7
┌ Warning: `getindex(o::PyObject, s::Symbol)` is deprecated in favor of dot overloading (`getproperty`) so elements should now be accessed as e.g. `o.s` instead of `o[:s]`.
│   caller = _versioninfo() at kernelstate.jl:33
└ @ Symata ~/.julia/packages/Symata/1oDeR/src/kernelstate.jl:33
SymPy version      1.4
```

Missing from Symata.jl are methods such as `Base.one`, which this
package requires. Hence one must start by retrofitting these, in
addition to calling `symatamath()`:

```@example symata
using Symata
import Base
symatamath()
Base.one(::Mxpr) = @symExpr 1
Base.zero(::Mxpr) = @symExpr 0
Base.inv(m::Mxpr) = mxpr(:Power, m, -1)
Base.sqrt(m::Mxpr) = mxpr(:Power, m, 1//2)
Base.abs(m::Mxpr) = sqrt(m*m)
Base.complex(x::Mxpr, y::Mxpr) = x + im * y
nothing # hide
```

The cavity to be designed is symmetric and esentially a [Fabry–Pérot
Interferometer](https://en.wikipedia.org/wiki/Fabry%E2%80%93P%C3%A9rot_interferometer)
with the exception of having curved rather than flat end mirrors. Both
mirrors shall have the same, fixed radius of curvature `roc` and the
distance `L` between them is a variable. Continuing from the
definition of `Base.one` and `Base.zero` above:

```@example symata
using ABCDBeamTrace
# 100 mm radius of curvature corresponding to 50 mm focal length
roc = 100e-3
# variable distance between mirrors (cavity length)
len = @sym L
if len isa Number # hide
    len = :L # hide
end # hide
# the optical system: a single pass through the cavity
singlepass = [FreeSpace(1.0len), Mirror(roc=roc)]
## a full cavity roundtrip (for this special case)
roundtrip = [singlepass; singlepass]
q = beamparameter(roundtrip)
# define a function to return the resonant mode
# as a Gaussian beam after substitution for length
function modeforlength(actuallength)
    setsymata(len, actuallength) # substitute for L, also in q
    return GaussianBeam(; q = symeval(q), λ = 632.8e-9)
end
println(1e6 * waistradius(modeforlength(75.0e-3)), " µm")
```

!!! note
    For substitution, `setsymata` and `symeval` was used. A nicer
    solution would be using a `@symExpr` involving `Replace` but that
    currently fails because of
    [this issue](https://github.com/jlapeyre/Symata.jl/issues/170).

The resulting on-mirror spotradius, intra-cavity waist radius, and the
intra-cavity Rayleigh range can easily be plotted over the range of
permissible cavity lengths (zero which is effectively co-planar to 0.2
[m] which is concentric).

```@example symata
using Plots
plot(
    z -> 1.0e6 * spotradius(modeforlength(z), none = NaN),
    label = "On-Mirror Spot Radius [µm]",
    xlims = (0.0, 200.0e-3),
    xlabel = "Mirror Distance [m]",
    ylims = (0.0, 400.0),
    ylabel = "1/e^2 Radius [µm] or Rayleigh Range [mm]",
    size = (800, 600)
)
plot!(
    z -> 1.0e6 * waistradius(modeforlength(z), none = NaN),
    label = "Intra-Cavity Waist Radius [µm]"
)
plot!(
    z -> 1.0e3 * rayleighrange(modeforlength(z)),
    label = "Intra-cavity Rayleigh Range [mm]"
)
savefig("plots-04.svg"); nothing # hide
```

![](plots-04.svg)

!!! info "To Do"
    It would be neat to have one or more inset plots of the beam
    profile, e.g. for the confocal case.

## Other Examples
More detailed examples are upcoming.  In particular, this package has been developed with the following tasks in mind:
* Cavity design
* Mode-matching

Jupyter notebooks implementing these tasks exist and are currently being cleaned so that they can be pushed to this repository to serve as examples.

Other interesting examples using this package are welcome!

## To Do
* Add examples of beam and ray tracing function calls
* Add detailed Jupyter notebooks showing examples
* Systematic tests of correctness and consistency

# API Reference

## Optical Elements
Optical elements are defined as types. All elements support extraction
of a sagittal pseudo-element via [`Sag`](@ref) and a tangential (aka
parallel) one via [`Tan`](@ref).

### Elements

```@docs
Element
ElementABCD
FreeSpace
Interface
ThinLens
Mirror
Tan
Sag
```

```@docs
Base.isapprox
```

!!! info "To Do"
    Implement `Base.isapprox` for beams as well.

## Ray Transfer Matrices
```@docs
Base.Matrix
transform
```

## Beam Tracing
### General

Beams are represented by two concrete types. Whilst they are addressed
in different Sections, namely in [Geometric Beams] and [Gaussian
Beams], both sections apply to all kinds of beams. Internally, an
abstract supertype [`AbstractBeam`](@ref) is used.

### Geometric Beams

The following is implemented primarily for geometric optics (using a
[`GeometricBeam`](@ref)) but as Gaussian (laser) optics shares the
same implementation, it is (partially) relevant for Gaussian optics
(using a [`GaussianBeam`](@ref)), too. In particular, the original
main functionality of this package is achieved by [`beamtrace`](@ref).

!!! note
    Unlike suggested by the section headings, the functionality for
    [Gaussian Beams](@ref) largely works for a
    [`GeometricBeam`](@ref) as well, in the simultaneous limits of
    vanishing Rayleigh range, wavelength, and beam waist.

```@docs
GeometricBeam
location
ior
radialpos
slope
discretize
beamtrace
```

### Gaussian Beams

```@docs
GaussianBeam
beamparameter
beamparameterproduct
wavefrontroc
rayleighrange
waistlocation
waistdistance
waistradius
spotradius
spotradiusfunc
```

## Internals

```@docs
AbstractBeam
```

### Accessor Functions
```@docs
η
dz
```

### Miscellanea
```@docs
color
```

## Index
```@index
```
