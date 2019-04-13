import Pkg;
# use a separate Project for building docs; see
# https://discourse.julialang.org/t/psa-use-a-project-for-building-your-docs/14974
Pkg.activate("docs")

# to find parent package when executing this file using
# `include("docs/make.jl")`
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end

# to find parent package from examples, doctests, etc. (executed in
# working directory docs/build/)
if !("../.." in LOAD_PATH)
    push!(LOAD_PATH, "../..")
end

using Documenter, ABCDBeamTrace

makedocs(
    sitename = "ABCDBeamTrace",
    format = Documenter.HTML(
        prettyurls = false # put everything in one file
    ),
    modules = [ABCDBeamTrace]
)

# re-activate parent Project
Pkg.activate(".")
