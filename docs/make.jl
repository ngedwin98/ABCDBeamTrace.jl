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

# if necessary, repeat the document generation because the example
# code for Symata.jl only works on the second iteration; add a third
# run that won't fail (not in strict mode) to aid debugging
# documentation
for i = 1:3
    bestrict = i < 3
    local docsbuilt
    try
        @info "Starting document build process"
        makedocs(
            sitename = "ABCDBeamTrace",
            format = Documenter.HTML(
                prettyurls = false # put everything in one file
            ),
            # Uncomment the following only when you know what you're
            # doing!  (it will make all doctests pass by overwritting
            # the expected output)
            #doctest = :fix,
            modules = [ABCDBeamTrace],
            strict = bestrict
        )
    catch error
        @error "Documentation build process failed with error" error
    end
    docsbuilt = isfile("docs/build/index.html")
    if docsbuilt
        break
    else
        if i == 1
            @info "Rebuilding docs because first attempt failed"
        elseif i == 2
            @warn string(
                "Turning off strict mode: ",
                "The documentation has a problem!"
            )
            bestrict = false
        end
    end
end

# re-activate parent Project
Pkg.activate(".")

if !isfile("docs/build/index.html")
    error("Building documents failed")
end
