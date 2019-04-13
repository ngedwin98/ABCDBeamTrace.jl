"""

Compare two vectors of elements using Base.isapprox for each element's
ray matrix (ABCD entries). Does consequently not consider one
discretization of element FreeSpace different from another, or one
realization of an imaging system from another as long as both achieve
(within tolerances) the same imaging.

!!! note
    Optional keyword arguments are not supported.

!!! warning
    This is an experimental feature. Use with caution. A particular
    risky point is that the maximum length scale is used for comparing
    matrix elements B and C—which might lead to near infinite
    tolerances and hence undesired approximate equality for vastly
    different elements.

"""
function Base.isapprox(
    a::Union{Element,Vector{<:Element}},
    b::Union{Element,Vector{<:Element}}
)
    ma = Matrix(a)
    mb = Matrix(b)
    # distance scale: max of abs(B), abs(1/C)
    maxdist1 = max(abs(ma[1,2]), abs(mb[1,2]))
    maxdist2 = max(abs(inv(ma[2,1])), abs(inv(mb[2,1])))
    maxdist = max(maxdist1, maxdist2)
    if isinf(maxdist1)
        maxdist = maxdist2
    elseif isinf(maxdist2)
        maxdist = maxdist1
    end
    # compare A; scale one
    Atol = sqrt(eps(one(ma[1,1])))
    if !isapprox(ma[1,1], mb[1,1]; atol=Atol)
        return false
    end
    # compare B; scale distance
    Btol = sqrt(eps(maxdist))
    if !isapprox(ma[1,2], mb[1,2]; atol=Btol)
        return false
    end
    # compare C; scale inverse distance
    Ctol = sqrt(eps(inv(maxdist)))
    if !isapprox(ma[2,1], mb[2,1]; atol=Ctol)
        return false
    end
    # compare D; scale one (reuse Atol from A)
    if !isapprox(ma[2,2], mb[2,2]; rtol=Atol)
        return false
    end
    return true
end
