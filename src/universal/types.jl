# Basic Types
"""
    SCMat{T<:Number}

Structure constants of an algebra. 
"""
struct SCMat{T<:Number} # <: AbstractArray
    d::Int
    # Matrix{SparseVector}
    mat::Matrix{SparseVector{T, Int}}
    syms::Vector{Num} # symbolics basis used to show the elements
end

function SCMat(mat::AbstractMatrix{K}) where K <: AbstractSparseVector{T} where T <: Number
    d = size(mat, 1)
    d == size(mat, 2) || throw(DimensionMismatch("The matrix is not square"))
    syms = eval(Meta.parse("@variables " * join(["x$(map_subscripts(i))" for i in 1:d], " ")))
    return SCMat{T}(d, mat, syms)
end

"""
    LieElem{T<:Number}

An element of the Lie algebra.
"""
struct LieElem{T<:Number}
    scmat::SCMat{T}
    # vector of coefficients of $\{x_1,\cdots,x_n\}$
    element::SparseVector{T}
end

"""
    EnvElem{T<:Number}

An element of the universal enveloping algebra of a Lie algebra.
"""
struct EnvElem{T<:Number}
    scmat::SCMat{T}
    # The decomposition basis of element
    # `keys(element)` is a tuple of indexes
    element::Dict{Tuple, T} # base => coefficience
end
"""zero element"""
EnvElem(scmat::SCMat{T}) where T<:Number = EnvElem{T}(scmat, Dict{Tuple, T}())
"""Initialize by a Lie element"""
EnvElem(x::LieElem{T}) where T<:Number = EnvElem{T}(x.scmat, sparse2dict(x.element))

"""
    AlgebraBySC{T<:Number}

Define an algebra by structure constants.
"""
struct AlgebraBySC{T<:Number}
    scmat::SCMat{T}
    basis::Vector{LieElem{T}}
end

function AlgebraBySC(scmat::SCMat{T}) where T <: Number
    d = dim(scmat)
    basis = [LieElem(scmat, sparsevec([i], ones(T, 1), d)) for i in 1:d]
    return AlgebraBySC{T}(scmat, basis)
end