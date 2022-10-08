# Basic Types
"""
    SCMat{T<:Number}

Structure constants of a Lie algebra. 
"""
struct SCMat{T<:Number} # <: AbstractArray
    dim::Int
    # Matrix{SparseVector}
    mat::AbstractMatrix
    syms::Vector{Num} # symbolics basis used to show the elements
end

function SCMat(mat::AbstractMatrix{K}) where K <: AbstractSparseVector{T} where T <: Number
    dim = size(mat, 1)
    dim == size(mat, 2) || throw(DimensionMismatch("The matrix is not square"))
    syms = eval(Meta.parse("@variables " * join(["x$(map_subscripts(i))" for i in 1:dim], " ")))
    return SCMat{T}(dim, mat, syms)
end

"""
    LieElement{T<:Number}

An element of the Lie algebra.
"""
struct LieElement{T<:Number}
    scmat::SCMat{T}
    # vector of coefficients of $\{x_1,\cdots,x_n\}$
    element::SparseVector{T}
end

"""
    EnvElement{T<:Number}

An element of the universal enveloping algebra.
"""
struct EnvElement{T<:Number}
    scmat::SCMat{T}
    # The decomposition basis of element
    # `keys(element)` is a tuple of indexes
    element::Dict{Tuple, T} # base => coefficience
end
"""zero element"""
EnvElement(scmat::SCMat{T}) where T<:Number = EnvElement{T}(scmat, Dict{Tuple, T}())
"""Initialize by a Lie element"""
EnvElement(x::LieElement{T}) where T<:Number = EnvElement{T}(x.scmat, sparse2dict(x.element))

"""
    AlgebraBySC{T<:Number}

Define an algebra by structure constants.
"""
struct AlgebraBySC{T<:Number}
    dim::Int
    scmat::SCMat{T}
    basis::Vector{LieElement{T}}
end

function AlgebraBySC(scmat::SCMat{T}) where T <: Number
    dim = scmat.dim
    basis = [LieElement(scmat, sparsevec([i], ones(T, 1), dim)) for i in 1:dim]
    return AlgebraBySC{T}(dim, scmat, basis)
end