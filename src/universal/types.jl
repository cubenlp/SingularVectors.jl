
# Basic Types
"""
    SCMat{T<:Number}

Structure constants of a Lie algebra.
"""
struct SCMat # <: AbstractArray
    dim::Int
    # Matrix{SparseVector}
    mat::AbstractMatrix
end

"""
    LieElement{T<:Number}

An element of the Lie algebra.
"""
struct LieElement
    scmat::SCMat
    # vector of coefficients of $\{x_1,\cdots,x_n\}$
    element::SparseVector
end

"""
    EnvElement{T<:Number}

An element of the universal enveloping algebra.
"""
struct EnvElement
    scmat::SCMat
    # The decomposition basis of element
    # `keys(element)` is a tuple of indexes
    element::Dict{Tuple, Int}
    
    """zero element"""
    EnvElement(scmat::SCMat) = new(scmat, Dict{Tuple, Int}())
    EnvElement(scmat::SCMat, element::Dict{Tuple, Int}) = new(scmat, element)
    """Initialize by a Lie element"""
    EnvElement(x::LieElement) = new(x.scmat, sparse2dict(x.element))
end

"""
    LieAlgebra

Data type of Lie algebras.
"""
struct LieAlgebra
    dim::Int
    scmat::SCMat
    basis::Vector{LieElement}
    function LieAlgebra(scmat::SCMat)
        dim = scmat.dim
        basis = [LieElement(scmat, sparsevec([i], [1], dim)) for i in 1:dim]
        return new(dim, scmat, basis)
    end
end