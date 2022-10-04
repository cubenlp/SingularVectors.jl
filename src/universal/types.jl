
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
    # NTuple{N, LieElement{T}} where N
    element::Dict{Tuple, Int}
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
        basis = [LieElement(scmat, sparsevec([i], [1])) for i in 1:scmat.dim]
        return new(scmat.dim, scmat, basis)
    end
end