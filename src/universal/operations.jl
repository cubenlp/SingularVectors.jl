# Operations of algebraic structures

# Basic operations
import Base: +, -, *, ÷, show, getindex, ==, iszero

"""
    ObjMatchError <: Exception

Operations on objects of different algebras.
"""
struct ObjMatchError <: Exception
    msg::String
end

## as vector space
"""
    +(x::LieElement, y::LieElement)

Addition of two elements of the Lie algebra.
"""
function +(x::LieElement, y::LieElement)
    x.scmat === y.scmat || throw(ObjMatchError("Addition of elements of different algebras"))
    return LieElement(x.scmat, x.element + y.element)
end

"""
    -(x::LieElement, y::LieElement)

Subtraction of two elements of the Lie algebra.
"""
function -(x::LieElement, y::LieElement)
    x.scmat === y.scmat || throw(ObjMatchError("Subtraction of elements of different algebras"))
    return LieElement(x.scmat, x.element - y.element)
end

"""
    -(x::LieElement)

Negation of an element of the Lie algebra.
"""
-(x::LieElement) = LieElement(x.scmat, -x.element)

"""
    *(a::Int, x::LieElement)

Multiplication of an element of the Lie algebra by an integer.
"""
function *(a::Int, x::LieElement)
    return LieElement(x.scmat, a * x.element)
end

"""
    ÷(x::LieElement, a::Int)

Division of an element of the Lie algebra by an integer.

**Note: the remainder is discarded without any warning.**
"""
function ÷(x::LieElement, a::Int)
    return LieElement(x.scmat, x.element .÷ a)
end

## as algebra
"""
    *(x::LieElement, y::LieElement)

Multiplication of two elements of the Lie algebra.
"""
function *(x::LieElement, y::LieElement)
    x.scmat === y.scmat || throw(ObjMatchError("Multiplication of elements of different algebras"))
    z = sum(val * val2 * x.scmat[ind, ind2] 
        for (ind, val) in zip(x.element.nzind, x.element.nzval)
            for (ind2, val2) in zip(y.element.nzind, y.element.nzval))
    return LieElement(x.scmat, z)
end

==(x::LieElement, y::LieElement) = x.scmat === y.scmat && x.element == y.element
iszero(x::LieElement) = iszero(x.element)

show(io::IO, alg::LieAlgebra) = print(io, "Lie algebra of dimension $(alg.dim)")
show(io::IO, x::LieElement) = print(io, "Lie element $(x.element)")
show(io::IO, scmat::SCMat) = print(io, scmat.mat)

getindex(scmat::SCMat, i::Int, j::Int) = scmat.mat[i, j]