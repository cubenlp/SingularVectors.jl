# Operations of algebraic structures

"""
    ObjMatchError <: Exception

Operations on objects of different algebras.
"""
struct ObjMatchError <: Exception
    msg::String
end

# Basic operations
import Base: +, -, *, ÷, show, getindex, ==, iszero, keys

"""
    getindex(scmat::SCMat, i::Int, j::Int)

Return the entry at position (i, j) in `scmat.mat`.
"""
getindex(scmat::SCMat, i::Int, j::Int) = scmat.mat[i, j]

"""
    sparse2dict(x::SparseVector)

Convert a sparse vector to a dictionary.
"""
function sparse2dict(x::SparseVector)
    return Dict((i,)=>x[i] for i in findall(!iszero, x))
end

## Operations for LieElement
## as vector space
"""
    +(x::LieElement, y::LieElement)

Addition of two elements of the Lie algebra.
"""
function +(x::LieElement{T}, y::LieElement{T}) where T<:Number
    x.scmat === y.scmat || throw(ObjMatchError("Addition of elements of different algebras"))
    return LieElement{T}(x.scmat, x.element + y.element)
end

"""
    -(x::LieElement, y::LieElement)

Subtraction of two elements of the Lie algebra.
"""
function -(x::LieElement{T}, y::LieElement{T}) where T<:Number
    x.scmat === y.scmat || throw(ObjMatchError("Subtraction of elements of different algebras"))
    return LieElement{T}(x.scmat, x.element - y.element)
end

"""
    -(x::LieElement)

Negation of an element of the Lie algebra.
"""
-(x::LieElement{T}) where T<:Number = LieElement{T}(x.scmat, -x.element)

"""
    *(a::Int, x::LieElement)

Multiplication of an element of the Lie algebra by an integer.
"""
function *(a::T, x::LieElement{T}) where T<:Number
    return LieElement{T}(x.scmat, a * x.element)
end
*(x::LieElement, a::T) where T<:Number = a * x

# """
#     ÷(x::LieElement, a::Int)

# Division of an element of the Lie algebra by an integer.

# **Note: the remainder is discarded without any warning.**
# """
# function ÷(x::LieElement{T}, a::Int)
#     return LieElement{T}(x.scmat, x.element .÷ a)
# end

raw"""
    getindex(x::LieElement, i::Int)

Return the coefficient of $x_i$ in `x`.
"""
getindex(x::LieElement, i::Int) = x.element[i]

## as algebra
"""
    *(x::LieElement, y::LieElement)

Multiplication of two elements of the Lie algebra.
"""
function *(x::LieElement{T}, y::LieElement{T}) where T<:Number
    x.scmat === y.scmat || throw(ObjMatchError("Multiplication of elements of different algebras"))
    scmat = x.scmat
    z = sum(x[ind] * y[ind2] * scmat[ind, ind2] 
            for ind in x.element.nzind 
                for ind2 in y.element.nzind)
    return LieElement{T}(scmat, z)
end

==(x::LieElement, y::LieElement) = iszero(x - y)
iszero(x::LieElement) = iszero(x.element)

## Operations for EnvElement
"""
    keys(x::EnvElement)

Return the decomposition basis of `x`.
"""
keys(x::EnvElement) = keys(x.element)

"""
    getindex(x::EnvElement, key::Tuple)

Return the coefficient of `key` in `x`.
"""
getindex(x::EnvElement, key::Tuple) = get(x.element, key, 0)

"""
    zero(x::EnvElement)

Return the zero element of the same algebra as `x`.
"""
zero(x::EnvElement{T}) where T<:Number = EnvElement(x.scmat)

"""
    unit(x::EnvElement)

Return the identity element of the same algebra as `x`.
"""
unit(x::EnvElement{T}) where T<:Number = EnvElement{T}(x.scmat, Dict{Tuple, T}(()=>1))

"""
    +(x::EnvElement, y::EnvElement)

Addition of two elements of the universal enveloping algebra.
"""
function +(x::EnvElement{T}, y::EnvElement{T}) where T<:Number
    x.scmat === y.scmat || throw(ObjMatchError("Addition of elements of different enveloping algebras"))
    ele = Dict{Tuple, T}()
    for k in union(keys(x), keys(y))
        val = x[k] + y[k]
        if !iszero(val)
            ele[k] = val
        end
    end
    return EnvElement{T}(x.scmat, ele)
end

"""
    -(x::EnvElement, y::EnvElement)

Subtraction of two elements of the universal enveloping algebra.
"""
function -(x::EnvElement{T}, y::EnvElement{T}) where T<:Number
    x.scmat === y.scmat || throw(ObjMatchError("Subtraction of elements of different enveloping algebras"))
    ele = Dict{Tuple, T}()
    for k in union(keys(x), keys(y))
        val = x[k] - y[k]
        if !iszero(val)
            ele[k] = val
        end
    end
    return EnvElement(x.scmat, ele)
end

"""
    -(x::EnvElement)

Negation of an element of the universal enveloping algebra.
"""
-(x::EnvElement{T}) where T<:Number = EnvElement{T}(x.scmat, Dict{Tuple, T}(k => -v for (k, v) in x.element))

"""
    *(a::Int, x::EnvElement)

Multiplication of an element of the universal enveloping algebra by an integer.
"""
function *(a::T, x::EnvElement{T}) where T<:Number
    return EnvElement(x.scmat, Dict{Tuple, T}(k => a * v for (k, v) in x.element))
end
*(x::EnvElement{T}, a::T) where T<:Number = a * x

"""
    ÷(x::EnvElement, a::Int)

Division of an element of the universal enveloping algebra by an integer.
"""
# function ÷(x::EnvElement, a::Int)
#     return EnvElement(x.scmat, Dict{Tuple, Int}(k => v ÷ a for (k, v) in x.element))
# end

"""
    ==(x::EnvElement, y::EnvElement)

Test if two elements of the universal enveloping algebra are equal.
"""
==(x::EnvElement, y::EnvElement) = iszero(x - y)

"""
    iszero(x::EnvElement)

Test if an element of the universal enveloping algebra is zero.
"""
iszero(x::EnvElement) = all(iszero, values(x.element))

## **PBW basis**

# @inline tuplejoin(x) = x
@inline tuplejoin(x, y) = (x..., y...)
# @inline tuplejoin(x, y, z...) = tuplejoin(tuplejoin(x, y), z...)

"""
    mult(scmat::SCMat, x::Tuple{LieElement}, y::Tuple{LieElement})

Simplify the product of two basis elements of the universal enveloping algebra.
"""
mult(scmat::SCMat, x::Tuple, y::Tuple) = simplify(scmat, tuplejoin(x, y))

"""
    simplify(x::Tuple{LieElement})

Reduce a sequence of `LieElement` to a standard `EnvElement`.

Basic rule: xᵢxⱼ = xⱼxᵢ + [xᵢ, xⱼ]

Example: (1, 3, 2, 1) => (1, 2, 3, 1) + (1, [3, 2], 1)
"""
function simplify(scmat::SCMat{T}, x::Tuple) where T<:Number
    ind = findfirst(i -> x[i] > x[i+1], 1:length(x)-1)
    isnothing(ind) && return EnvElement{T}(scmat, Dict{Tuple, T}(x=>1))
    # decrease the inversion number by 1
    reducecase = simplify(scmat, (x[1:ind-1]..., x[ind+1], x[ind], x[ind+2:end]...))
    # decrease the length by 1
    ele = scmat[x[ind], x[ind+1]]
    return sum([ele[i] * simplify(scmat, (x[1:ind-1]..., i, x[ind+2:end]...)) 
                for i in ele.nzind]; init=reducecase)
end

"""
    *(x::EnvElement, y::EnvElement)

Multiplication of two elements of the universal enveloping algebra.
"""
function *(x::EnvElement{T}, y::EnvElement{T}) where T<:Number
    x.scmat === y.scmat || throw(ObjMatchError("Multiplication of elements of different enveloping algebras"))
    ele = zero(x)
    for k1 in keys(x), k2 in keys(y)
        coef = x[k1] * y[k2]
        if !iszero(coef)
            ele += coef * mult(x.scmat, k1, k2)
        end
    end
    return ele
end

# Treat Lie algebra as a subalgebra of the Enveloping algebra.
*(x::LieElement{T}, y::EnvElement{T}) where T<:Number = EnvElement(x) * y
*(x::EnvElement{T}, y::LieElement{T}) where T<:Number = x * EnvElement(y)
+(x::EnvElement{T}, y::LieElement{T}) where T<:Number = x + EnvElement(y)
+(x::LieElement{T}, y::EnvElement{T}) where T<:Number = EnvElement(x) + y
-(x::EnvElement{T}, y::LieElement{T}) where T<:Number = x - EnvElement(y)
-(x::LieElement{T}, y::EnvElement{T}) where T<:Number = EnvElement(x) - y

# should be replaced by Symbolic.jl
show(io::IO, alg::AlgebraBySC) = print(io, "Lie algebra of dimension $(alg.dim)")
show(io::IO, x::LieElement) = print(io, sum(x.element .* x.scmat.syms))
show(io::IO, scmat::SCMat) = print(io, "Structure constants of length $(scmat.dim)")
show(io::IO, x::EnvElement) = print(io, "Env element:\n $(x.element)")