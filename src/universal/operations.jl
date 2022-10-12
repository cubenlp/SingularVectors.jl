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
function sparse2dict(x::AbstractSparseVector{T}) where T<:Number
    return Dict{Tuple, T}((i,)=>x[i] for i in findall(!iszero, x))
end

## Operations for LieElem
## as vector space
"""
    +(x::LieElem, y::LieElem)

Addition of two elements of the Lie algebra.
"""
function +(x::LieElem{T}, y::LieElem{T}) where T<:Number
    x.scmat === y.scmat || throw(ObjMatchError("Addition of elements of different algebras"))
    return LieElem{T}(x.scmat, x.element + y.element)
end

"""
    -(x::LieElem, y::LieElem)

Subtraction of two elements of the Lie algebra.
"""
function -(x::LieElem{T}, y::LieElem{T}) where T<:Number
    x.scmat === y.scmat || throw(ObjMatchError("Subtraction of elements of different algebras"))
    return LieElem{T}(x.scmat, x.element - y.element)
end

"""
    -(x::LieElem)

Negation of an element of the Lie algebra.
"""
-(x::LieElem{T}) where T<:Number = LieElem{T}(x.scmat, -x.element)

"""
    *(a::Int, x::LieElem)

Multiplication of an element of the Lie algebra by an integer.
"""
function *(a::T, x::LieElem{T}) where T<:Number
    return LieElem{T}(x.scmat, a * x.element)
end
*(x::LieElem, a::T) where T<:Number = a * x

raw"""
    getindex(x::LieElem, i::Int)

Return the coefficient of $x_i$ in `x`.
"""
getindex(x::LieElem, i::Int) = x.element[i]

## as algebra
"""
    *(x::LieElem, y::LieElem)

Multiplication of two elements of the Lie algebra.
"""
function *(x::LieElem{T}, y::LieElem{T}) where T<:Number
    x.scmat === y.scmat || throw(ObjMatchError("Multiplication of elements of different algebras"))
    scmat = x.scmat
    z = spzeros(T, scmat.d)
    for i in x.element.nzind, j in y.element.nzind
        val = x[i] * y[j]
        iszero(val) && continue
        z += val * scmat[i, j]
    end
    return LieElem(scmat, z)
end

==(x::LieElem, y::LieElem) = iszero(x - y)
iszero(x::LieElem) = iszero(x.element)

## Operations for EnvElem
"""
    keys(x::EnvElem)

Return the decomposition basis of `x`.
"""
keys(x::EnvElem) = keys(x.element)

"""
    getindex(x::EnvElem, key::Tuple)

Return the coefficient of `key` in `x`.
"""
getindex(x::EnvElem, key::Tuple) = get(x.element, key, 0)

"""
    zero(x::EnvElem)

Return the zero element of the same algebra as `x`.
"""
zero(x::EnvElem{T}) where T<:Number = EnvElem(x.scmat)

"""
    unit(x::EnvElem)

Return the identity element of the same algebra as `x`.
"""
unit(x::EnvElem{T}) where T<:Number = EnvElem{T}(x.scmat, Dict{Tuple, T}(()=>1))

"""
    +(x::EnvElem, y::EnvElem)

Addition of two elements of the universal enveloping algebra.
"""
function +(x::EnvElem{T}, y::EnvElem{T}) where T<:Number
    x.scmat === y.scmat || throw(ObjMatchError("Addition of elements of different enveloping algebras"))
    ele = Dict{Tuple, T}()
    for k in union(keys(x), keys(y))
        val = x[k] + y[k]
        if !iszero(val)
            ele[k] = val
        end
    end
    return EnvElem{T}(x.scmat, ele)
end

"""Addition without creating new memory"""
function add!(x::EnvElem, y::EnvElem)
    for k in keys(y)
        val = x[k] + y[k]
        if iszero(val)
            delete!(x.element, k)
        else
            x.element[k] = val
        end
    end
    return x
end

"""
    -(x::EnvElem, y::EnvElem)

Subtraction of two elements of the universal enveloping algebra.
"""
function -(x::EnvElem{T}, y::EnvElem{T}) where T<:Number
    x.scmat === y.scmat || throw(ObjMatchError("Subtraction of elements of different enveloping algebras"))
    ele = Dict{Tuple, T}()
    for k in union(keys(x), keys(y))
        val = x[k] - y[k]
        if !iszero(val)
            ele[k] = val
        end
    end
    return EnvElem(x.scmat, ele)
end

"""
    -(x::EnvElem)

Negation of an element of the universal enveloping algebra.
"""
-(x::EnvElem{T}) where T<:Number = EnvElem{T}(x.scmat, Dict{Tuple, T}(k => -v for (k, v) in x.element))

"""
    *(a::Int, x::EnvElem)

Multiplication of an element of the universal enveloping algebra by an integer.
"""
function *(a::T, x::EnvElem{T}) where T<:Number
    return EnvElem(x.scmat, Dict{Tuple, T}(k => a * v for (k, v) in x.element))
end
*(x::EnvElem{T}, a::T) where T<:Number = a * x

"""multiplication without creating new memory"""
function times!(x::EnvElem{T}, a::T) where T<:Number
    for k in keys(x.element)
        x.element[k] *= a
    end
    return x
end    

# """
#     ÷(x::EnvElem, a::Int)

# Division of an element of the universal enveloping algebra by an integer.
# """
# function ÷(x::EnvElem, a::Int)
#     return EnvElem(x.scmat, Dict{Tuple, Int}(k => v ÷ a for (k, v) in x.element))
# end

"""
    ==(x::EnvElem, y::EnvElem)

Test if two elements of the universal enveloping algebra are equal.
"""
==(x::EnvElem, y::EnvElem) = iszero(x - y)

"""
    iszero(x::EnvElem)

Test if an element of the universal enveloping algebra is zero.
"""
iszero(x::EnvElem) = all(iszero, values(x.element))

## **PBW basis**

# @inline tuplejoin(x) = x
@inline tuplejoin(x, y) = (x..., y...)
# @inline tuplejoin(x, y, z...) = tuplejoin(tuplejoin(x, y), z...)

"""
    mult(scmat::SCMat, x::Tuple{LieElem}, y::Tuple{LieElem})

Simplify the product of two basis elements of the universal enveloping algebra.
"""
mult(scmat::SCMat, x::Tuple, y::Tuple) = simplify(scmat, tuplejoin(x, y))

"""
    simplify(x::Tuple{LieElem})

Reduce a sequence of `LieElem` to a standard `EnvElem`.

Basic rule: xᵢxⱼ = xⱼxᵢ + [xᵢ, xⱼ]

Example: (1, 3, 2, 1) => (1, 2, 3, 1) + (1, [3, 2], 1)
"""
function simplify(scmat::SCMat{T}, x::Tuple) where T<:Number
    ind = findfirst(i -> x[i] > x[i+1], 1:length(x)-1)
    isnothing(ind) && return EnvElem(scmat, Dict{Tuple, T}(x=>1))
    # decrease the inversion number by 1
    reducecase = simplify(scmat, (x[1:ind-1]..., x[ind+1], x[ind], x[ind+2:end]...))
    # decrease the length by 1
    ele = scmat[x[ind], x[ind+1]]
    for i in ele.nzind
        add!(reducecase, times!(simplify(scmat, (x[1:ind-1]..., i, x[ind+2:end]...)), ele[i]))
    end
    return reducecase
end

"""
    *(x::EnvElem, y::EnvElem)

Multiplication of two elements of the universal enveloping algebra.
"""
function *(x::EnvElem{T}, y::EnvElem{T}) where T<:Number
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

"""
    issortedbypbw(x::EnvElem)

Test if the key values of the `EnvElem` is sorted by PBW basis.
"""
issortedbypbw(x::EnvElem) = all(issorted, keys(x.element))

# Treat Lie algebra as a subalgebra of the Enveloping algebra.
*(x::LieElem{T}, y::EnvElem{T}) where T<:Number = EnvElem(x) * y
*(x::EnvElem{T}, y::LieElem{T}) where T<:Number = x * EnvElem(y)
+(x::EnvElem{T}, y::LieElem{T}) where T<:Number = x + EnvElem(y)
+(x::LieElem{T}, y::EnvElem{T}) where T<:Number = EnvElem(x) + y
-(x::EnvElem{T}, y::LieElem{T}) where T<:Number = x - EnvElem(y)
-(x::LieElem{T}, y::EnvElem{T}) where T<:Number = EnvElem(x) - y

# element type
eltype(::EnvElem{T}) where T = T
eltype(::LieElem{T}) where T = T
eltype(::SCMat{T}) where T = T
eltype(::AlgebraBySC{T}) where T = T

# dimension
dim(x::EnvElem) = x.scmat.d
dim(x::LieElem) = x.scmat.d
dim(x::SCMat) = x.d
dim(x::AlgebraBySC) = x.scmat.d

# should be replaced by Symbolic.jl
function show(io::IO, alg::AlgebraBySC)
    txt = "Lie algebra of dimension $(dim(alg)) with basis: " * join(string.(alg.basis), ", ")
    print(io,  txt)
end

show(io::IO, x::LieElem) = print(io, sum(x.element .* x.scmat.syms))

show(io::IO, scmat::SCMat) = print(io, "Structure constants of length $(scmat.d)")

"""
    show(io::IO, x::EnvElem)

Show an element of the universal enveloping algebra.

Note: the output is sorted by dictionary order, not by PBW basis. However, the 
user should know that it represents the element defined by PBW ordering.

Also, one can use `issortedbypbw` to test if the basis is a standard PBW basis.
"""
function show(io::IO, x::EnvElem)
    # warning: nonstandard ordering !
    dict, syms = x.element, x.scmat.syms
    key2sym(key) = prod(syms[i] for i in key)
    print(io, sum(val * key2sym(k) for (k, val) in dict))
end