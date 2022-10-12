# structure of Lie superalgebra
"""
    SupSCMat{T<:Number}

Structure constants of a Lie superalgebra. 
"""
struct SupSCMat{T<:Number} <: AbstractSCMat{T}
    d::Tuple{Int, Int}
    mat::Matrix{SparseVector{T, Int}}
    syms::Vector{Num} # symbolics basis used to show the elements
end

function SupSCMat(mat::AbstractMatrix{K}, even::Int) where K <: AbstractSparseVector{T} where T <: Number
    (dim = size(mat, 1)) == size(mat, 2) || throw(DimensionMismatch("Should be a square matrix"))
    odd = dim - even
    syms = eval(Meta.parse("@variables " * join(
        vcat(["x$(map_subscripts(i))" for i in 1:even], 
             ["y$(map_subscripts(i))" for i in 1:odd]), ' ')))
    return SupSCMat{T}((even, odd), mat, syms)
end

function SupSCMat(mat::AbstractMatrix, even::Int, odd::Int)
    (dim = size(mat, 1)) == even + odd || throw(DimensionMismatch(
        "Sum of the even and odd parts should be equal to the dimension"))
    return SupSCMat(mat, even)
end

"""
    SupElem{T<:Number}

LieElem type of a supalgebra.
"""
struct SupElem{T<:Number} <: AbstractElem{T}
    scmat::SupSCMat{T}
    # vector of coefficients of $\{x_1,\cdots,x_n\}$
    element::SparseVector{T}
end

"""
    SupEnvElem{T<:Number}

An element of the universal enveloping algebra of a Lie supalgebra.
"""
struct SupEnvElem{T<:Number}
    scmat::SupSCMat{T}
    # The decomposition basis of element
    # `keys(element)` is a tuple of indexes
    element::Dict{Tuple, T} # base => coefficience
end
"""zero element"""
SupEnvElem(scmat::SupSCMat{T}) where T<:Number = SupEnvElem{T}(scmat, Dict{Tuple, T}())
"""Initialize by a Lie element"""
SupEnvElem(x::SupElem{T}) where T<:Number = SupEnvElem{T}(x.scmat, sparse2dict(x.element))
Base.convert(::Type{SupEnvElem{T}}, x::SupElem{T}) where T = SupEnvElem(x)

"""
    SupAlgebraBySC{T<:Number}

Define a Lie superalgebra by structure constants.
"""
struct SupAlgBySC{T<:Number}
    scmat::SupSCMat{T}
    evenbasis::Vector{SupElem{T}}
    oddbasis::Vector{SupElem{T}}
end

function SupAlgBySC(scmat::SupSCMat{T}) where T <: Number
    even, odd = scmat.d
    dim = even + odd
    evenbasis = [SupElem(scmat, sparsevec([i], ones(T, 1), dim)) for i in 1:even]
    oddbasis = [SupElem(scmat, sparsevec([i], ones(T, 1), dim)) for i in (even+1):dim]
    return SupAlgBySC{T}(scmat, evenbasis, oddbasis)
end

## operations
"""
    simplify(x::Tuple{LieElem})

Reduce a sequence of `LieElem` to a standard `EnvElem`.

Basic rule: xᵢxⱼ = (-1)ⁱʲ⋅xⱼxᵢ + [xᵢ, xⱼ]

Example: (1, 3, 2, 1) => (-1)^{|2|⋅|3|}(1, 2, 3, 1) + (1, [3, 2], 1)
"""
function simplify(scmat::SupSCMat{T}, x::Tuple) where T<:Number
    ind = findfirst(i -> x[i] > x[i+1], 1:length(x)-1)
    isnothing(ind) && return SupEnvElem(scmat, Dict{Tuple, T}(x=>1))
    # decrease the inversion number by 1
    reducecase = simplify(scmat, (x[1:ind-1]..., x[ind+1], x[ind], x[ind+2:end]...))
    even, _ = scmat.d
    if x[ind] > even && x[ind+1] > even # both odd
        times!(reducecase, -1)
    end
    # decrease the length by 1
    ele = scmat[x[ind], x[ind+1]]
    for i in ele.nzind
        add!(reducecase, times!(simplify(scmat, (x[1:ind-1]..., i, x[ind+2:end]...)), ele[i]))
    end
    return reducecase
end

# should be replaced by Symbolic.jl
function show(io::IO, alg::SupAlgBySC)
    txt = "Lie superalgebra of dimension $(dim(alg)) with basis: " * join(string.(alg.basis), ", ")
    print(io,  txt)
end

dim(alg::SupAlgBySC) = alg.scmat.d
eltype(alg::SupAlgBySC{T}) where T = T