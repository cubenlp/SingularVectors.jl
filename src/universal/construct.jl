# Construction of Lie algebras

"""
    commutative_liealgebra(::Type{T}, n::Int) where T<:Number

Construct a commutative Lie algebra of dimension `n`.
"""
function commutative_liealgebra(::Type{T}, n::Int) where T<:Number
    scmat = [spzeros(T, n) for _ in 1:n, _ in 1:n]
    return AlgebraBySC(SCMat{T}(n, scmat))
end
commutative_liealgebra(n::Int) = commutative_liealgebra(Int, n)

"""Structure constants of sl2 over Integer fields"""
const sl2scmat = SCMat([spzeros(Int, 3) for _ in 1:3, _ in 1:3])
sl2scmat[1, 2][1], sl2scmat[1, 3][2] = -2, 1
sl2scmat[2, 1][1], sl2scmat[2, 3][3] = 2, -2
sl2scmat[3, 1][2], sl2scmat[3, 2][3] = -1, 2

raw"""
Lie algebra $sl_2$ with chevalley basis $\{e, h, f\}$.
"""
const sl2 = AlgebraBySC(sl2scmat)
