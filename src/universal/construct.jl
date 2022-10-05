# Construction of Lie algebras

"""
    commutative_liealgebra(n::Int)

Construct a commutative Lie algebra of dimension `n`.
"""
function commutative_liealgebra(n::Int)
    scmat = [spzeros(Int, n) for _ in 1:n, _ in 1:n]
    return AlgebraBySC(SCMat(n, scmat))
end

const sl2scmat = [spzeros(Int, 3) for _ in 1:3, _ in 1:3]
sl2scmat[1, 2][1], sl2scmat[1, 3][2] = -2, 1
sl2scmat[2, 1][1], sl2scmat[2, 3][3] = 2, -2
sl2scmat[3, 1][2], sl2scmat[3, 2][3] = -1, 2

raw"""
Lie algebra $sl_2$ with chevalley basis $\{e, h, f\}$.
"""
const sl2 = AlgebraBySC(SCMat(3, sl2scmat))
