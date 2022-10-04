# Construction of Lie algebras

"""
    commutative_liealgebra(n::Int)

Construct a commutative Lie algebra of dimension `n`.
"""
function commutative_liealgebra(n::Int)
    scmat = [spzeros(Int, n) for _ in 1:n, _ in 1:n]
    return LieAlgebra(SCMat(n, scmat))
end

const scmat = [spzeros(Int, 3) for _ in 1:3, _ in 1:3]
scmat[1, 2], scmat[1, 3] = sparsevec([1], [-2], 3), sparsevec([2], [1], 3)
scmat[2, 1], scmat[2, 3] = sparsevec([1], [2], 3), sparsevec([3], [-2])
scmat[3, 1], scmat[3, 2] = sparsevec([2], [-1], 3), sparsevec([3], [2])

raw"""
Lie algebra $sl_2$ with chevalley basis $\{e, h, f\}$.
"""
const sl2 = LieAlgebra(SCMat(3, scmat))
