module SingularVectors

# Data types
export SCMat, LieElem, EnvElem, AlgebraBySC, ObjMatchError

# Constructions
export commutative_liealgebra, sl2, sl2scmat, unit

# tools
export sparse2dict, issortedbypbw

include("universalenv.jl")
using .UniversalEnvelope

end
