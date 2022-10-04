module SingularVectors

# Data types
export SCMat, LieElement, EnvElement, LieAlgebra, ObjMatchError

# Constructions
export commutative_liealgebra, sl2

# tools
export sparse2dict

include("universalenv.jl")
using .UniversalEnvelope

end
