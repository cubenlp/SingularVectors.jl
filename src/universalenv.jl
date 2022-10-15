# Universal Enveloping algebra
# outline
## 1. basic operations

module UniversalEnvelope
using SparseArrays
using Symbolics: map_subscripts, @variables, Num

# Data types
export SCMat, LieElem, EnvElem, AlgebraBySC
export ObjMatchError

# Constructions
export commutative_liealgebra, sl2, sl2scmat, unit

# tools
export sparse2dict, issortedbypbw

# Data types
include("universal/types.jl")

# Operations
include("universal/operations.jl")

# Construction
include("universal/construct.jl")


end