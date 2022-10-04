# Universal Enveloping algebra
# outline
## 1. basic operations

module UniversalEnvelope
using SparseArrays

# Data types
export SCMat, LieElement, EnvElement, LieAlgebra

# Constructions
export commutative_liealgebra, sl2

# Data types
include("universal/types.jl")

# Operations
include("universal/operations.jl")

# Construction
include("universal/construct.jl")


end