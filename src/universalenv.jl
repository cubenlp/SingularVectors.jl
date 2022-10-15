# Universal Enveloping algebra
# outline
## 1. basic operations

module UniversalEnvelope
using SparseArrays
using Symbolics: map_subscripts, @variables, Num

export 
    # element types
    AbstractElem, LieElem, SupElem,
    AbstractEnvElem, LieEnvElem, SupEnvElem,
    # structure constants
    AbstractSCMat, SCMat, SupSCMat,
    # construction of algebras
    AlgebraBySC, commutative_liealgebra,
    # tools
    sparse2dict, issortedbypbw,
    # others
    ObjMatchError, sl2, sl2scmat, unit

# Data types
include("universal/types.jl")

# Operations
include("universal/operations.jl")

# Construction
include("universal/construct.jl")

include("universal/suptypes.jl")
end