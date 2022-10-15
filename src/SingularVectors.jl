module SingularVectors

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

include("universalenv.jl")
using .UniversalEnvelope

end
