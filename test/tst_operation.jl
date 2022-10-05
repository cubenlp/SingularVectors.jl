# Test for operations
@testset "Operations for Lie algebra" begin
    # Test for +, -, *, ÷
    e, h, f = sl2.basis
    @test e + f == f + e
    @test e * f == h
    @test f * e == -h
    @test all(iszero(x * x) for x in [e, h, f])
    @test e * h == -2 * e
    @test h * e == 2 * e
    @test h * e ÷ 2 == e
    @test f * h == 2 * f
    @test h * f == -2 * f
    @test all(iszero(x - x) for x in [e, h, f])
    # discard the remainder without warning
    @test iszero(e ÷ 2)

    ## test show
    show(stdout, sl2) # AlgebraBySC
    show(stdout, e) # LieElement
    show(stdout, sl2.scmat) # SCMat
    @test true

    ## commutative algebra
    comalg = commutative_liealgebra(3)
    x, y, z = comalg.basis
    @test all(iszero, [x * y, y * x, x * z, z * x, y * z, z * y,
                       x * x, y * y, z * z])
    @test_throws ObjMatchError x * e
end

@testset "Universal Envelope" begin
    e, h, f = sl2.basis
    scmat = sl2.scmat
    ee, hh, ff = EnvElement.(sl2.basis)
    @test ee + ff == ff + ee
    @test ee - hh == -(hh - ee)
    @test ee + hh + (-hh) == ee
    @test ee * (ff + hh) == ee * ff + ee * hh
    @test (ee + hh) * ff == ee * ff + hh * ff

    @test ee * ff - ff * ee == hh
    @test ee * hh - hh * ee == -2 * ee
    @test hh * ee - ee * hh == 2 * ee
    @test hh * ff - ff * hh == -2 * ff
    @test ff * hh - hh * ff == 2 * ff
    @test all(iszero, [xx- xx for xx in [ee, hh, ff]])
    # (ff, ee) = (ee, ff) + ([ff, ee],)
    @test ff * ee == ee * ff + EnvElement(f * e)
end