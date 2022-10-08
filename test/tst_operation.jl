# Test for operations

@testset "Operations for Universal Envelope" begin
    alg = AlgebraBySC(sl2scmat)
    e, h, f = alg.basis
    scmat = alg.scmat
    ee, hh, ff = EnvElement.(alg.basis)
    # basic operations
    @test ee + ff == ff + ee == f + ee == ee + f
    @test ee - hh == -(hh - ee) == -(h - ee) == ee - h
    @test ee + hh + (-hh) == ee
    @test ee * (ff + hh) == ee * ff + ee * hh
    @test (ee + hh) * ff == ee * ff + hh * ff
    @test 2 * (ee + hh) == 2 * ee + 2 * hh == (ee + hh) * 2
    @test unit(ee) * ee == ee == ee * unit(ee)

    # multiplication
    @test ee * ff - ff * ee == hh == EnvElement(e * f)
    @test ee * hh - hh * ee == -2 * ee
    @test hh * ee - ee * hh == 2 * ee
    @test hh * ff - ff * hh == -2 * ff
    @test ff * hh - hh * ff == 2 * ff
    @test ff * ee * hh == ff * (ee * hh)
    @test all(iszero, [xx- xx for xx in [ee, hh, ff]])
    # (ff, ee) = (ee, ff) + ([ff, ee],)
    @test ff * ee == ee * ff + f * e
    # (hh, ff, ee) = (hh, ee, ff) + (hh, [ff, ee])
    #              = (ee, hh, ff) + 2 * (ee, ff) - (hh, hh)
    @test hh * ff * ee == hh * ee * ff + hh * (f * e)
    @test hh * ff * ee == ee * hh * ff + 2 * ee * ff - hh * hh
    @test ee * f - f * ee == hh
    @test all(issortedbypbw, [ee * ff - ff * ee, hh, ee * hh - hh * ee, -2 * ee,
          hh * ff - ff * hh, ff * ee * hh, ee * ff + f * e, hh * ff * ee])
end

@testset "Operations for Lie algebra" begin
    # Test for +, -, *, ÷
    e, h, f = sl2.basis
    @test e + f == f + e
    @test e * f == h
    @test f * e == -h
    @test all(iszero(x * x) for x in [e, h, f])
    @test e * h == -2 * e
    @test h * e == 2 * e
    @test f * h == 2 * f
    @test h * f == -2 * f
    @test all(iszero(x - x) for x in [e, h, f])

    ## test show
    show(stdout, sl2) # AlgebraBySC
    show(stdout, e) # LieElement
    show(stdout, sl2.scmat) # SCMat
    show(stdout, EnvElement(e))
    @test true

    ## commutative algebra
    comalg = commutative_liealgebra(3)
    x, y, z = comalg.basis
    @test all(iszero, [x * y, y * x, x * z, z * x, y * z, z * y,
                       x * x, y * y, z * z])
    @test_throws ObjMatchError x * e
end

@testset "Algebra with parameters" begin
    scmat = SCMat([spzeros(Num, 3) for _ in 1:3, _ in 1:3])
    scmat[1, 2][1], scmat[1, 3][2] = -2, 1
    scmat[2, 1][1], scmat[2, 3][3] = 2, -2
    scmat[3, 1][2], scmat[3, 2][3] = -1, 2
    alg = AlgebraBySC(scmat)
    e, h, f = alg.basis
    ee, hh, ff = EnvElement.(alg.basis)
    @test ee * ff - ff * ee == hh == EnvElement(e * f)
    @variables a b c
    @test a * e * h == -2a * e
    @test h * f * b == -2b * f
    @test a * e * h * f * b == -2a * b * h
end