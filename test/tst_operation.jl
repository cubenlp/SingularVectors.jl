# Test for operations
@testset "Basic operations" begin
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

    ## test show
    show(stdout, sl2)
    show(stdout, sl2.basis)
    show(stdout, sl2.basis[1])
    show(stdout, sl2.scmat)
    @test true

    ## commutative algebra
    comalg = commutative_liealgebra(3)
    x, y, z = comalg.basis
    @test all(iszero, [x * y, y * x, x * z, z * x, y * z, z * y,
                       x * x, y * y, z * z])
    @test_throws ObjMatchError x * e
end