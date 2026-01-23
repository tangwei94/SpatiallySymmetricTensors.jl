@testset "defined point groups: invalid irreps" begin
    @test_throws ArgumentError SpatiallySymmetricTensors.irrep_chars(C3v(), :BadRep)
    @test_throws ArgumentError SpatiallySymmetricTensors.irrep_chars(C4v(), :BadRep)
    @test_throws ArgumentError SpatiallySymmetricTensors.irrep_chars(C6v(), :BadRep)
    @test_throws ArgumentError SpatiallySymmetricTensors.irrep_chars(C4(), :BadRep)
    @test_throws ArgumentError SpatiallySymmetricTensors.irrep_chars(D2(), :BadRep)
end

@testset "defined point groups: basic perm structure" begin
    # basic sanity: counts and permutation structure per point group/operator
    groups = [
        (C3v(), 4),
        (C4v(), 5),
        (C6v(), 7),
        (C4(), 5),
        (D2(), 5),
    ]
    for (g, nlegs) in groups
        ops = SpatiallySymmetricTensors.group_elements(g)
        for (name, perm) in ops
            @test perm[1] == (1,)
            @test sort(collect(perm[2])) == collect(2:nlegs)
        end
    end
end

@testset "defined point groups: perm algebra" begin
    # check involution for reflections and group order/inverse for rotations
    function compose_perm(p1, p2)
        a = p1[2]
        b = p2[2]
        composed = ntuple(i -> a[b[i] - 1], length(a))
        return (p1[1], composed)
    end

    function is_identity_perm(p, nlegs)
        return p[1] == (1,) && p[2] == ntuple(i -> i + 1, nlegs - 1)
    end

    function perm_power(p, k)
        res = p
        for _ in 2:k
            res = compose_perm(res, p)
        end
        return res
    end

    # reflection two times should be the identity
    for g in (C3v(), C4v(), C6v(), D2())
        nlegs = g isa C6v ? 7 : (g isa C3v ? 4 : 5)
        ops = SpatiallySymmetricTensors.group_elements(g)
        for (name, perm) in ops
            if occursin("σ", String(name))
                @test is_identity_perm(compose_perm(perm, perm), nlegs)
            end
            if occursin("C2", String(name))
                @test is_identity_perm(compose_perm(perm, perm), nlegs)
            end
        end
    end

    # for Cnv, rotation n times should be the identity
    for g in (C4v(), C4())
        ops = SpatiallySymmetricTensors.group_elements(g)
        for (name, perm) in ops
            if occursin("R", String(name)) 
                @test is_identity_perm(perm_power(perm, 4), 5)
            end
        end
    end
    for g in (C6v(),)
        ops = SpatiallySymmetricTensors.group_elements(g)
        for (name, perm) in ops
            if occursin("R", String(name)) 
                @test is_identity_perm(perm_power(perm, 6), 7)
            end
        end
    end
    for g in (C3v(),)
        ops = SpatiallySymmetricTensors.group_elements(g)
        for (name, perm) in ops
            if occursin("R", String(name))
                @test is_identity_perm(perm_power(perm, 3), 4)
            end
        end
    end
end

@testset "defined point groups: representations " begin

    group = C4v()

    for irrep in (:A1, :A2, :B1, :B2, :E)
        rep = SpatiallySymmetricTensors.irrep_rep(group, irrep)
        Id_mat = rep[:Id]
        _is_identity(M) = norm(M - Id_mat) < 1e-12
        @test _is_identity(rep[:R1]^4)
        @test _is_identity(rep[:R3]^4)
        @test _is_identity(rep[:C2]^2)
        @test _is_identity(rep[:σv1]^2)
        @test _is_identity(rep[:σv2]^2)
        @test _is_identity(rep[:σd1]^2)
        @test _is_identity(rep[:σd2]^2)
    end
end
