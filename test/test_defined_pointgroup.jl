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
    function is_identity_rep(M, Id_mat)
        return norm(M - Id_mat) < 1e-12
    end

    # reflections / C2 elements square to identity
    # rotations: order checks
    for (group, irreps, rotation_order) in [
        (C4v(), (:A1, :A2, :B1, :B2, :E), 4),
        (C4(), (:A, :B), 4),
        (C6v(), (:A1, :A2, :B1, :B2), 6),
        (C3v(), (:A1, :A2), 3),
    ]
        rep_all = SpatiallySymmetricTensors.group_elements(group)
        for irrep in irreps
            rep = SpatiallySymmetricTensors.irrep_rep(group, irrep)
            Id_mat = rep[:Id]
            for name in keys(rep_all)
                if occursin("σ", String(name)) || occursin("C2", String(name))
                    @test is_identity_rep(rep[name]^2, Id_mat)
                end
                if occursin("R", String(name))
                    @test is_identity_rep(rep[name]^rotation_order, Id_mat)
                end
            end
        end
    end
end
