import HaloArrays: HaloArray,
                   nhalo,
                   comm,
                   reqs,
                   origin,
                   LEFT,
                   RIGHT

using Test
using MPI

MPI.Init()

@testset "constructor and accessors              " begin
    # should work
    a = HaloArray{Float64}(MPI.COMM_WORLD,
                           (   1,    2,    1),
                           (true, true, true),
                           (  64,   64,   64),
                           (   1,    2,    3))
    @test size(a)  == (64, 64, 64)
    @test nhalo(a) == (1, 2, 3)
    @test size(parent(a)) == (66, 68, 70)
    @test origin(a) == (2, 3, 4)
    @test length(reqs(a)) == 3
    @test sum(length, reqs(a)) == 12
    @test length(reqs(a, 1)) == 4
    @test length(reqs(a, 2)) == 4
    @test length(reqs(a, 3)) == 4
    @test all(dim -> reqs(a)[dim] === reqs(a, dim), 1:3)
    @test all(r -> eltype(r) === MPI.Request, reqs(a))

    # these should error
    # invalid topology
    @test_throws ArgumentError HaloArray{Float64}(MPI.COMM_WORLD,
                                                  (   1,    2,    3),
                                                  (true, true, true),
                                                  (  64,   64,   64),
                                                  (   1,    2,    3))

    # no halo points along periodic direction
    @test_throws ArgumentError HaloArray{Float64}(MPI.COMM_WORLD,
                                                  (   1,    2,    1),
                                                  (true, true, true),
                                                  (  64,   64,   64),
                                                  (   0,    1,    1))

    # no halo points along direction with more than one proc
    @test_throws ArgumentError HaloArray{Float64}(MPI.COMM_WORLD,
                                                  (   1,    2,    1),
                                                  (true, true, true),
                                                  (  64,   64,   64),
                                                  (   1,    0,    1))

    # negative halo region
    @test_throws ArgumentError HaloArray{Float64}(MPI.COMM_WORLD,
                                                  (   1,    2,    1),
                                                  (true, true, true),
                                                  (  64,   64,   64),
                                                  (  -1,    1,    1))
end

@testset "similar and copy                       " begin

    # should work
    a = HaloArray{Float64}(MPI.COMM_WORLD,
                           (   1,    2,    1),
                           (true, true, true),
                           (   1,    1,    1),
                           (   1,    1,    1))

    a[1, 1, 1] = 1
    @test a[1, 1, 1] == 1

    b = copy(a)
    @test b[1, 1, 1] == 1
    @test nhalo(b) == (1, 1, 1)

    c = similar(a)
    @test nhalo(c) == (1, 1, 1)

    parent(a) .= reshape(1:length(parent(a)), size(parent(a)))
    b = copy(a)
    @test parent(b) == parent(a)
end

@testset "similar preserves communication options" begin
    a = HaloArray{Float64}(MPI.COMM_WORLD,
                           (   1,    2,    1),
                           (true, true, true),
                           (   1,    1,    1),
                           (   1,    1,    1); economic=false, safe=false)

    b = similar(a)
    @test b.economic == false
    @test b.safe == false
    @test typeof(reqs(b)) == typeof(reqs(a))
    @test all(r -> eltype(r) === MPI.UnsafeRequest, reqs(b))
    @test comm(b) == comm(a)
end


@testset "halo point indexing                    " begin
    # locally the actual data is a 3x3x3 array and we
    # index it with indices from 0 to 2
    a = HaloArray{Float64}(MPI.COMM_WORLD,
                           (   1,    2,    1),
                           (true, true, true),
                           (   1,    1,    1),
                           (   1,    1,    1))

    for k = 0:2, j = 0:2, i = 0:2
        a[i, j, k] = i*j*k
        @test a[i, j, k] == i*j*k
    end
    @test_throws BoundsError a[-1, 0, 1]
end

MPI.Finalize()
