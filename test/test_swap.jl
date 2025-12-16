import HaloArrays: HaloArray,
                   nhalo,
                   comm,
                   reqs,
                   source_dest_ranks,
                   origin,
                   LEFT,
                   RIGHT,
                   haloswap!,
                   swapregions

using Test
using MPI

MPI.Init()

@testset "source_dest_ranks                      " begin
    
    # test in 1D
    # PROC 3 | PROC 0 | PROC 1 | PROC 2 | PROC 3 | PROC 0
    a = HaloArray{Float64}(MPI.COMM_WORLD,
                           (4,), (true,), (2,), (1,))

    if MPI.Comm_rank(comm(a)) == 3
        @test source_dest_ranks(a, (LEFT,))  == (0, 2)
        @test source_dest_ranks(a, (RIGHT,)) == (2, 0)
    end
                           
    if MPI.Comm_rank(comm(a)) == 2
        @test source_dest_ranks(a, (LEFT,))  == (3, 1)
        @test source_dest_ranks(a, (RIGHT,)) == (1, 3)
    end

    if MPI.Comm_rank(comm(a)) == 1
        @test source_dest_ranks(a, (LEFT,))  == (2, 0)
        @test source_dest_ranks(a, (RIGHT,)) == (0, 2)
    end

    if MPI.Comm_rank(comm(a)) == 0
        @test source_dest_ranks(a, (LEFT,))  == (1, 3)
        @test source_dest_ranks(a, (RIGHT,)) == (3, 1)
    end

    # test in 1D
    # PROC 0 | PROC 1 | PROC 2 | PROC 3
    a = HaloArray{Float64}(MPI.COMM_WORLD,
                           (4,), (false,), (2,), (1,))
    
    if MPI.Comm_rank(comm(a)) == 3
        @test source_dest_ranks(a, (LEFT,))  == (MPI.PROC_NULL, 2)
        @test source_dest_ranks(a, (RIGHT,)) == (2, MPI.PROC_NULL)
    end
end

@testset "economic = true                        " begin
    @testset "swap 1                                 " begin
        # before the swap
        # + - - - + - - - +
        # | 0 0 0 | 1 1 1 |
        # | 0 0 0 | 1 1 1 |
        # | 0 0 0 | 1 1 1 |
        # + - - - + - - - +
        # | 2 2 2 | 3 3 3 |
        # | 2 2 2 | 3 3 3 |
        # | 2 2 2 | 3 3 3 |
        # + - - - + - - - +
        # after the swap
        # + - - - + - - - +
        # | 0 2 0 | 1 3 1 |
        # | 1 0 1 | 0 1 0 |
        # | 0 2 0 | 1 3 1 |
        # + - - - + - - - +
        # | 2 0 2 | 3 1 3 |
        # | 3 2 3 | 2 3 2 |
        # | 2 0 2 | 3 1 3 |
        # + - - - + - - - +
        rank = MPI.Comm_rank(MPI.COMM_WORLD)
        a = HaloArray{Float64}(MPI.COMM_WORLD, 
                               (2,    2), 
                               (true, true),
                               (1,    1),
                               (1,    1))
        parent(a) .= rank
        haloswap!(a, rand([true, false]))
        MPI.Waitall(reqs(a))
        rank == 0 && @test parent(a) == [0 2 0; 
                                         1 0 1; 
                                         0 2 0]
        rank == 1 && @test parent(a) == [1 3 1; 
                                         0 1 0; 
                                         1 3 1]
        rank == 2 && @test parent(a) == [2 0 2; 
                                         3 2 3; 
                                         2 0 2]
        rank == 3 && @test parent(a) == [3 1 3; 
                                         2 3 2; 
                                         3 1 3]
    end

    @testset "swap 2                                 " begin
        # before the swap
        # | 0 0 0 | 1 1 1 |
        # | 0 0 0 | 1 1 1 |
        # | 0 0 0 | 1 1 1 |
        # + - - - + - - - +
        # | 2 2 2 | 3 3 3 |
        # | 2 2 2 | 3 3 3 |
        # | 2 2 2 | 3 3 3 |

        # after the swap
        # | 0 0 0 | 1 1 1 |
        # | 1 0 1 | 0 1 0 |
        # | 0 2 0 | 1 3 1 |
        # + - - - + - - - +
        # | 2 0 2 | 3 1 3 |
        # | 3 2 3 | 2 3 2 |
        # | 2 2 2 | 3 3 3 |
        rank = MPI.Comm_rank(MPI.COMM_WORLD)
        a = HaloArray{Float64}(MPI.COMM_WORLD, 
                               (2,     2), 
                               (false, true),
                               (1,     1),
                               (1,     1))
        parent(a) .= rank
        haloswap!(a, rand([true, false]))
        MPI.Waitall(reqs(a))
        rank == 0 && @test parent(a) == [0 0 0; 
                                         1 0 1; 
                                         0 2 0]
        rank == 1 && @test parent(a) == [1 1 1; 
                                         0 1 0; 
                                         1 3 1]
        rank == 2 && @test parent(a) == [2 0 2; 
                                         3 2 3; 
                                         2 2 2]
        rank == 3 && @test parent(a) == [3 1 3; 
                                         2 3 2; 
                                         3 3 3]
    end

    @testset "swap 3                                 " begin
        # before the swap
        # 0 0 0 | 1 1 1
        # 0 0 0 | 1 1 1
        # 0 0 0 | 1 1 1
        # - - - + - - -
        # 2 2 2 | 3 3 3
        # 2 2 2 | 3 3 3
        # 2 2 2 | 3 3 3

        # after the swap
        # 0 0 0 | 1 1 1
        # 0 0 1 | 0 1 1
        # 0 2 0 | 1 3 1
        # - - - + - - -
        # 2 0 2 | 3 1 3
        # 2 2 3 | 2 3 3
        # 2 2 2 | 3 3 3
        rank = MPI.Comm_rank(MPI.COMM_WORLD)
        a = HaloArray{Float64}(MPI.COMM_WORLD,
                               (2,     2),
                               (false, false),
                               (1,     1),
                               (1,     1))
        parent(a) .= rank
        haloswap!(a, rand([true, false]))
        MPI.Waitall(reqs(a))
        rank == 0 && @test parent(a) == [0 0 0; 
                                         0 0 1; 
                                         0 2 0]
        rank == 1 && @test parent(a) == [1 1 1; 
                                         0 1 1; 
                                         1 3 1]
        rank == 2 && @test parent(a) == [2 0 2; 
                                         2 2 3; 
                                         2 2 2]
        rank == 3 && @test parent(a) == [3 1 3; 
                                         2 3 3; 
                                         3 3 3]
    end

    @testset "swap 4                                 " begin
        # before the swap
        # 0 0 0 | 1 1 1 | 2 2 2 | 3 3 3
        # 0 0 0 | 1 1 1 | 2 2 2 | 3 3 3
        # 0 0 0 | 1 1 1 | 2 2 2 | 3 3 3

        # after the swap
        # 0 0 0 | 1 1 1 | 2 2 2 | 3 3 3
        # 0 0 1 | 0 1 2 | 1 2 3 | 2 3 3
        # 0 0 0 | 1 1 1 | 2 2 2 | 3 3 3
        rank = MPI.Comm_rank(MPI.COMM_WORLD)
        a = HaloArray{Float64}(MPI.COMM_WORLD,
                               (1,     4),
                               (false, false),
                               (1,     1),
                               (1,     1))
        parent(a) .= rank
        haloswap!(a, rand([true, false]))
        MPI.Waitall(reqs(a))
        rank == 0 && @test parent(a) == [0 0 0; 
                                         0 0 1; 
                                         0 0 0]
        rank == 1 && @test parent(a) == [1 1 1; 
                                         0 1 2; 
                                         1 1 1]
        rank == 2 && @test parent(a) == [2 2 2; 
                                         1 2 3; 
                                         2 2 2]
        rank == 3 && @test parent(a) == [3 3 3; 
                                         2 3 3; 
                                         3 3 3]
    end

    @testset "swap 5                                 " begin
        # before the swap
        # 0 0 0 | 1 1 1 | 2 2 2 | 3 3 3
        # 0 0 0 | 1 1 1 | 2 2 2 | 3 3 3

        # after the swap
        # 0 0 1 | 0 1 2 | 1 2 3 | 2 3 3
        # 0 0 1 | 0 1 2 | 1 2 3 | 2 3 3
        rank = MPI.Comm_rank(MPI.COMM_WORLD)
        a = HaloArray{Float64}(MPI.COMM_WORLD,
                               (1,     4),
                               (false, false),
                               (2,     1),
                               (0,     1))
        parent(a) .= rank
        haloswap!(a, rand([true, false]))
        MPI.Waitall(reqs(a))
        rank == 0 && @test parent(a) == [0 0 1;
                                         0 0 1]
        rank == 1 && @test parent(a) == [0 1 2;
                                         0 1 2]
        rank == 2 && @test parent(a) == [1 2 3;
                                         1 2 3]
        rank == 3 && @test parent(a) == [2 3 3;
                                         2 3 3]
    end

    @testset "swap 6                                 " begin
        # before the swap
        # -------- + -------- + -------- + --------
        #  0  0  0 |  0  0  0 |  0  0  0 |  0  0  0
        #  0  1  0 |  0  3  0 |  0  5  0 |  0  7  0
        #  0  2  0 |  0  4  0 |  0  6  0 |  0  8  0
        #  0  0  0 |  0  0  0 |  0  0  0 |  0  0  0
        # -------- + -------- + -------- + --------

        # after the swap
        # -------- + -------- + -------- + --------
        #  0  2  0 |  0  4  0 |  0  6  0 |  0  8  0
        #  0  1  3 |  1  3  5 |  3  5  7 |  5  7  0
        #  0  2  4 |  2  4  6 |  4  6  8 |  6  8  0
        #  0  1  0 |  0  3  0 |  0  5  0 |  0  7  0
        # -------- + -------- + -------- + --------

        rank = MPI.Comm_rank(MPI.COMM_WORLD)
        a = HaloArray{Float64}(MPI.COMM_WORLD,
                               (1,     4),
                               (true,  false),
                               (2,     1),
                               (1,     1); economic=true)
        a .= [2*rank + 1; 
              2*rank + 2]

        haloswap!(a, rand([true, false]))
        MPI.Waitall(reqs(a))
        rank == 0 && @test parent(a) == [0  2  0;
                                         0  1  3;
                                         0  2  4;
                                         0  1  0]
                                         
        rank == 1 && @test parent(a) == [0  4  0;
                                         1  3  5;
                                         2  4  6;
                                         0  3  0]
                                         
        rank == 2 && @test parent(a) == [0  6  0;
                                         3  5  7;
                                         4  6  8;
                                         0  5  0]
                                         
        rank == 3 && @test parent(a) == [0  8  0;
                                         5  7  0;
                                         6  8  0;
                                         0  7  0]
    end
end

@testset "economic = false                       " begin
    @testset "swap 1                                 " begin
        # before the swap
        # + - - - + - - - +
        # | 0 0 0 | 1 1 1 |
        # | 0 0 0 | 1 1 1 |
        # | 0 0 0 | 1 1 1 |
        # + - - - + - - - +
        # | 2 2 2 | 3 3 3 |
        # | 2 2 2 | 3 3 3 |
        # | 2 2 2 | 3 3 3 |
        # + - - - + - - - +
        # after the swap
        # + - - - + - - - +
        # | 3 2 3 | 2 3 2 |
        # | 1 0 1 | 0 1 0 |
        # | 3 2 3 | 2 3 2 |
        # + - - - + - - - +
        # | 1 0 1 | 0 1 0 |
        # | 3 2 3 | 2 3 2 |
        # | 1 0 1 | 0 1 0 |
        # + - - - + - - - +
        rank = MPI.Comm_rank(MPI.COMM_WORLD)
        a = HaloArray{Float64}(MPI.COMM_WORLD, 
                               (2,    2), 
                               (true, true),
                               (1,    1),
                               (1,    1); economic=false)
        parent(a) .= rank
        # * non-blocking comms don't work with corner points
        haloswap!(a)
        rank == 0 && @test parent(a) == [3 2 3; 
                                         1 0 1; 
                                         3 2 3]
        rank == 1 && @test parent(a) == [2 3 2; 
                                         0 1 0; 
                                         2 3 2]
        rank == 2 && @test parent(a) == [1 0 1; 
                                         3 2 3; 
                                         1 0 1]
        rank == 3 && @test parent(a) == [0 1 0; 
                                         2 3 2; 
                                         0 1 0]
    end

    @testset "swap 2                                 " begin
        # before the swap
        # 0 0 0 | 1 1 1 | 2 2 2 | 3 3 3
        # 0 0 0 | 1 1 1 | 2 2 2 | 3 3 3
        # 0 0 0 | 1 1 1 | 2 2 2 | 3 3 3

        # after the swap
        # 0 0 1 | 0 1 2 | 1 2 3 | 2 3 3
        # 0 0 1 | 0 1 2 | 1 2 3 | 2 3 3
        # 0 0 1 | 0 1 2 | 1 2 3 | 2 3 3
        rank = MPI.Comm_rank(MPI.COMM_WORLD)
        a = HaloArray{Float64}(MPI.COMM_WORLD,
                               (1,     4),
                               (false, false),
                               (1,     1),
                               (1,     1); economic=false)
        parent(a) .= rank
        haloswap!(a, rand([true, false]))
        MPI.Waitall(reqs(a))
        rank == 0 && @test parent(a) == [0 0 1; 
                                         0 0 1; 
                                         0 0 1]
        rank == 1 && @test parent(a) == [0 1 2; 
                                         0 1 2; 
                                         0 1 2]
        rank == 2 && @test parent(a) == [1 2 3; 
                                         1 2 3; 
                                         1 2 3]
        rank == 3 && @test parent(a) == [2 3 3; 
                                         2 3 3; 
                                         2 3 3]
    end
end

MPI.Finalize()
