import HaloArrays: HaloArray,
                   comm

using BenchmarkTools
using Test
using MPI

function mysum(a)
    S = zero(eltype(a))
    for k = axes(a, 3)
        for j = axes(a, 2)
            @simd for i = axes(a, 1)
                @inbounds S += a[i, j, k]
            end
        end
    end
    return S
end

MPI.Init()

@testset "sum                                    " begin
    # use smaller sizes when #33927 is solved
    a = HaloArray{Float64}(MPI.COMM_WORLD,
                          (    1,     1,     1),
                          (false, false, false),
                          (  256,   256,   256),
                          (    0,     0,     0))
    pa = parent(a)
    if MPI.Comm_rank(comm(a)) == 0
        t1 = @belapsed $mysum($a)
        t2 = @belapsed $mysum($pa)
        @test (t1-t2)/t1 < 0.05
    end
end

MPI.Finalize()
