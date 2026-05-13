import HaloArrays: HaloArray,
                   comm

using BenchmarkTools
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

a = HaloArray{Float64}(MPI.COMM_WORLD,
                      (    1,     1,     1),
                      (false, false, false),
                      (  256,   256,   256),
                      (    0,     0,     0))
pa = parent(a)

if MPI.Comm_rank(comm(a)) == 0
    t_haloarray = @belapsed $mysum($a)
    t_parent = @belapsed $mysum($pa)
    println("HaloArray: ", t_haloarray)
    println("parent:    ", t_parent)
    println("ratio:     ", t_haloarray / t_parent)
end

MPI.Finalize()
