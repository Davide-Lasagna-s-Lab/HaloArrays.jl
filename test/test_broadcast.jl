import HaloArrays: HaloArray

using Test
using MPI

MPI.Init()

@testset "broadcasting allocation                " begin
    a = HaloArray{Float64}(MPI.COMM_WORLD,
                           (   1,    1,    1),
                           (true, true, true),
                           (   1,    1,    1),
                           (   1,    1,    1))
    b = copy(a)
    c = copy(a)
    d = copy(a)
    e = copy(a)
    f = copy(a)

    # test broadcasting does not allocate
    foo(a, b, c, d, e, f) = (@allocated a .= 2.0.*b .+ 3.0.*c .+ 4.0.*d .+ 5.0.*e)

    @test foo(a, b, c, d, e, f) == 0
end

MPI.Finalize()
