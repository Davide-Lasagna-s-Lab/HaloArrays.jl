function _test_indexing_kernel!(out, A::CUDAExt.HaloDeviceArray{T, N, NHALO, SIZE}) where {T, N, NHALO, SIZE}
    idx = threadIdx().x

    rem = idx - 1i32
    i   = rem % (SIZE[1] + 2*NHALO[1]) + 1i32
    rem = rem ÷ (SIZE[1] + 2*NHALO[1])
    j   = rem % (SIZE[2] + 2*NHALO[2]) + 1i32
    rem = rem ÷ (SIZE[2] + 2*NHALO[2])
    k   = rem % (SIZE[3] + 2*NHALO[3]) + 1i32
    rem = rem ÷ (SIZE[3] + 2*NHALO[3])
    I = (i, j, k)

    @inbounds out[I...] = A[(I .- NHALO)...]

    return nothing
end

@testset "CUDAExt - HaloArray adaptation         " begin
    a = HaloArray{Float64}(comm,
                           (   1,    1,    1),
                           (true, true, true),
                           (   1,    1,    1),
                           (   1,    1,    1))
    parent(a) .= randn(Float64, size(parent(a)))

    ad       = @test_nowarn CUDA.cu(a)
    ad_adapt = @test_nowarn Adapt.adapt(CUDA.KernelAdaptor(), ad)

    @test ad       isa         HaloArray{      Float32, 3, (1,    1,    1),    (1,    1,    1),    <:CuArray{Float32}}
    for val in values(ad.buffers); @test val.data isa CuArray; end
    @test ad_adapt isa CUDAExt.HaloDeviceArray{Float32, 3, (1i32, 1i32, 1i32), (1i32, 1i32, 1i32), <:CuDeviceArray{Float32}}
    @test isbits(ad_adapt)

    @test similar(ad) isa typeof(ad)
    @test    copy(ad) isa typeof(ad)

    CUDA.allowscalar() do
        for k = 0:2, j = 0:2, i = 0:2
            @test ad[i, j, k] == Float32(a[i, j, k])
        end
        @test_throws BoundsError ad[-1, 0, 1]
    end

    out = CUDA.zeros(Float32, size(parent(a)))
    @cuda threads=length(parent(a)) _test_indexing_kernel!(out, ad)
    @test Array(out) == Array(parent(ad))
end

MPI.Finalize()
