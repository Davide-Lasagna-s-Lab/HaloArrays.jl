module CUDAExt

using Adapt,
      CUDA

import CUDA: i32

using HaloArrays

struct HaloDeviceArray{T, N, NHALO, SIZE, A} <: HaloArrays.AbstractHaloArray{T, N, NHALO, SIZE}
    data::A

    HaloDeviceArray{NHALO, SIZE}(data::A) where {T, N, NHALO, SIZE, A<:AbstractArray{T, N}} =
        new{T, N, NHALO, SIZE, A}(data)
end

Base.similar(a::HaloDeviceArray{T}) where {T} = similar(a, T)
Base.similar(a::HaloDeviceArray, ::Type{T}) where {T} =
    HaloDeviceArray{nhalo(a), size(a)}(similar(parent(a), T))


function Adapt.adapt_structure(to, a::HaloArray{T, N, NHALO, SIZE}) where {T, N, NHALO, SIZE}
    data = Adapt.adapt_structure(to, a.data)
    return HaloDeviceArray{Int32.(NHALO), Int32.(SIZE)}(data)
end

function CUDA.cu(a::HaloArray{T, N, NHALO, SIZE}) where {T, N, NHALO, SIZE}
    data_d = CUDA.cu(a.data)
    return HaloArray{NHALO, SIZE}(data_d, a.haloregions, a.comm, a.reqs, a.economic)
end

end
