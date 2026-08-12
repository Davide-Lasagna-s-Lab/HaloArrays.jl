# CUDA Extension

HaloArrays.jl provides a package extension (`HaloArraysCUDAExt`) that adds
GPU support via [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl). The extension
is loaded automatically when both `HaloArrays` and `CUDA` are present in your
environment.

## Moving data to the GPU

A `HaloArray` can be moved to the GPU with `CUDA.cu`:

```julia
using HaloArrays, CUDA

A = HaloArray(rand(128, 128), nhalo=(1, 1), comm)
A_gpu = cu(A)
```

`cu(A)` returns a `HaloArray` whose underlying storage and communication
buffers are backed by `CuArray`s. All other fields (the MPI communicator and
request objects) are left untouched.

### CUDA-aware MPI

Halo exchange on a GPU-backed `HaloArray` uses CUDA-aware MPI, i.e. `MPI.Isend`/
`MPI.Irecv!` are called directly on device buffers without staging through the
host. This requires an MPI implementation and system configuration that
support CUDA-aware transport. Consult your MPI implementation's documentation
and verify your environment before relying on this path — if CUDA-aware
support is not correctly configured, communication calls may fail or silently
corrupt data.

## Using HaloArrays in kernels

A `HaloArray` can be passed directly as an argument to a `@cuda` kernel. It is
automatically converted, via `Adapt.jl`, into a lightweight device-side
representation exposing only the underlying data and halo metadata:

```julia
function deriv_kernel!(out, A, h)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if 1 <= i <= size(A, 1)
        @inbounds out[i] = (A[i + 1] - A[i - 1]) / (2h)
    end
    return
end

@cuda threads=128 deriv_kernel!(out, A_gpu, h)
```

Indexing into the device-side array is offset by the number of halo points,
so index `1` refers to the first local (interior) element, and indices `<= 0`
resolve into the neighboring ghost region. No manual offset bookkeeping is
required inside kernel code.
