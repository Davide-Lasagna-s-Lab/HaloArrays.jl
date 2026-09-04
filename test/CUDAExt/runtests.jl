using Test,
      MPI,
      CUDA,
      Adapt

using HaloArrays

import CUDA: i32

# setup MPI environment
MPI.Init()
comm = MPI.COMM_WORLD
np   = MPI.Comm_size(comm)
rank = MPI.Comm_rank(comm)

CUDAExt = Base.get_extension(HaloArrays, :CUDAExt)

function cuda_available()
    try
        CUDA.functional()
    catch
        false
    end
end

if !cuda_available()
    @warn "Skipping GPU tests - CUDA not functional"
else
    include("test_cudaext.jl")
end
