# HaloArrays.jl

`HaloArrays.jl` provides MPI-aware Julia arrays with halo, or ghost, cells for
structured-grid stencil codes.

Each MPI rank owns a dense local block. A [`HaloArray`](@ref) stores that
interior block together with symmetric halo cells around it, and
[`haloswap!`](@ref) exchanges boundary data with neighbouring ranks in a
Cartesian MPI topology.

## Main Features

- Rank-local dense arrays with halo storage on both sides of each dimension.
- Ordinary Julia iteration and broadcast over the local interior.
- Scalar indexing that can address halo cells with indices such as `a[0, j]`
  and `a[nx + 1, j]`.
- Blocking and non-blocking halo exchange through MPI.jl.
- Face-only exchange for common nearest-neighbour stencils.
- Wider staged exchange for stencils that need edge or corner halo values.

## What HaloArrays.jl Does Not Do

`HaloArrays.jl` is intentionally small. It does not provide global indexing,
mesh generation, load balancing, or a global distributed-array abstraction. It
focuses on the local storage and communication pattern used by many
finite-difference, finite-volume, lattice, and structured-grid methods.

## A Minimal Example

```julia
using MPI
using HaloArrays

MPI.Init()
try
    a = HaloArray{Float64}(
        MPI.COMM_WORLD,
        (2, 2),          # Cartesian process grid
        (true, true),    # periodic in both dimensions
        (64, 64),        # local interior size per rank
        (1, 1),          # one halo cell per side
    )

    a .= MPI.Comm_rank(comm(a))
    haloswap!(a)

    nx, ny = size(a)
    out = similar(a)
    for j in 1:ny, i in 1:nx
        out[i, j] = a[i - 1, j] + a[i + 1, j] +
                    a[i, j - 1] + a[i, j + 1] -
                    4a[i, j]
    end
finally
    MPI.Finalize()
end
```

Run scripts like this with MPI:

```bash
mpiexec -n 4 julia --project=. stencil.jl
```

## Manual

Start with [Getting started](manual/getting_started.md), then read
[Halo indexing](manual/indexing.md) and [Halo exchange](manual/exchange.md) for
the communication model and the `economic` keyword.
