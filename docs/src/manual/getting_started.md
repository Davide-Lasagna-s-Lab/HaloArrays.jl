# Getting Started

## Installation

During development, install the package from the repository:

```julia
using Pkg
Pkg.add(url="https://github.com/Davide-Lasagna-s-Lab/HaloArrays.jl")
```

`HaloArrays.jl` depends on MPI.jl, so your scripts should initialize MPI before
constructing arrays and finalize it before exiting.

## Constructing A HaloArray

The most common constructor creates a Cartesian communicator from an existing
communicator:

```julia
a = HaloArray{Float64}(
    MPI.COMM_WORLD,
    (2, 2),          # number of ranks in each Cartesian direction
    (true, false),   # periodicity in each direction
    (128, 64),       # local interior size on each rank
    (1, 1),          # halo width in each direction
)
```

If you already have a Cartesian communicator, pass it directly:

```julia
cart = MPI.Cart_create(MPI.COMM_WORLD, (2, 2);
                       periodic=(true, false), reorder=false)

a = HaloArray{Float64}(cart, (128, 64), (1, 1))
```

The communicator is shared with the array. It is not duplicated by
[`similar`](@ref) or [`copy`](@ref), and `HaloArray` does not free it.

## Local Computation

`size(a)` and `axes(a)` describe the local interior only. That means broadcast
and normal Julia iteration work naturally on the values owned by the rank:

```julia
a .= 0.0

for I in CartesianIndices(a)
    a[I] = 1.0
end
```

Stencil loops can use scalar indexing to read halo values after a halo exchange:

```julia
haloswap!(a)

nx, ny = size(a)
out = similar(a)
for j in 1:ny, i in 1:nx
    out[i, j] = a[i - 1, j] + a[i + 1, j] +
                a[i, j - 1] + a[i, j + 1] -
                4a[i, j]
end
```

## Boundary Conditions

Periodic boundaries are handled by the Cartesian communicator. For non-periodic
directions, MPI has no neighbour outside the domain; set those halo cells in
your boundary-condition code:

```julia
nx, ny = size(a)

# Example: a left physical boundary on this rank.
for j in 1:ny
    a[0, j] = 0.0
end
```

