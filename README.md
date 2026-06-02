<p align="center">
  <img src="assets/logo.svg" alt="HaloArrays.jl" width="720">
</p>

# HaloArrays.jl

`HaloArrays.jl` provides MPI-aware Julia arrays with halo, or ghost, cells for
block-structured stencil codes.

Each MPI rank owns a dense local block of data. `HaloArray` stores that local
interior together with halo cells around it, and `haloswap!` fills those halo
cells from neighbouring ranks in a Cartesian MPI topology. This is the common
data layout used in finite-difference, finite-volume, lattice, and other
nearest-neighbour methods on structured grids.

Use `HaloArrays.jl` when you want:

- local array code that looks like ordinary stencil code,
- halo cells addressable with natural indices such as `a[0, j]` or
  `a[nx + 1, j]`,
- blocking or non-blocking MPI halo exchange over Cartesian process grids,
- explicit control over physical boundary conditions on non-periodic edges.

It is not a global distributed-array abstraction: there is no global indexing,
load balancing, or mesh management. The package focuses on the rank-local array
and the halo exchange.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Davide-Lasagna-s-Lab/HaloArrays.jl")
```

## A First Stencil

Run this kind of script with MPI, for example:

```bash
mpiexec -n 4 julia --project=. stencil.jl
```

```julia
using MPI
using HaloArrays

MPI.Init()
try
    a = HaloArray{Float64}(
        MPI.COMM_WORLD,
        (2, 2),          # Cartesian process grid
        (true, true),    # periodic in x and y
        (64, 64),        # local interior size on each rank
        (1, 1),          # one halo cell on each side
    )

    rank = MPI.Comm_rank(comm(a))
    a .= rank            # broadcast touches the local interior only

    haloswap!(a)         # fill halo cells from neighbouring ranks

    lap = similar(a)
    nx, ny = size(a)
    for j in 1:ny, i in 1:nx
        lap[i, j] = a[i - 1, j] + a[i + 1, j] +
                    a[i, j - 1] + a[i, j + 1] -
                    4 * a[i, j]
    end
finally
    MPI.Finalize()
end
```

The loop only writes interior points, but it can read halo values at `i == 0`,
`i == nx + 1`, `j == 0`, and `j == ny + 1`.

## Core Model

A `HaloArray` separates the logical interior from the padded local storage:

```julia
size(a)       # local interior size
axes(a)       # interior axes, used by iteration and broadcast
nhalo(a)      # halo width in each dimension
parent(a)     # dense local storage, including halo cells
```

For a one-dimensional array with `size(a) == (5,)` and `nhalo(a) == (2,)`:

```text
logical index:  -1   0   1   2   3   4   5   6   7
region:          H   H   I   I   I   I   I   H   H
```

Interior Julia operations stay on the interior:

```julia
a .= 0
for I in CartesianIndices(a)
    # visits interior points only
end
```

Scalar indexing can also reach halo cells:

```julia
a[1, j]          # first interior point in dimension 1
a[0, j]          # first left halo point
a[-1, j]         # second left halo point, if nhalo(a)[1] >= 2
a[size(a, 1) + 1, j]  # first right halo point
```

If a halo index would fall outside the local padded storage, normal Julia bounds
checking throws a `BoundsError`.

## Construction

The most common constructor creates a Cartesian communicator from an existing
communicator:

```julia
a = HaloArray{Float64}(
    MPI.COMM_WORLD,
    (2, 2),          # number of ranks in each Cartesian direction
    (true, false),   # periodicity in each direction
    (128, 64),       # local interior size
    (1, 1),          # halo width
)
```

If you already have a Cartesian communicator, pass it directly:

```julia
cart = MPI.Cart_create(MPI.COMM_WORLD, (2, 2);
                       periodic=(true, false), reorder=false)

a = HaloArray{Float64}(cart, (128, 64), (1, 1))
```

Constructor options:

| Option | Default | Meaning |
| --- | --- | --- |
| `economic` | `true` | Exchange only face halos. This is enough for axis-aligned nearest-neighbour stencils such as a 5-point or 7-point stencil. |

The communicator is shared with the array. `HaloArray` does not duplicate it for
`similar` or `copy`, and it does not free it.

## Halo Exchange

For most codes, the blocking call is the right starting point:

```julia
haloswap!(a)
```

This exchanges all configured halo regions with neighbouring ranks. In
non-periodic directions, boundary ranks have no neighbour on one side; those
halo cells are left for your physical boundary condition code to set.

### Overlapping Work

With the default `economic=true`, a non-blocking exchange over all dimensions is
available:

```julia
requests = haloswap!(a, false)

# Do work that does not read the in-flight halos.

foreach(MPI.Waitall, requests)
```

`requests` is a tuple containing one MPI request vector per Cartesian dimension.
Complete those requests before starting another non-blocking exchange on the
same array.

### Stencils That Need Corners

The default exchange updates faces, not corners. That is deliberate: many common
stencils never read diagonal halo values.

Use `economic=false` when your stencil needs edge or corner halo values, for
example a 9-point stencil in 2D or a 27-point stencil in 3D:

```julia
a = HaloArray{Float64}(
    MPI.COMM_WORLD,
    (2, 2),
    (true, true),
    (64, 64),
    (1, 1);
    economic=false,
)

haloswap!(a)    # blocking exchange, including corner data
```

For non-blocking corner exchanges, stage the exchange dimension by dimension and
wait between stages:

```julia
r = haloswap!(a, 1; block=false)
# Work that does not need dimension-1 halos.
MPI.Waitall(r)

r = haloswap!(a, 2; block=false)
# Work that does not need dimension-2 halos.
MPI.Waitall(r)
```

The staged version is the non-blocking equivalent of the blocking
`haloswap!(a)`: each stage can carry halo data produced by previous stages, so
corner values are propagated correctly.

## Boundary Conditions

Periodic directions are handled by the Cartesian communicator. For non-periodic
directions, a rank at the physical boundary does not receive data from outside
the domain. Set those halo cells yourself before the stencil reads them:

```julia
nx, ny = size(a)

# Example: left physical boundary on this rank.
for j in 1:ny
    a[0, j] = 0.0
end
```

You can mix physical boundary conditions and MPI halo exchange in the usual way:
exchange processor boundaries with `haloswap!`, then fill the physical boundary
halos for the ranks that lie on the outside of the global domain.

## Allocation And Copying

`similar(a)` creates a new `HaloArray` with the same communicator, local size,
halo width, and exchange options:

```julia
rhs = similar(a)
```

`copy(a)` copies the full local padded storage, including halo cells:

```julia
b = copy(a)
parent(b) == parent(a)
```

Broadcast assignment, by contrast, follows the normal array interface and
updates the interior only:

```julia
b .= a
```

Use `parent(a)` when you deliberately want to inspect or initialize the complete
local storage, including halos.

## Public API

| Function | Purpose |
| --- | --- |
| `HaloArray{T}(comm, nprocesses, isperiodic, localsize, nhalo; kwargs...)` | Build a halo array and create a Cartesian communicator. |
| `HaloArray{T}(cart_comm, localsize, nhalo; kwargs...)` | Build a halo array on an existing Cartesian communicator. |
| `haloswap!(a)` | Blocking halo exchange over all configured dimensions. |
| `haloswap!(a, false)` | Non-blocking halo exchange over all dimensions. Supported for the default face-only exchange. |
| `haloswap!(a, dim; block=false)` | Start a non-blocking exchange in one Cartesian dimension. Useful for staged exchanges. |
| `reqs(a)` / `reqs(a, dim)` | Access reusable MPI request storage for all dimensions or for one dimension. |
| `nhalo(a)` | Halo width tuple. |
| `comm(a)` | Cartesian communicator used by the array. |
| `origin(a)` | Parent-storage index of the first interior point. |
| `parent(a)` | Dense local storage, including halo cells. |

## Limitations

- Halo widths are symmetric: each dimension has the same width on the left and
  right side.
- Halo exchange currently targets Cartesian MPI topologies with up to three
  exchanging dimensions.
- Non-blocking request storage is reused. Always wait on returned requests
  before reusing the same array and dimension.
- Communicator lifetime remains the caller's responsibility.

## Development

Run the test suite with:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
