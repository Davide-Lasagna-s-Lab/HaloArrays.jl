# Internals

This page describes the implementation model used by `HaloArrays.jl`. It is
intended for contributors and for users who need to reason about communication
costs, buffer layout, or staged halo exchanges.

## Local Storage Model

A [`HaloArray`](@ref) is a rank-local dense parent array plus metadata that
describes how the parent is split into interior and halo regions.

```text
1D parent storage with nhalo(a) == (2,) and size(a) == (5,)

logical index:   -1   0 | 1   2   3   4   5 | 6   7
parent index:     1   2 | 3   4   5   6   7 | 8   9
storage role:    HL  HL | I   I   I   I   I | HR  HR
```

The public array size is the local interior size. Scalar indexing shifts
logical indices by `nhalo(a)` before indexing the parent storage, so halo cells
can be addressed with zero, negative, or `size(a, d) + k` indices when those
indices still map inside the parent array.

Broadcast, `axes(a)`, and `CartesianIndices(a)` use the interior domain only.
Use [`parent`](@ref) when code deliberately needs the complete local storage,
including halos.

## Region Tuples

Halo exchange regions are represented by `NTuple{N, Region}` values. Each tuple
has exactly one exchanged dimension, marked by `LEFT` or `RIGHT`. The remaining
dimensions describe the transverse extent of the message:

```text
2D face-only exchange

    (LEFT,  CENTER)     left face
    (RIGHT, CENTER)     right face
    (CENTER, LEFT)      bottom face
    (CENTER, RIGHT)     top face

2D wider staged exchange

    (LEFT,  ALL)        left strip, including transverse halos
    (RIGHT, ALL)        right strip, including transverse halos
    (ALL,   LEFT)       bottom strip, including transverse halos
    (ALL,   RIGHT)      top strip, including transverse halos
```

`CENTER` means the interior span of one direction. `ALL` means the complete
parent span in that direction, including both halo sides.

The helper [`HaloArrays.regions_to_swap`](@ref) chooses these tuples from the
Cartesian process grid, periodicity, and `economic` option. Dimensions are
included when there is more than one process in that direction or when the
direction is periodic.

`economic` is expected to be a concrete Boolean policy selected by the caller.
Higher-level distributed-grid wrappers should pass and store that value
explicitly instead of deriving it from local details such as stencil width,
derivative order, or the number of decomposed dimensions.

## Send And Receive Slices

The same `Region` value selects different parent slices depending on whether it
is used for sending or receiving. In 1D with one halo cell:

```text
parent index:   1 | 2   3   4   5   6 | 7
storage role:  HL | I   I   I   I   I | HR

region = LEFT
    SEND -> 2        interior strip next to the left halo
    RECV -> 1        left halo storage

region = RIGHT
    SEND -> 6        interior strip next to the right halo
    RECV -> 7        right halo storage
```

[`HaloArrays.subarray_slices`](@ref) computes these slices. The constructor
wraps each slice in an `MPI.Buffer` and stores it in
`a.buffers[(region, intent)]`, where `intent` is `SEND` or `RECV`.

The buffer dictionary is currently typed as `MPI.Buffer`, the abstract MPI.jl
buffer supertype. This is intentional for now: `MPI.Buffer(view(...))` may
return different concrete buffer types depending on whether the view is
contiguous. A more specialized representation could reduce type instability,
but it would need to be justified by profiling.

## Blocking Exchange

Blocking exchange uses `MPI.Sendrecv!` for each configured region. For a region
on rank `p`, HaloArrays sends from that region to the destination rank returned
by `MPI.Cart_shift`, and receives into the opposite region:

```text
send region:    region
receive region: opposite(region)

example:
    send    (LEFT, CENTER)
    receive (RIGHT, CENTER)
```

Using `Sendrecv!` keeps the blocking path deadlock-free while preserving a
simple region-by-region implementation.

## Non-Blocking Exchange

Non-blocking exchange uses one reusable `MPI.MultiRequest` per Cartesian
dimension. `reqs(a, dim)` returns the request storage for a staged exchange in
one dimension; `reqs(a)` returns the tuple for all dimensions.

Each region in a dimension consumes two request slots: one receive and one send.
The slots are reused by later non-blocking calls, so callers must complete the
requests before starting another exchange that uses the same dimension:

```julia
r = haloswap!(a, 1, false)
MPI.Waitall(r)
```

All-at-once non-blocking exchange is allowed for the default face-only exchange.
For `economic=false` with multiple exchanged dimensions, users should stage the
exchange manually and wait between dimensions so corner values can be
propagated in a defined order.

## Corners And Current Tradeoffs

HaloArrays does not send explicit corner or edge regions such as
`(LEFT, LEFT)`. Corner values are propagated by staged one-dimensional
exchanges. A later stage can send halo values received by an earlier stage:

```text
2D economic=false, staged dim 1 then dim 2

dim 1 stage fills horizontal halo strips
dim 2 stage sends vertical strips that include those dim-1 halo values
```

This avoids diagonal neighbours and keeps the communication model close to
axis-aligned Cartesian exchange. The tradeoff is that corner propagation has an
ordering dependency: the later stage needs the earlier stage to have completed.

The current `economic=false` implementation is conservative. It uses `ALL` for
every transverse direction in every staged direction. That is correct, but not
minimal: early stages may send transverse halo values that have not yet been
updated in the current halo exchange. A future stage-aware implementation could
use `ALL` only in dimensions that have already completed, and `CENTER` in
dimensions still waiting for a later stage.

## Communicator Ownership

The communicator stored in a `HaloArray` is shared with the caller. The array
does not own it and does not free it. Arrays produced by [`similar`](@ref) and
[`copy`](@ref) reuse the same communicator.

## Internal API

These helpers are documented for contributors. They are not exported and should
not be treated as a stable public API.

```@docs
HaloArrays.regions_to_swap
HaloArrays.opposite
HaloArrays.subarray_slices
HaloArrays.source_dest_ranks
HaloArrays.exchange_dim
HaloArrays.checkdim
HaloArrays.haloregions
HaloArrays._haloswap!
HaloArrays._Ihaloswap!
```
