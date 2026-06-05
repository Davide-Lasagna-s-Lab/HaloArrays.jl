# Halo Exchange

Halo exchange is performed with [`haloswap!`](@ref). It communicates with
neighbouring ranks in the Cartesian communicator returned by [`comm`](@ref).

## Blocking Exchange

The blocking call is the simplest and most robust option:

```julia
haloswap!(a)
```

It returns after all configured halo regions have been exchanged.

## Non-Blocking Exchange

For the default face-only exchange, use:

```julia
requests = haloswap!(a, false)

# Do work that does not read in-flight halo values.

foreach(MPI.Waitall, requests)
```

`reqs(a)` returns one `MPI.MultiRequest` per Cartesian dimension. `reqs(a, dim)`
returns the reusable request storage for one staged dimension. Wait on the
returned requests before starting another non-blocking exchange that reuses the
same storage.

## The Economic Keyword

The `economic` keyword controls the size of each exchanged strip in directions
transverse to the direction being exchanged.

With `economic=true`, the default, transverse directions use `CENTER`. The
message contains only the face data needed by axis-aligned nearest-neighbour
stencils:

```text
economic=true, union of 2D send footprints

    +---+---+---+
    | . | M | . |
    +---+---+---+
    | M | I | M |
    +---+---+---+
    | . | M | . |
    +---+---+---+

I: local interior
M: sent face footprint
.: not sent
```

Use this for stencils such as:

```text
    .   x   .
    x   u   x
    .   x   .
```

With `economic=false`, transverse directions use `ALL`. Messages include halo
values in the transverse direction, so staged exchanges can propagate edge and
corner values:

```text
economic=false, union of 2D send footprints

    +---+---+---+
    | M | M | M |
    +---+---+---+
    | M | I | M |
    +---+---+---+
    | M | M | M |
    +---+---+---+
```

Use this for stencils such as:

```text
    x   x   x
    x   u   x
    x   x   x
```

Blocking all-dimension exchange supports `economic=false`:

```julia
a = HaloArray{Float64}(
    MPI.COMM_WORLD, (2, 2), (true, true), (64, 64), (1, 1);
    economic=false,
)

haloswap!(a)
```

For non-blocking corner-capable exchange, stage the dimensions and wait between
stages:

```julia
r = haloswap!(a, 1, false)
MPI.Waitall(r)

r = haloswap!(a, 2, false)
MPI.Waitall(r)
```

The later stage sends halo values produced by the earlier stage, which is how
corner data moves through a sequence of one-dimensional exchanges.

