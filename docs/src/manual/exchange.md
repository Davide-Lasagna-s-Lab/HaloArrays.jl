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

The diagrams below use one notation throughout:

```text
H    halo storage, not communicated in this diagram
I    interior storage, not communicated in this diagram
[H]  halo storage that is communicated
[I]  interior storage that is communicated
```

The letter is always the storage class. Brackets mark the cells touched by the
message or exchange being shown.

## Choosing `economic`

`economic` is a Boolean communication policy. It should be chosen by the user
or by the higher-level distributed-grid object and then passed through to
`HaloArray` construction. It should not be `nothing`, and it should not be
inferred mechanically from stencil width, derivative order, or the number of
decomposed dimensions.

Choose `economic=true` when the algorithm only reads axis-aligned face halos
after a swap. Choose `economic=false` when edge or corner halo values must also
be valid. For example, a wide one-dimensional finite-difference stencil along a
single decomposed axis may still only need face halos, while a compact 2D or 3D
stencil that reads diagonal neighbours needs corners.

With `economic=true`, the default, transverse directions use `CENTER`. For a
single `LEFT` send in dimension 1, only the interior span of dimension 2 is
sent:

```text
economic=true: region = (LEFT, CENTER)

                         dimension 2
               halo        interior        halo
            +---------+----------------+---------+
dim 1 halo  |    H    |       H        |    H    |
            +---------+----------------+---------+
LEFT strip  |    H    |      [I]       |    H    |
            +---------+----------------+---------+
interior    |    H    |       I        |    H    |
            +---------+----------------+---------+
dim 1 halo  |    H    |       H        |    H    |
            +---------+----------------+---------+
```

That is enough for axis-aligned nearest-neighbour stencils:

```text
    .   x   .
    x   u   x
    .   x   .
```

With `economic=false`, transverse directions use `ALL`. The same `LEFT` send
includes transverse halo storage:

```text
economic=false: region = (LEFT, ALL)

                         dimension 2
               halo        interior        halo
            +---------+----------------+---------+
dim 1 halo  |    H    |       H        |    H    |
            +---------+----------------+---------+
LEFT strip  |   [H]   |      [I]       |   [H]   |
            +---------+----------------+---------+
interior    |    H    |       I        |    H    |
            +---------+----------------+---------+
dim 1 halo  |    H    |       H        |    H    |
            +---------+----------------+---------+
```

The wider message is useful for stencils that eventually need edge or corner
halo values:

```text
    x   x   x
    x   u   x
    x   x   x
```

After all dimensions have been exchanged, the received halo storage differs as
follows:

```text
economic=true

    +---+---+---+
    | H |[H]| H |
    +---+---+---+
    |[H]| I |[H]|
    +---+---+---+
    | H |[H]| H |
    +---+---+---+

economic=false, staged

    +---+---+---+
    |[H]|[H]|[H]|
    +---+---+---+
    |[H]| I |[H]|
    +---+---+---+
    |[H]|[H]|[H]|
    +---+---+---+
```

With `economic=true`, only face halos are filled. With `economic=false`, staged
exchanges can fill edge and corner halos by forwarding halo values from one
dimension through the next.

## Corner Halos and Staged Exchange

`HaloArrays.jl` does not send corner or edge regions as explicit diagonal
messages. For example, in 2D it does not construct regions such as
`(LEFT, LEFT)` or communicate directly with diagonal MPI ranks. Corner values
are instead propagated by a sequence of one-dimensional exchanges:

```text
2D economic=false, staged order: dimension 1, then dimension 2

after dim 1 waits             after dim 2 waits

    +---+---+---+                 +---+---+---+
    |[H]|[H]|[H]|                 |[H]|[H]|[H]|
    +---+---+---+                 +---+---+---+
    | H | I | H |                 |[H]| I |[H]|
    +---+---+---+                 +---+---+---+
    |[H]|[H]|[H]|                 |[H]|[H]|[H]|
    +---+---+---+                 +---+---+---+
```

The first stage updates halos in dimension 1. The second stage sends strips in
dimension 2 that include the dimension-1 halo values, so diagonal data reaches
the corner without a direct diagonal send.

This staged approach is common in structured-grid halo exchange. It avoids
communicating with diagonal neighbours and keeps the implementation close to
axis-aligned neighbour exchange, but it introduces an ordering dependency:
later dimensions cannot be started until the previous stage has completed if
corner values are required. Explicit corner messages can avoid this dependency
and may reduce synchronization, but they add more neighbour pairs, buffer
descriptions, tags, and message matching work. The right choice is
machine- and stencil-dependent. MiniGhost, for example, was designed to explore
these boundary-exchange strategy tradeoffs, including aggregation of boundary
data to reduce message count.

The current `economic=false` implementation is conservative rather than fully
optimal. It uses `ALL` in every transverse dimension for every staged direction.
That guarantees that a later stage has the halo values it may need, but earlier
stages can send transverse halo values that have not yet been refreshed by any
previous stage. A leaner implementation could make the staged regions depend on
the exchange order: dimensions that have already been exchanged would use
`ALL`, while dimensions still waiting for a later stage would use `CENTER`.
This would preserve corner propagation while reducing the amount of halo data
sent by early stages.

Blocking all-dimension exchange supports `economic=false` and performs the
necessary staged exchange internally:

```julia
a = HaloArray{Float64}(
    MPI.COMM_WORLD, (2, 2), (true, true), (64, 64), (1, 1);
    economic=false,
)

haloswap!(a)
```

All-at-once non-blocking exchange is rejected for `economic=false` when more
than one dimension is exchanged. Users should stage the exchange along
dimensions and wait between stages:

```julia
r = haloswap!(a, 1, false)
MPI.Waitall(r)

r = haloswap!(a, 2, false)
MPI.Waitall(r)
```

The wait between stages is part of the algorithm, not just an implementation
detail. The later stage sends halo values produced by the earlier stage, which
is how corner data moves through a sequence of one-dimensional exchanges.

## Further Reading

- Kjolstad and Snir, [Ghost Cell Pattern](https://doi.org/10.1145/1953611.1953615),
  ParaPLoP 2010.
- Barrett, Vaughan, and Heroux, [MiniGhost: A Miniapp for Exploring Boundary
  Exchange Strategies Using Stencil Computations in Scientific Parallel
  Computing](https://doi.org/10.2172/1039405), Sandia National Laboratories,
  2012.
- Huenich and Knuepfer, [A Halo abstraction for distributed n-dimensional
  structured grids within the C++ PGAS library
  DASH](https://doi.org/10.7717/peerj-cs.1203), PeerJ Computer Science, 2023.
- Brodtkorb and Saetra, [Simulating the Euler equations on multiple GPUs using
  Python](https://doi.org/10.3389/fphy.2022.985440), Frontiers in Physics,
  2022. Section 2.4 describes transferring corners by staged north-south then
  east-west exchange, which avoids diagonal transfers at the cost of an extra
  synchronization point.
