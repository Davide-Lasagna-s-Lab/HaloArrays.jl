"""
    HaloArrays

MPI-distributed dense arrays with halo cells.

`HaloArray` stores a local interior array together with symmetric halo cells on
both sides of each Cartesian direction. The public `size` and default Julia
iteration domain describe only the interior region. Scalar indexing is
halo-aware: indices such as `0` or negative values are accepted when, after
shifting by `nhalo(a)`, they still refer to storage inside `parent(a)`.

In one dimension, a local array with five interior points and one halo cell on
each side is laid out as:

```text
logical index:    0 | 1   2   3   4   5 | 6
storage role:    HL | I   I   I   I   I | HR
parent index:     1 | 2   3   4   5   6 | 7

HL, HR: left/right halo storage
I:      local interior storage
```

This lets stencil code write expressions such as `a[i-1, j]` at interior
boundaries while keeping broadcast and ordinary array iteration focused on the
interior.
"""
module HaloArrays

#=
Resources:

[1] https://pdfs.semanticscholar.org/9f38/75a9a639d2371708e64b15b9556fe1eb4039.pdf
[2] http://www.hector.ac.uk/cse/distributedcse/reports/cp2k/cp2k/node9.html

=#

using MPI

export HaloArray,
       nhalo,
       comm,
       reqs,
       origin,
       haloswap!

# ------------------------------------------------------------------------------
# Region descriptors
# ------------------------------------------------------------------------------

"""
    Intent

Internal enum describing whether a halo region is used as a `SEND` region or a
`RECV` region when constructing MPI buffers.
"""
@enum Intent SEND RECV

"""
    Region

Internal enum describing the position of a subarray in each Cartesian direction.

`LEFT` and `RIGHT` denote halo-adjacent faces, `CENTER` denotes the interior
range in that direction, and `ALL` denotes the full parent range including halo
cells. A tuple of `Region` values uniquely identifies a halo face, edge, or
corner once the local size and halo width are known.

The enum is spelled `CENTER` in code. It denotes the centre/interior interval of
one Cartesian direction:

```text
one dimension of parent(a)

       left side        interior         right side
    +-------------+----------------+-------------+
    |    LEFT     |     CENTER     |    RIGHT    |
    +-------------+----------------+-------------+

    ALL covers the complete line, including both halo sides.
```

`Region` values are combined across dimensions. In 2D, for example,
`(LEFT, CENTER)` is the left face and `(LEFT, LEFT)` is the lower-left corner.
For halo exchange the code normally constructs face regions such as
`(LEFT, CENTER)` or wider regions such as `(LEFT, ALL)`.
"""
@enum Region LEFT RIGHT CENTER ALL

"""
    regions_to_swap(nprocesses::NTuple{N, Int},
                isperiodic::NTuple{N, Bool},
                economic::Bool) where {N}

Construct an array of `N`-tuples of the enum type `Region` specifying the
regions where halo exchange between processors or within the same processor
should happen. The exchange should happen along directions `dim` where
`nprocesses[dim] > 1` or `isperiodic[dim]` is true. At most three exchange
directions are currently supported.

The `economic` argument controls the transverse extent of each message. The
exchanged dimension is always `LEFT` or `RIGHT`; the other dimensions are either
`CENTER` (`economic=true`) or `ALL` (`economic=false`).

`economic` is a communication policy chosen by the caller. It should be a
concrete `Bool`, not `nothing`, and higher-level wrappers should pass and store
that choice explicitly. Do not infer it mechanically from stencil width,
derivative order, or the number of decomposed dimensions: choose
`economic=true` when the algorithm only needs face halos, and `economic=false`
when edge or corner halo values must be valid after the exchange.

For a 2D decomposition, consider one `LEFT` send in dimension 1. The
source region is the `LEFT` interior strip in dimension 1. Dimension 2 is
transverse to the exchange.

```text
Legend:
    H    halo storage, not in this message
    I    interior storage, not in this message
    [H]  halo storage included in this message
    [I]  interior storage included in this message

The letter is always the storage class. Brackets mark the message.

economic = true: region = (LEFT, CENTER)

                              dimension 2
                    halo        interior        halo
                 +---------+----------------+---------+
dim 1 halo       |    H    |       H        |    H    |
                 +---------+----------------+---------+
LEFT strip       |    H    |      [I]       |    H    |
                 +---------+----------------+---------+
interior         |    H    |       I        |    H    |
                 +---------+----------------+---------+
dim 1 halo       |    H    |       H        |    H    |
                 +---------+----------------+---------+

Only the interior span of dimension 2 is sent. This is economical because it
avoids sending transverse halo data. It fills face halos, but not corner halos.

economic = false: region = (LEFT, ALL)

                              dimension 2
                    halo        interior        halo
                 +---------+----------------+---------+
dim 1 halo       |    H    |       H        |    H    |
                 +---------+----------------+---------+
LEFT strip       |   [H]   |      [I]       |   [H]   |
                 +---------+----------------+---------+
interior         |    H    |       I        |    H    |
                 +---------+----------------+---------+
dim 1 halo       |    H    |       H        |    H    |
                 +---------+----------------+---------+

The full span of dimension 2 is sent, including halo values. After staged
exchanges this allows data to travel through edges and corners.

Destination halo storage filled after a complete 2D exchange

economic=true                  economic=false, staged

    +---+---+---+                  +---+---+---+
    | H |[H]| H |                  |[H]|[H]|[H]|
    +---+---+---+                  +---+---+---+
    |[H]| I |[H]|                  |[H]| I |[H]|
    +---+---+---+                  +---+---+---+
    | H |[H]| H |                  |[H]|[H]|[H]|
    +---+---+---+                  +---+---+---+

This last summary shows receive locations, not source regions. With
`economic=true`, only face halos are filled. With `economic=false`, staged
exchanges can fill edge and corner halos by forwarding halo values from one
dimension through the next.

HaloArrays does not create explicit corner or edge send regions such as
`(LEFT, LEFT)`. This follows the common structured-grid approach of propagating
corner values through ordered one-dimensional exchanges rather than sending
directly to diagonal neighbours. The tradeoff is fewer neighbour directions and
simpler buffer descriptions, at the cost of a wait between staged non-blocking
exchanges when corner values are needed.

The current `economic=false` regions are conservative: every staged direction
uses `ALL` in every transverse dimension. This is correct but not minimal,
because early stages can send transverse halo values that have not yet been
refreshed by a previous stage. A more optimal staged implementation would use
`ALL` only in dimensions that have already been exchanged, and `CENTER` in
dimensions still waiting for a later stage.
```
"""
function regions_to_swap(nprocesses::NTuple{N, Int},
                         isperiodic::NTuple{N, Bool}, economic::Bool) where {N}
    # vector of dimensions along which exchange must happen
    exdims = findall( (nprocesses .!= 1) .| isperiodic )

    # number of nprocesses over which exchange must happen
    ndims = length(exdims)

    # we do not expect 4D transfers
    ndims > 3 && throw(ArgumentError("invalid argument"))

    # init
    regions = NTuple{N, Region}[]

    # Efficient transfers of areas. With `economic=false` this deliberately
    # uses `ALL` in every transverse dimension, which is correct but not fully
    # minimal for staged corner propagation.
    #
    # TODO: make `economic=false` regions stage-aware. A leaner exchange would
    # use `ALL` only for dimensions that have already completed earlier stages
    # and `CENTER` for dimensions that will be exchanged later.
    AREA = economic == true ? CENTER : ALL

    if ndims ≥ 1
        push!(regions, ntuple(i -> i == exdims[1] ? LEFT  : AREA, N))
        push!(regions, ntuple(i -> i == exdims[1] ? RIGHT : AREA, N))
        if ndims ≥ 2
            push!(regions, ntuple(i -> i == exdims[1] ? AREA :
                                       i == exdims[2] ? LEFT   : AREA, N))
            push!(regions, ntuple(i -> i == exdims[1] ? AREA :
                                       i == exdims[2] ? RIGHT  : AREA, N))
            if ndims == 3
                push!(regions, ntuple(i -> i == exdims[1] ? AREA :
                                           i == exdims[2] ? AREA :
                                           i == exdims[3] ? LEFT  : AREA, N))
                push!(regions, ntuple(i -> i == exdims[1] ? AREA :
                                           i == exdims[2] ? AREA :
                                           i == exdims[3] ? RIGHT  : AREA, N))
            end
        end
    end

    # also add the opposite regions
    return regions
end


"""
    opposite(regions::NTuple{N, Region}) where {N}

Return the region tuple on the opposite side of the local array.

`LEFT` is changed to `RIGHT`, `RIGHT` is changed to `LEFT`, and `CENTER`/`ALL`
are left unchanged. This is used when a send buffer for one side of the local
domain is paired with the receive buffer on the opposite side.

```text
opposite((LEFT,  CENTER)) == (RIGHT, CENTER)
opposite((RIGHT, ALL))    == (LEFT,  ALL)
opposite((ALL,   LEFT))   == (ALL,   RIGHT)
```
"""
opposite(regions::NTuple{N, Region}) where {N} =
    ntuple(i -> regions[i] == LEFT  ? RIGHT :
                regions[i] == RIGHT ? LEFT  : regions[i], N)

"""
    subarray_slices(localsize::NTuple{N, Int},
                        nhalo::NTuple{N, Int},
                       region::NTuple{N, Region},
                       intent::Intent) where {N}

Calculate the slice corresponding to the subarray of the array of the total halo array
(including halo points), corresponding to the halo region specified by `region` and
`intent`, where the array has `nhalo` halo points on the left and right boundaries of
the `N` cartesian directions.

The returned tuple indexes the parent storage, not the logical interior indexing
used by `getindex(::HaloArray, ...)`.

The same `Region` selects different storage depending on the `Intent`.
With one halo cell and five interior cells:

```text
parent index:   1 | 2   3   4   5   6 | 7
storage role:  HL | I   I   I   I   I | HR

region = LEFT
    RECV -> 1        fill the left halo
    SEND -> 2        send the interior strip next to the left halo

region = RIGHT
    SEND -> 6        send the interior strip next to the right halo
    RECV -> 7        fill the right halo

region = CENTER -> 2:6
region = ALL    -> 1:7
```
"""
subarray_slices(localsize::NTuple{N, Int},
                    nhalo::NTuple{N, Int},
                   region::NTuple{N, Region},
                   intent::Intent) where {N} =
    ntuple(N) do i
        region[i] == LEFT   && return (intent == RECV ? (1:nhalo[i]) : (nhalo[i]+1:2*nhalo[i]))
        region[i] == RIGHT  && return (intent == RECV ? (localsize[i]+nhalo[i]+1:localsize[i]+2*nhalo[i]) : (localsize[i]+1:localsize[i]+nhalo[i]))
        region[i] == CENTER && return nhalo[i]+1:localsize[i]+nhalo[i]
        region[i] == ALL    && return 1:localsize[i]+2*nhalo[i]
    end

# ------------------------------------------------------------------------------
# HaloArray type and constructors
# ------------------------------------------------------------------------------

"""
    HaloArray{T, N, NHALO, SIZE, A} <: DenseArray{T, N}

Distributed dense array with local halo cells.

The parent storage has size `SIZE .+ 2 .* NHALO`. The logical `size(a)` is
`SIZE`, which describes only the interior region. Scalar indexing is shifted by
`nhalo(a)`, so `a[1, ...]` refers to the first interior point and values such as
`a[0, ...]` or `a[-1, ...]` can access halo cells when they still lie inside the
parent storage.

Type parameters:

- `T`: element type.
- `N`: number of Cartesian dimensions.
- `NHALO`: tuple of symmetric halo widths.
- `SIZE`: tuple containing the local interior size.
- `A`: parent dense array type.

Non-blocking communication uses one `MPI.MultiRequest` per Cartesian
dimension (`reqs(a, dim)`). `MultiRequest` is MPI.jl's batched-request
primitive: it owns a contiguous `Vector{MPI_Request}` whose slots are written
in place by `Irecv!`/`Isend`, and `Testall`/`Waitall` operate on that array
directly. The wait path therefore allocates nothing in steady state, and there
is no separate request cache to keep in sync.

The field `economic` records the construction option so that `similar(a)`
preserves halo-exchange behavior.
"""
struct HaloArray{T, N, NHALO, SIZE, A<:DenseArray{T, N}} <: DenseArray{T, N}
           data::A
        # ! buffers::Dict{Tuple{NTuple{N, Region}, Intent}, MPI.Buffer{A}}
        buffers::Dict{Tuple{NTuple{N, Region}, Intent}, MPI.Buffer}
    haloregions::Vector{NTuple{N, Region}}
           comm::MPI.Comm
           reqs::NTuple{N, MPI.MultiRequest}
       economic::Bool
    """
        HaloArray{T}(comm::MPI.Comm,
                localsize::NTuple{N, Int},
                    nhalo::NTuple{N, Int};
                 economic::Bool=true) where {T}

    Construct a `N` dimensional distributed array over the MPI communicator `comm`,
    with halo regions along the `N` cartesian dimension of size specified by `nhalo`.
    The local size of the array, i.e. excluding the halo points, is defined by
    `localsize`. The global size of the distributed array is implicit by the
    topology of the communicator, which must have a Cartesian topology of
    dimension `N`.

    `economic` is a Boolean communication policy, not a stencil-width or
    derivative-order parameter. It should be chosen explicitly by the caller or
    by any higher-level distributed-grid wrapper. Use `economic=true` when the
    algorithm only needs axis-aligned face halos. Use `economic=false` when
    edge or corner halo values must be valid after the exchange. All-dimension
    non-blocking swaps are rejected for `economic=false` when multiple exchange
    dimensions are active; use staged dimension swaps in that case.

    The practical choice is determined by the stencil footprint:

    ```text
    economic=true, for face-only stencils

        .   x   .
        x   u   x
        .   x   .

        u reads axis neighbours only. The four face halos must be exchanged;
        diagonal corner halos are not needed.

    economic=false, for corner-reading stencils

        x   x   x
        x   u   x
        x   x   x

        u may read diagonal neighbours. Corner data must propagate through
        staged exchanges, so messages carry transverse halo values.
    ```

    The communicator is shared with the caller; `HaloArray` does not own or free it.
    """
    function HaloArray{T}(comm::MPI.Comm,
                     localsize::NTuple{N, Int},
                         nhalo::NTuple{N, Int}; economic::Bool=true) where {N, T}

        # nhalo should be non-negative
        minimum(nhalo) ≥ 0 ||
            throw(ArgumentError("invalid halo specification"))

        # check N is right
        MPI.Cartdim_get(comm) == N ||
            throw(ArgumentError("incompatible communicator cartesian topology"))

        # obtain cartesian nprocesses information and convert to N tuples
        stuff = MPI.Cart_get(comm)
        nprocesses = tuple( Int.(stuff[1])...)
        isperiodic = tuple(Bool.(stuff[2])...)

        # if nhalo == 0, nprocesses can't be periodic and have more then one proc
        for dim = 1:N
            if nhalo[dim] == 0 && (isperiodic[dim] || nprocesses[dim] != 1)
                throw(ArgumentError("invalid halo specification"))
            end
        end

        # size of actual array include internal region size plus twice the halo
        data = zeros(T, localsize .+ nhalo .+ nhalo)

        # Build a dict of MPI buffers, one per (halo region, intent) pair.
        # Note that a HaloArray may carry halo regions that are not exchanged
        # with neighbours, but are used as ghost points for boundary
        # calculations, e.g. at a wall.
        # TODO: fix the type instability.
        # The dict is currently typed as `MPI.Buffer` (the abstract supertype)
        # because `MPI.Buffer(view)` returns different concrete subtypes
        # depending on whether the view is contiguous in memory. Resolving the
        # type instability would require either using a single derived-datatype
        # path (which needs a finalizer to call `MPI.Types.free!`) or wrapping
        # the buffer constructor with `invoke`. Neither has been shown to be a
        # measurable win in profiles so far; revisit if it turns up.
        buffers = Dict{Tuple{NTuple{N, Region}, Intent}, MPI.Buffer}()
        haloregions = regions_to_swap(nprocesses, isperiodic, economic)
        for region in haloregions
            for intent in (SEND, RECV)
                sub = view(data, subarray_slices(localsize, nhalo, region, intent)...)
                buffers[(region, intent)] = MPI.Buffer(sub)
            end
        end

        # Preallocate one MPI.MultiRequest per Cartesian dimension, sized to
        # two slots per halo region of that dimension (one send + one recv).
        # `MultiRequest` owns the contiguous `Vector{MPI_Request}` that the C
        # calls read and write directly, so reuse across non-blocking
        # haloswap! calls is allocation-free.
        nreqs_per_dim = zeros(Int, N)
        for region in haloregions
            nreqs_per_dim[exchange_dim(region)] += 2
        end
        reqs = ntuple(d -> MPI.MultiRequest(nreqs_per_dim[d]), N)

        return new{T, N, nhalo,
                   localsize, typeof(data)}(data, buffers, haloregions,
                                            comm, reqs, economic)
    end

    """
        HaloArray{T}(comm::MPI.Comm,
               nprocesses::NTuple{N, Int},
               isperiodic::NTuple{N, Bool},
                localsize::NTuple{N, Int},
                    nhalo::NTuple{N, Int};
                 economic::Bool=true) where {T}

    Create a Cartesian topology with number of processes `nprocesses` along each
    Cartesian direction with periodicity defined by `isperiodic`, then construct
    a `HaloArray` on the resulting communicator.

    `prod(nprocesses)` must match `MPI.Comm_size(comm)`. The newly created
    Cartesian communicator is stored by the array and shared by arrays derived
    from it with `similar` or `copy`; `HaloArray` does not free communicators.

    See also [`HaloArray{T}(comm, localsize, nhalo)`](@ref) for the meaning of
    `economic`, `localsize`, and `nhalo`.
    """
    function HaloArray{T}(comm::MPI.Comm,
                    nprocesses::NTuple{N, Int},
                    isperiodic::NTuple{N, Bool},
                     localsize::NTuple{N, Int},
                         nhalo::NTuple{N, Int}; economic::Bool=true) where {N, T}
        # checks
        prod(nprocesses) == MPI.Comm_size(comm) ||
            throw(ArgumentError("incompatible nprocesses specification"))

        # construct new communicator
        new_comm = MPI.Cart_create(comm, nprocesses, periodic=isperiodic, reorder=false)

        return HaloArray{T}(new_comm, localsize, nhalo; economic=economic)
    end
end

# ------------------------------------------------------------------------------
# Public accessors
# ------------------------------------------------------------------------------

"""
    nhalo(a::HaloArray{T, N}) where {T, N}

Return a `N`-tuple of `Int`s with the size of the halo region along the `N`
cartesian directions of the array. The size of the halo region is the same
on the "left" and "right" boundaries of the array.

For example, if `nhalo(a) == (2,)`, the first interior element is addressed as
`a[1]`, while the two left halo cells are addressed as `a[0]` and `a[-1]`.
"""
nhalo(::HaloArray{T, N, NHALO}) where {T, N, NHALO} = NHALO

"""
    comm(a::HaloArray)

Return the MPI communicator underlying the data distribution of the array.

The communicator is shared. `HaloArray` does not take ownership of it or free it.
"""
comm(a::HaloArray) = a.comm

"""
    reqs(a::HaloArray)
    reqs(a::HaloArray, dim::Integer)

Return reusable request storage for non-blocking communication.

`reqs(a)` returns a tuple of one `MPI.MultiRequest` per Cartesian dimension.
`reqs(a, dim)` returns the `MPI.MultiRequest` reserved for one Cartesian
dimension.

Each `MultiRequest` is reused by non-blocking [`haloswap!`](@ref). Callers must
complete outstanding requests, for example with `MPI.Waitall(reqs(a, dim))`,
before starting another non-blocking swap for the same dimension.
"""
reqs(a::HaloArray) = a.reqs
reqs(a::HaloArray, dim::Integer) = a.reqs[checkdim(a, dim)]

"""
    origin(a::HaloArray{T, N}) where {T< N}

Return the parent-array indices corresponding to the first interior element.

For an array with `nhalo(a) == (1, 2, 3)`, `origin(a) == (2, 3, 4)`, meaning
that `a[1, 1, 1]` is stored at `parent(a)[2, 3, 4]`.
"""
origin(a::HaloArray) = nhalo(a) .+ 1

# ------------------------------------------------------------------------------
# Array interface
# ------------------------------------------------------------------------------

"""
    parent(a::HaloArray)

Return the dense local storage backing `a`.

The parent includes both interior and halo cells. Its size is
`size(a) .+ 2 .* nhalo(a)`, and parent indices are ordinary Julia one-based
indices.
"""
@inline Base.parent(a::HaloArray) = a.data

"""
    Base.size(a::HaloArray)

Return the size of the array `a`. This excludes the contribution of the halo points.

This is the local interior size. Consequently, Julia iteration APIs such as
`axes(a)` and broadcast operate over the interior region, while scalar
`getindex` and `setindex!` additionally accept halo indices when the shifted
index lies inside `parent(a)`.
"""
Base.size(::HaloArray{T, N, NHALO, SIZE}) where {T, N, NHALO, SIZE} = SIZE

"""
    IndexStyle(::Type{<:HaloArray})

Declare that `HaloArray` uses Cartesian indexing.
"""
Base.IndexStyle(::HaloArray) = Base.IndexCartesian()

"""
    similar(a::HaloArray)
    similar(a::HaloArray, ::Type{T})

Construct a new `HaloArray` with the same communicator, local interior size,
halo width, and `economic` option as `a`.

The communicator is shared with `a`; it is not duplicated. The parent storage is
newly allocated and initialized to zeros by the constructor.
"""
Base.similar(a::HaloArray{T}) where {T} = similar(a, T)
Base.similar(a::HaloArray, ::Type{T}) where {T} =
    HaloArray{T}(comm(a), size(a), nhalo(a); economic=a.economic)

"""
    copy(a::HaloArray)

Return a new `HaloArray` with the same metadata and a copy of the full parent
storage, including halo cells.

This differs from broadcasting assignment, which only covers the interior
iteration domain.
"""
Base.copy(a::HaloArray) = (b = similar(a); parent(b) .= parent(a); b)

"""
    Base.elsize(::Type{<:HaloArray{T}}) where {T}

Return the element size used by Julia's array and broadcast machinery.
"""
Base.elsize(::Type{<:HaloArray{T}}) where {T} = sizeof(T)

"""
    _offset(nhalo, idxs)

Shift logical halo-aware indices into parent-array indices.

For example, with `nhalo == (2,)`, logical index `-1` maps to parent index `1`,
logical index `1` maps to parent index `3`, and logical index `size(a, 1) + 2`
maps to the final right halo cell.
"""
@inline _offset(nhalo, idxs) = nhalo .+ idxs

"""
    getindex(a::HaloArray, idxs::Vararg{Int, N}) where {N}

Read one element using halo-aware logical indices.

The interior domain is indexed as usual from `1:size(a, d)`. Halo cells are
available at indices outside that range when `_offset(nhalo(a), idxs)` remains
inside `parent(a)`. Thus `size(a)` still describes only the interior, but scalar
indexing can address valid halo cells. Bounds are checked against the shifted
parent indices.
"""
Base.@propagate_inbounds @inline function Base.getindex(a::HaloArray{T, N},
                                                     idxs::Vararg{Int, N}) where {T, N}
    pidxs = _offset(nhalo(a), idxs)
    @boundscheck checkbounds(parent(a), pidxs...)
    @inbounds v = parent(a)[pidxs...]
    return v
end

"""
    setindex!(a::HaloArray, v, idxs::Vararg{Int, N}) where {N}

Store one element using halo-aware logical indices.

See [`getindex(::HaloArray, ::Vararg{Int})`](@ref) for the indexing convention.
"""
Base.@propagate_inbounds @inline function Base.setindex!(a::HaloArray{T, N},
                                                   v, idxs::Vararg{Int, N}) where {T, N}
    pidxs = _offset(nhalo(a), idxs)
    @boundscheck checkbounds(parent(a), pidxs...)
    @inbounds parent(a)[pidxs...] = v
    return v
end

# ------------------------------------------------------------------------------
# Broadcast support
# ------------------------------------------------------------------------------

const HAStyle = Broadcast.ArrayStyle{HaloArray}

"""
    BroadcastStyle(::Type{<:HaloArray})

Use a custom broadcast style so broadcasted operations with `HaloArray`
arguments allocate `HaloArray` results.
"""
Base.BroadcastStyle(::Type{<:HaloArray}) = HAStyle()

# define proper broadcasting behaviour
"""
    similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{HaloArray}}, ::Type{T})

Allocate a broadcast result matching the first `HaloArray` found in the
broadcast expression.
"""
Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{HaloArray}}, ::Type{T}) where {T} = similar(find_ha(bc), T)

"""
    find_ha(x)

Return the first `HaloArray` contained in a broadcast expression or argument
tuple.
"""
find_ha(bc::Broadcast.Broadcasted)    = find_ha(bc.args)
find_ha(args::Tuple)                  = find_ha(find_ha(args[1]), Base.tail(args))
find_ha(u::HaloArray, rest)           = u
find_ha(::Any, rest)                  = find_ha(rest)
find_ha(x)                            = x
find_ha(::Tuple{})                    = nothing

# ------------------------------------------------------------------------------
# Low-level pointer conversion
# ------------------------------------------------------------------------------

"""
    unsafe_convert(::Type{Ptr{T}}, a::HaloArray{T}) where {T}

Return a pointer to the parent storage.

This mirrors dense-array pointer conversion and is intended only for low-level
code that knows it is operating on the full local parent storage, including halo
cells.
"""
Base.unsafe_convert(::Type{Ptr{T}}, a::HaloArray{T}) where {T} =
    Base.unsafe_convert(Ptr{T}, parent(a))

# ------------------------------------------------------------------------------
# MPI topology helpers
# ------------------------------------------------------------------------------

"""
    source_dest_ranks(a::HaloArray{T, N}, region::NTuple{N, Region}) where {T, N}

Return ranks of neighbouring processors at region `region`.

The returned tuple is `(source_rank, dest_rank)` as provided by
`MPI.Cart_shift`. `region` must identify exactly one `LEFT` or `RIGHT`
direction.
"""
function source_dest_ranks(a::HaloArray{T, N}, region::NTuple{N, Region}) where {T, N}
    # Determine the direction (dimension) of the shift from the
    # specification of the region. This is the index of the first
    # and unique index of `region` which is either right or left.
    # The displacement argument to Cart_shift depends on whether
    # we consider the LEFT or RIGHT boundaries. When there is no
    # neighbour, MPI returns MPI_PROC_NULL which is typically a
    # negative number.
    direction = exchange_dim(region)
    disp = region[direction] == LEFT ? -1 : 1

    # call the MPI api and return the source and destination ranks
    return MPI.Cart_shift(comm(a), direction - 1, disp)
end

"""
    exchange_dim(region::NTuple{N, Region}) where {N}

Return the Cartesian dimension exchanged by `region`.

Each halo exchange region must contain exactly one `LEFT` or `RIGHT` entry. That
entry identifies the dimension shifted by `MPI.Cart_shift`.
"""
exchange_dim(region::NTuple{N, Region}) where {N} =
    findfirst(r -> r in (LEFT, RIGHT), region)

"""
    checkdim(a::HaloArray, dim::Integer)

Validate a Cartesian dimension index for `a` and return it as an `Int`.
"""
function checkdim(::HaloArray{T, N}, dim::Integer) where {T, N}
    1 <= dim <= N || throw(ArgumentError("invalid exchange dimension"))
    return Int(dim)
end

"""
    haloregions(a::HaloArray, dim::Integer)

Return the configured halo exchange regions for Cartesian dimension `dim`.
"""
haloregions(a::HaloArray, dim::Integer) =
    filter(region -> exchange_dim(region) == checkdim(a, dim), a.haloregions)

# ------------------------------------------------------------------------------
# Halo exchange
# ------------------------------------------------------------------------------

"""
    haloswap!(a::HaloArray{T, N}[, block::Bool=true]) where {T, N}
    haloswap!(a::HaloArray, dim::Integer, block::Bool=true)

Exchange halo data between adjacent ranks in the Cartesian communicator.

If `block` is `true`, communication is performed with `MPI.Sendrecv!` and the
call returns after the selected halo cells have been updated. If `block` is
`false`, communication is started with `MPI.Irecv!`/`MPI.Isend`; the internal
requests stored in [`HaloArray`](@ref) can be accessed with
[`reqs(::HaloArray)`](@ref).

In one dimension, with one halo cell on each side:

```text
before exchange

rank p-1                 rank p                   rank p+1
+----+---------+----+    +----+---------+----+    +----+---------+----+
| HL | interior| HR |    | HL | interior| HR |    | HL | interior| HR |
+----+---------+----+    +----+---------+----+    +----+---------+----+

for region = RIGHT on rank p

rank p sends its right interior strip ------------------------> rank p+1
rank p receives into its left halo <--------------------------- rank p-1

for region = LEFT on rank p

rank p sends its left interior strip  ------------------------> rank p-1
rank p receives into its right halo <-------------------------- rank p+1

after both LEFT and RIGHT regions have been exchanged

rank p
+----------------------+---------+----------------------+
| data from rank p-1   | interior| data from rank p+1   |
+----------------------+---------+----------------------+
```

In multiple dimensions, the same rule is applied to each configured
`Region` tuple. For example, `(LEFT, CENTER)` sends a vertical face to the left
neighbour and receives the matching right-neighbour data into
`opposite((LEFT, CENTER)) == (RIGHT, CENTER)`.

For a non-blocking swap, callers must complete the requests before reusing them,
for example with `foreach(MPI.Waitall, reqs(a))`.

All-at-once non-blocking swaps are rejected for `economic=false` when multiple
exchange dimensions are active. Corner/edge data must be propagated in a
defined order: one dimension is exchanged, those requests are completed, and the
next dimension can then send the halo values produced by the previous stage.

For `economic=false`, users should stage the non-blocking exchange explicitly:

```julia
r = haloswap!(a, 1, false)
MPI.Waitall(r)

r = haloswap!(a, 2, false)
MPI.Waitall(r)
```

The requests for each stage are also available with `reqs(a, dim)`. This lets
callers wait between dimensions and overlap each stage with useful work that
does not read the in-flight halo values.

For `economic=false`, staged exchanges propagate corner values by allowing each
later dimension to send halo data filled by earlier dimensions:

```text
2D economic=false, staged order dim 1 then dim 2

Legend:
    I  local interior
    H  halo not yet valid
    1  halo filled by the dimension-1 stage
    2  halo filled by the dimension-2 stage
    C  corner filled by the dimension-2 stage using data carried by dim 1

initial local storage

    +----+----+----+
    | H  | H  | H  |
    +----+----+----+
    | H  | I  | H  |
    +----+----+----+
    | H  | H  | H  |
    +----+----+----+

stage dim 1: exchange (LEFT, ALL) and (RIGHT, ALL), then wait

    +----+----+----+
    | 1  | 1  | 1  |
    +----+----+----+
    | H  | I  | H  |
    +----+----+----+
    | 1  | 1  | 1  |
    +----+----+----+

stage dim 2: exchange (ALL, LEFT) and (ALL, RIGHT), then wait

    +----+----+----+
    | C  | 1  | C  |
    +----+----+----+
    | 2  | I  | 2  |
    +----+----+----+
    | C  | 1  | C  |
    +----+----+----+

The dim-2 messages use `ALL` in dim 1, so they include values produced by the
dim-1 stage. That is how diagonal corner data moves through two
one-dimensional exchanges.
```

Blocking swaps return `nothing`. Non-blocking swaps return the
`MPI.MultiRequest` for a staged swap, or the tuple of `MultiRequest`s for an
all-dimension swap.
"""
function haloswap!(a::HaloArray, block::Bool=true)
    if !block && !a.economic && length(a.haloregions) > 2
        throw(ArgumentError("non-blocking halo swaps are not supported with economic=false for multiple exchange dimensions"))
    end
    return block ? _haloswap!(a) : _Ihaloswap!(a)
end

function haloswap!(a::HaloArray, dim::Integer, block::Bool=true)
    return block ? _haloswap!(a, dim) : _Ihaloswap!(a, dim)
end

"""
    _haloswap!(a::HaloArray)
    _haloswap!(a::HaloArray, dim::Integer)

Perform a blocking halo exchange using `MPI.Sendrecv!`.

Each configured halo region is sent to the destination rank returned by
[`source_dest_ranks`](@ref), while the opposite receive region is filled from
the corresponding source rank.
"""
_haloswap!(a) = _haloswap!(a, a.haloregions)
_haloswap!(a, dim::Integer) = _haloswap!(a, haloregions(a, dim))

function _haloswap!(a, regions)
    # We do not need to care about the periodicity and take care of boundary conditions
    # see https://www.mpi-forum.org/docs/mpi-1.1/mpi-11-html/node53.html
    #
    # For every halo region, we need to SEND it to the adjacent processor and
    # RECV from the processor adjacent to the opposite region. This corresponds
    # to a standard shifting of the data across processors and avoids deadlock.
    for region in regions
        tag = findfirst(==(region), a.haloregions)
        source_rank, dest_rank = source_dest_ranks(a, region)
        MPI.Sendrecv!(a.buffers[(         region,  SEND)],
                      a.buffers[(opposite(region), RECV)],
                      comm(a),
                      dest=dest_rank,
                      sendtag=tag,
                      source=source_rank,
                      recvtag=tag)
    end
    return nothing
end

"""
    _Ihaloswap!(a::HaloArray)
    _Ihaloswap!(a::HaloArray, dim::Integer)

Start a non-blocking halo exchange using `MPI.Irecv!` and `MPI.Isend`.

The requests are stored in `reqs(a)` for all-dimension swaps and in
`reqs(a, dim)` for staged swaps. They must be completed by the caller before the
same request vector is reused.
"""
function _Ihaloswap!(a)
    for dim in eachindex(reqs(a))
        _Ihaloswap!(a, dim)
    end
    return reqs(a)
end

function _Ihaloswap!(a, dim::Integer)
    # Same as haloswap! except using non-blocking communication.
    requests = reqs(a, dim)
    for (i, region) in enumerate(haloregions(a, dim))
        tag = findfirst(==(region), a.haloregions)
        source_rank, dest_rank = source_dest_ranks(a, region)
        MPI.Irecv!(a.buffers[(opposite(region), RECV)],
                   comm(a),
                   requests[2*i],
                   source=source_rank,
                   tag=tag)
        MPI.Isend(a.buffers[(region, SEND)],
                  comm(a),
                  requests[2*i-1],
                  dest=dest_rank,
                  tag=tag)
    end
    return requests
end

end
