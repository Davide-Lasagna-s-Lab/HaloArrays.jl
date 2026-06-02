"""
    HaloArrays

MPI-distributed dense arrays with halo cells.

`HaloArray` stores a local interior array together with symmetric halo cells on
both sides of each Cartesian direction. The public `size` and default Julia
iteration domain describe only the interior region. Scalar indexing is
halo-aware: indices such as `0` or negative values are accepted when, after
shifting by `nhalo(a)`, they still refer to storage inside `parent(a)`.

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
"""
@enum Region LEFT RIGHT CENTER ALL

"""
    swapregions(nprocesses::NTuple{N, Int},
                isperiodic::NTuple{N, Bool},
                economic::Bool) where {N}

Construct an array of `N`-tuples of the enum type `Region` specifying the
regions where halo exchange between processors or within the same processor
should happen. The exchange should happen along directions `dim` where
`nprocesses[dim] > 1` and `isperiodic[dim]` is true. At most three exchange
directions are currently supported.

If `economic` is `true`, only face regions are exchanged and transverse
directions are restricted to `CENTER`. If `economic` is `false`, transverse
directions are set to `ALL`, so edge and corner values are included in the
messages.
"""
function swapregions(nprocesses::NTuple{N, Int},
                     isperiodic::NTuple{N, Bool}, economic::Bool) where {N}
    # vector of dimensions along which exchange must happen
    exdims = findall( (nprocesses .!= 1) .| isperiodic )

    # number of nprocesses over which exchange must happen
    ndims = length(exdims)

    # we do not expect 4D transfers
    ndims > 3 && throw(ArgumentError("invalid argument"))

    # init
    regions = NTuple{N, Region}[]

    # efficient transfers of areas
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

    `economic=true` exchanges only face regions and avoids corner/edge values.
    `economic=false` exchanges wider regions that include transverse halo
    values. All-dimension non-blocking swaps are rejected for `economic=false`
    when multiple exchange dimensions are active; use staged dimension swaps in
    that case.

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
        #
        # The dict is currently typed as `MPI.Buffer` (the abstract supertype)
        # because `MPI.Buffer(view)` returns different concrete subtypes
        # depending on whether the view is contiguous in memory. Resolving the
        # type instability would require either using a single derived-datatype
        # path (which needs a finalizer to call `MPI.Types.free!`) or wrapping
        # the buffer constructor with `invoke`. Neither has been shown to be a
        # measurable win in profiles so far; revisit if it turns up.
        buffers = Dict{Tuple{NTuple{N, Region}, Intent}, MPI.Buffer}()
        haloregions = swapregions(nprocesses, isperiodic, economic)
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
    haloswap!(a::HaloArray, dim::Integer; block::Bool=true)

Executes a blocking halo swap between adjacent processes. If `block` is `true`
then the communication is blocking, for `false` the communication is non-
blocking which updates the internal requests stored in [`HaloArray`](@ref) which
can be accessed with [`reqs(::HaloArray)`](@ref).

For a non-blocking swap, callers must complete the requests before reusing them,
for example with `foreach(MPI.Waitall, reqs(a))`.

Non-blocking swaps are rejected for `economic=false` when multiple exchange
dimensions are active, because those swaps include corner/edge data whose update
ordering is not currently implemented safely.

Use `haloswap!(a, dim; block=false)` to stage a non-blocking exchange in a
single Cartesian dimension. The requests for that stage are available with
`reqs(a, dim)`. This lets callers wait between dimensions and overlap each
stage with useful work.

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

function haloswap!(a::HaloArray, dim::Integer; block::Bool=true)
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
