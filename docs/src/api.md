# API Reference

## Public API

```@docs
HaloArrays
HaloArray
nhalo
comm
reqs
origin
haloswap!
Base.parent(::HaloArray)
Base.size(::HaloArray)
Base.similar(::HaloArray)
Base.copy(::HaloArray)
```

## Internal Region Helpers

These helpers are mostly useful when reading or extending the implementation.

```@docs
HaloArrays.Region
HaloArrays.Intent
HaloArrays.regions_to_swap
HaloArrays.opposite
HaloArrays.subarray_slices
HaloArrays.source_dest_ranks
```
