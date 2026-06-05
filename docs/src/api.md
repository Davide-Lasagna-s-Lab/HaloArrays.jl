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

## Exchange Region Enums

```@docs
HaloArrays.Region
HaloArrays.LEFT
HaloArrays.RIGHT
HaloArrays.CENTER
HaloArrays.ALL
HaloArrays.Intent
HaloArrays.SEND
HaloArrays.RECV
```

Implementation-oriented helpers are covered in the
[Internals](manual/internals.md) page.
