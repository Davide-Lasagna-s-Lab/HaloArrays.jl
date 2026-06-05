# Halo Indexing

A `HaloArray` separates the public array domain from the padded parent storage.

```julia
size(a)       # local interior size
axes(a)       # local interior axes
nhalo(a)      # halo width tuple
parent(a)     # dense local storage, including halos
origin(a)     # parent index of a[1, 1, ...]
```

For a one-dimensional array with `size(a) == (5,)` and `nhalo(a) == (2,)`:

```text
logical index:   -1   0 | 1   2   3   4   5 | 6   7
storage role:    HL  HL | I   I   I   I   I | HR  HR
parent index:     1   2 | 3   4   5   6   7 | 8   9

HL, HR: left/right halo storage
I:      local interior storage
```

Scalar indexing shifts logical indices by the halo width before reading or
writing `parent(a)`. This is why the first interior point is `a[1]`, while
left halo values are available at `a[0]`, `a[-1]`, and so on when the halo width
allows it.

```julia
a[1, j]              # first interior point in dimension 1
a[0, j]              # first left halo point
a[-1, j]             # second left halo point, if nhalo(a)[1] >= 2
a[size(a, 1) + 1, j] # first right halo point
```

Broadcast and iteration remain interior-only:

```julia
a .= 1.0

for I in CartesianIndices(a)
    # visits interior points only
end
```

Use [`parent`](@ref) when you deliberately want the complete local storage:

```julia
fill!(parent(a), 0.0) # clears interior and halos
```

`copy(a)` copies the full parent storage, including halo cells. Broadcast
assignment such as `b .= a` follows the interior iteration domain.

