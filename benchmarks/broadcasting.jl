# Compare performance of broadcast addition of HaloArray with standard dense arrays, and HaloArray with HaloArray.

using BenchmarkTools, MPI

using HaloArrays

MPI.Init()

comm = MPI.COMM_WORLD

# define constants
Ny    = 16
nhalo = 2
Nx    = 33
Nz    = 33
Nt    = 33
P     = 1

# construct arrays
h1  = HaloArray{Float64}(comm, (1, P, 1, 1), (false, false, false, false), (Ny, Nx, Nz, Nt), (0, nhalo, 0, 0))
h2  = similar(h1)
h2 .= randn(Float64, size(h2))
h3  = similar(h1)
h3 .= randn(Float64, size(h3))

d1 = zeros(Float64, div(Ny, P), Nx, Nz, Nt)
d2 = randn(Float64, div(Ny, P), Nx, Nz, Nt)
d3 = randn(Float64, div(Ny, P), Nx, Nz, Nt)

# broadcasting wrapper function
f(a, b, c) = (a .= 2.0.*b .+ c./5.0; a)

# check addition correctness
f(d1, h2, d3)
# @assert d1 == h2[1:end, :, :, :]*2 + d3/5

# benchmark
bm_dd = @benchmark f(d1, d2, d3)
bm_hh = @benchmark f(h1, h2, h3)
bm_hd = @benchmark f(d1, d2, h3)

println("Broadcast: DA & DA")
show(stdout, MIME("text/plain"), bm_dd)
println(); println()
println("Broadcast: HA & HA")
show(stdout, MIME("text/plain"), bm_hh)
println(); println()
println("Broadcast: HA & DA")
show(stdout, MIME("text/plain"), bm_hd)
println()
