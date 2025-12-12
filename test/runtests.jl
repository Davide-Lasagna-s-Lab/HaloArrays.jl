using MPI

# list of filenames and number of processors 
# use `0` for serial code
tests = [
          ("test_utils.jl",     0),
          ("test_haloarray.jl", 2),
          ("test_sum.jl",       1),
          ("test_broadcast.jl", 1),
          ("test_swap.jl",      4)
    ]

for (filename, nprocs) in tests
    if nprocs == 0
        run(`julia --startup-file=no $filename`)
        Base.with_output_color(:green, stdout) do io
            println(io, "\tSUCCESS: $filename - with $nprocs processors")
        end
    else
        try
            run(`mpiexecjl -n $nprocs --project=./ julia --startup-file=no $filename`)
            # using this form of call (recommended by MPI.jl) lows down the sum and broadcasting
            # tests significantly, as well as causing them to fail, for unknown reasons
            # run(`$(mpiexec()) -np $nprocs $(Base.julia_cmd()) $filename`)
            Base.with_output_color(:green, stdout) do io
                println(io, "\tSUCCESS: $filename - with $nprocs processors")
            end
        catch ex
            Base.with_output_color(:red, stderr) do io
                println(io, "\tERROR: $filename - with $nprocs processors")
                showerror(io, ex, backtrace())
            end
        end
    end
end
