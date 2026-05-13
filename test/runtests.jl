const TEST_DIR = @__DIR__
const PACKAGE_DIR = dirname(TEST_DIR)
const LOAD_PATH_SEPARATOR = Sys.iswindows() ? ';' : ':'
const TEST_LOAD_PATH = join((TEST_DIR, PACKAGE_DIR, "@stdlib"), LOAD_PATH_SEPARATOR)
const TEST_JULIA = `$(Base.julia_cmd()) --startup-file=no --project=$(TEST_DIR)`

using Pkg
Pkg.activate(TEST_DIR)
Pkg.instantiate()

using MPI

testcmd(filename) =
    addenv(`$(TEST_JULIA) $(joinpath(TEST_DIR, filename))`,
           "JULIA_LOAD_PATH" => TEST_LOAD_PATH)

mpitestcmd(filename, nprocs) =
    addenv(`$(mpiexec()) -n $nprocs $(TEST_JULIA) $(joinpath(TEST_DIR, filename))`,
           "JULIA_LOAD_PATH" => TEST_LOAD_PATH)

# list of filenames and number of processors 
# use `0` for serial code
const tests = [
          ("test_utils.jl",     0),
          ("test_haloarray.jl", 2),
          ("test_broadcast.jl", 1),
          ("test_swap.jl",      4)
    ]

for (filename, nprocs) in tests
    if nprocs == 0
        run(testcmd(filename))
        Base.with_output_color(:green, stdout) do io
            println(io, "\tSUCCESS: $filename - with $nprocs processors")
        end
    else
        run(mpitestcmd(filename, nprocs))
        # using this form of call (recommended by MPI.jl) lows down the sum and broadcasting
        # tests significantly, as well as causing them to fail, for unknown reasons
        # run(`$(mpiexec()) -np $nprocs $(Base.julia_cmd()) $filename`)
        Base.with_output_color(:green, stdout) do io
            println(io, "\tSUCCESS: $filename - with $nprocs processors")
        end
    end
end
