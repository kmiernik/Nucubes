using Nucubes
using HDF5
using Printf
using BenchmarkTools
using Profile
using Test
using Distributed
using InteractiveUtils


function create_dummy(file_name, n)
    println("* Creting test file ", file_name)
    h5open(file_name, "w") do fid
        g = create_group(fid, "GeGeGe")
        chunk_size = 10_000
        i = 0
        data = Array{UInt32, 2}(undef, (7, n + 2))
        while i < n
            for loc in [[1, 2, 3], [1, 3, 2], [2, 3, 1], 
                        [2, 1, 3], [3, 1, 2], [3, 2, 1]]
                i += 1
                # 1-2-3 "prompt"
                data[1, i] = loc[1]
                data[2, i] = loc[2]
                data[3, i] = loc[3]
                data[4, i] = 1
                data[5, i] = 1
                data[6, i] = 1
            end
        end

        # 1-2-10 "prompt" (out-of-range)
        data[1, end-1] = 1
        data[2, end-1] = 2
        data[3, end-1] = 10
        data[4, end-1] = 1
        data[5, end-1] = 1
        data[6, end-1] = 1

        # 1-2-3 "delayed"
        data[1, end] = 1
        data[2, end] = 2
        data[3, end] = 3
        data[4, end] = 3
        data[5, end] = 3
        data[6, end] = 3

        for m in 3:3
            s = @Printf.sprintf("%d", m)
            for i in 1:n+2
                data[7, i] = m * 100 + 7
            end
            dout = create_dataset(g, s, datatype(UInt32), dataspace(7, n+2),
                                  chunk=(7, chunk_size))
            dout[:, :] = data[:, :]
        end
    end
end

function benchmark(test_file)
    println("* Starting benchmarks ")
    g = Nucubes.GGGate([8], [1.0, 2.0], [2.0, 3.0], [1.0, 2.0], [3.0, 4.0],
                        "ppp", 1, 1)
    b = @benchmark Nucubes.process($test_file , [3, 3], $g, Nucubes.select_ggg!)
    println(b)

    #gates = [g for i in 1:6]
    #b = @benchmark pmap(gate->Nucubes.process($test_file, [3, 3], gate,
    #                                Nucubes.select_ggg!), $gates)
    #println(b)
    println()
end

function profile(test_file)
    println("* Starting profiling ")
    g = Nucubes.GGGate([8], [1.0, 2.0], [2.0, 3.0], [1.0, 2.0], [3.0, 4.0],
                        "ppp", 1, 1)
    # Dry run
    Nucubes.process(test_file, [3, 3], g, Nucubes.select_ggg!)
    # Actual run
    @profile Nucubes.process("dummy.h5" , [3, 3], g, Nucubes.select_ggg!)
    Profile.print()
    println()
end

function test_results(test_file, n)
    println("* Results testing ")
    @testset "Nucubes" begin
        g1 = Nucubes.GGGate([8], [1.0, 2.0], [2.0, 3.0], [1.0, 2.0], [3.0, 4.0],
                            "ppp", 1, 1)
        m = Nucubes.process(test_file , [3, 3], g1, Nucubes.select_ggg!) 
        @test m[4] == n
        g2 = Nucubes.GGGate([8], [1.0, 2.0], [2.0, 3.0], [1.0, 2.0], [3.0, 4.0],
                            "ddd", 1, 1)
        m = Nucubes.process(test_file , [3, 3], g2, Nucubes.select_ggg!)
        @test m[4] == 1
        @InteractiveUtils.code_warntype Nucubes.process(test_file, [3, 3], g2, Nucubes.select_ggg!)
    end
    println()
end

function types_warnings(test_file)
    println("* Types warnings ")
    println("  + GGGate ")
    @InteractiveUtils.code_warntype Nucubes.GGGate([8], [1.0, 2.0], [2.0, 3.0], [1.0, 2.0], [3.0, 4.0],
                                    "ddd", 1, 1)
    g = Nucubes.GGGate([8], [1.0, 2.0], [2.0, 3.0], [1.0, 2.0], [3.0, 4.0],
                                    "ddd", 1, 1)
    println("  + process ")
    @InteractiveUtils.code_warntype Nucubes.process(test_file, [3, 3], g, Nucubes.select_ggg!)
    println()
end

println("# Start tests ")
test_file = "dummy.h5"
create_dummy(test_file, 12_000_000)

benchmark(test_file)
profile(test_file)
test_results(test_file, 12_000_000)
types_warnings(test_file)

println("* Removing test file ", test_file)
rm(test_file, force=true)
println("# Finish tests ")
