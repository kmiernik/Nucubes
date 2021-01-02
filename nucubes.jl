"""
This is test of julia's speed of Ge-Ge-Ge gating
Unfortunately it looks like bottleneck is the HDF loading
and casting to DataFrame, speed is comparable of python's version
"""

using HDF5
using Dates
using DataFrames
using StatsBase
using Printf

"""
* Usage is similar to that of nucubes.py except of
* gate_m: if empty list, a whole range is used, if length is 1, it is 
          taken as lower limit, if length is 2 - low and high end, and
          any other, as a list of multiplicities
"""        
function gegege(fin::HDF5File, gate_z::Array{Float64}, gate_y::Array{Float64}, 
                gate_m::Array{Int64})::Array{Int64}
    D = 4096
    edges = [0:D;]
    E_unit = 100
    t_unit = 1000

    M = Int64[]
    group = fin["GeGeGe"]
    for m_set in names(group)
        append!(M, parse(Int64, m_set))
    end

    if size(gate_m)[1] == 0
        multi = M
    elseif size(gate_m)[1] == 1
        multi = [max(gate_m[1], first(M)):last(M);]
    elseif size(gate_m)[1] == 2
        low_m = max(gate_m[1], first(M))
        if low_m > last(M)
            low_m = last(M)
        end
        high_m = min(gate_m[2], last(M))
        multi = [low_m:high_m;]
    else
        multi = Int64[]
        for m in gate_m
            if first(M) < m < last(M)
                append!(multi, m)
            end
        end
    end
   
    matrix = zeros(Int64, D)
    n_all = 0
    n_processed = 0
    for m in multi
        dataset = group[string(m)]
        n_all += size(dataset)[2]
    end

    t0 = Dates.Time(Dates.now())

    for m in multi
        dataset = group[string(m)]
        chunk_size = get_chunk(dataset)[2]
        n = size(dataset)[2]
        left_pos = 1
        is_something_left = true
        while is_something_left
            right_pos = left_pos + chunk_size - 1
            if right_pos > n
                right_pos = n
                is_something_left = false
            end
            # Here is the bottleneck
            data = dataset[:, left_pos:right_pos]

            df = DataFrame(data', 
                           ["E0", "E1", "E2", "t0", "t1", "t2", "pattern"])

            # Multithreading has some small improvement in speed
            Threads.@threads for loc in [[1, 2, 3], [1, 3, 2], [2, 3, 1], 
                                         [2, 1, 3], [3, 1, 2], [3, 2, 1]]
                sub = df[(gate_y[1] .<= df[:, loc[2]] ./ E_unit .< gate_y[2]) .& (gate_z[1] .<= df[:, loc[1]] ./ E_unit .< gate_z[2]), :]
                h = fit(Histogram, sub[:, loc[3]] ./ E_unit, edges)
                matrix += h.weights
            end

            n_processed += right_pos - left_pos

            t1 = Dates.Time(Dates.now())
            dt = t1 - t0
            print("\r", round(n_processed / n_all * 100, digits=1), "% ",
                  round(dt.value * 1e-9, digits=1), " s (",
                  round(dt.value * 1e-9 * n_all / n_processed, digits=0), " s)")

            left_pos = right_pos + 1
        end
    end
    println()
    t1 = Dates.Time(Dates.now())
    dt = t1 - t0
    @printf("%12d %8.5f\n", n_all, dt.value * 1.0e-9)
    return matrix
end

# Notice that this is for test purposes only
fin = h5open("u238_m_ggg.h5", "r")
matrix = gegege(fin, [241.7, 243.7], [102.1, 104.1], [])
println(matrix)
