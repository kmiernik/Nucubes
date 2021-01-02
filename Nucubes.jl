"""
Module for nucubes analysis
"""
module Nucubes

using HDF5
using Dates
using Printf
using DataFrames
using StatsBase

export gegege
export gate


"""
Version with dataframes and histogram (not used)
"""
function gate_df!(data::Array{UInt32, 2}, 
                  y0::Int64, y1::Int64, 
                  z0::Int64, z1::Int64, 
                  E_unit::Int64, matrix::Array{Int64, 1})
    df = DataFrame(data', 
                    ["E0", "E1", "E2", "t0", "t1", "t2", "pattern"])
    D = size(matrix)[1]
    edges = [0:D;]
    for loc in [[1, 2, 3], [1, 3, 2], [2, 3, 1], 
                                    [2, 1, 3], [3, 1, 2], [3, 2, 1]]
        sub = df[
            (y0 .<= df[:, loc[2]] .< y1) .& (z0 .<= df[:, loc[1]] .< z1)
                , :]
        h = fit(Histogram, sub[:, loc[3]] ./ E_unit, edges)
        matrix += h.weights
    end
end

"""
Version with data smart indexing and histogram (not used)
"""
function gate_filter!(data::Array{UInt32, 2}, 
              y0::Int64, y1::Int64, 
              z0::Int64, z1::Int64, 
              E_unit::Int64, matrix::Array{Int64, 1})
    D = size(matrix)[1]
    edges = [0:D;]
    for loc in [[1, 2, 3], [1, 3, 2], [2, 3, 1], 
                                    [2, 1, 3], [3, 1, 2], [3, 2, 1]]

        sub = data[:, 
            (y0 .<= data[loc[2], :] .< y1) .& (z0 .<= data[loc[1], :] .< z1)]
        h = fit(Histogram, sub[:, loc[3]] ./ E_unit, edges)
        matrix += h.weights
    end

end


"""
Fastest version with simple iteration
c is conditions array:
    * y0, y1, z0, z1, t00, t01, t10, t11, t02, t02
       1  2   3   4    5    6   7    8    9    10
"""
function gate_iter!(data::Array{UInt32, 2}, c::Array{Int64, 1}, 
                    E_unit::Int64, matrix::Array{Int64, 1})
    @assert size(c)[1] >= 10
    D = size(matrix)[1]
    n = size(data)[2]
    for loc in [[1, 2, 3], [1, 3, 2], [2, 3, 1], 
                                    [2, 1, 3], [3, 1, 2], [3, 2, 1]]
        for i in 1:n
            if (   (c[1] <= data[loc[1], i] < c[2]) 
                && (c[3] <= data[loc[2], i] < c[4]) 
                && (c[5] <= data[loc[1]+3, i] <= c[6])
                && (c[7] <= data[loc[2]+3, i] <= c[8])
                && (c[9] <= data[loc[3]+3, i] <= c[10])
               )
                k = trunc(Int64, data[loc[3], i] / E_unit) + 1
                if k <= D
                    matrix[k] += 1
                end
            end
        end
    end
end


"""
* Usage is similar to that of nucubes.py except of
* gate_m: if empty list, a whole range is used, if length is 1, it is 
          taken as lower limit, if length is 2 - low and high end, and
          any other, as a list of multiplicities
"""        
function gegege(fin::HDF5File, gate_z::Array{Float64, 1}, 
                gate_y::Array{Float64, 1}, 
                prompt::Array{Float64, 1},
                delayed::Array{Float64, 1},
                gate_m::Array{Int64, 1}, 
                ttype::String)::Array{Int64}
    D = 4096
    E_unit = 100
    t_unit = 1000

    conditions = [round(Int, gate_y[1] * E_unit), 
                  round(Int, gate_y[2] * E_unit),
                  round(Int, gate_z[1] * E_unit),
                  round(Int, gate_z[2] * E_unit)]
    @assert length(ttype) >= 3
    for c in ttype
        if c == 'p'
            push!(conditions, Int64(prompt[1] * t_unit))
            push!(conditions, Int64(prompt[2] * t_unit))
        elseif c == 'd'
            push!(conditions, Int64(delayed[1] * t_unit))
            push!(conditions, Int64(delayed[2] * t_unit))
        else
            push!(conditions, Int64(0))
            push!(conditions, Int64(400 * t_unit))
        end
    end

    M = Int64[]
    group::HDF5Group = fin["GeGeGe"]
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
        dataset::HDF5Dataset = group[string(m)]
        n_all += size(dataset)[2]
    end

    t0 = Dates.Time(Dates.now())
    
    i_report = 0
    n_report = 1_000_000
    for m in multi
        dataset::HDF5Dataset = group[string(m)]
        n = size(dataset)[2]
        chunk_size = min(n, 10_000)
        try
            chunk_size = HDF5.get_chunk(dataset)[2]
        catch
        end
        left_pos = 1
        is_something_left = true
        while is_something_left
            right_pos = left_pos + chunk_size - 1
            if right_pos > n
                right_pos = n
                is_something_left = false
            elseif left_pos >= right_pos
                break
            end
            data = dataset[:, left_pos:right_pos]

            gate_iter!(data, conditions, E_unit, matrix)

            n_processed += right_pos - left_pos

            j_report = trunc(Int64, n_processed / n_report)
            if j_report > i_report
                i_report += 1
                t1 = Dates.Time(Dates.now())
                dt = t1 - t0
                print("\r", round(n_processed / n_all * 100, digits=1), "% ",
                    round(dt.value * 1e-9, digits=1), " s (",
                    round(dt.value * 1e-9 * n_all / n_processed, digits=0),
                    " s)          ")
            end

            left_pos = right_pos + 1
        end
    end
    println()
    t1 = Dates.Time(Dates.now())
    dt = t1 - t0
    @printf("%12d %8.5f\n", n_all, dt.value * 1.0e-9)
    return matrix
end


"""
Check time distribution of events
"""
function time_distribution(fin::HDF5File, 
                           gate_m::Array{Int64, 1})::Array{Int64, 2}
    DE = 4096
    DT = 600
    matrix = zeros(Int64, DE, DT)
    t_unit = 1000
    E_unit = 100

    M = Int64[]
    group::HDF5Group = fin["GeGeGe"]
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
   
    n_all = 0
    n_processed = 0
    for m in multi
        dataset::HDF5Dataset = group[string(m)]
        n_all += size(dataset)[2]
    end

    t0 = Dates.Time(Dates.now())
    
    i_report = 0
    n_report = 1_000_000
    for m in multi
        dataset::HDF5Dataset = group[string(m)]
        n = size(dataset)[2]
        chunk_size = min(n, 10_000)
        try
            chunk_size = HDF5.get_chunk(dataset)[2]
        catch
        end
        left_pos = 1
        is_something_left = true
        while is_something_left
            right_pos = left_pos + chunk_size - 1
            if right_pos > n
                right_pos = n
                is_something_left = false
            elseif left_pos >= right_pos
                break
            end

            data = dataset[:, left_pos:right_pos]
            
            n_chunk = size(data)[2]
            for i in 1:n_chunk
                for j in 1:3
                    iE = trunc(Int64, data[j, i] / E_unit) + 1
                    iT = trunc(Int64, data[j+3, i] / t_unit) + 1
                    if 1 <= iE <= DE && 1 <= iT <= DT
                        matrix[iE, iT] += 1
                    end
                end
            end

            n_processed += right_pos - left_pos

            j_report = trunc(Int64, n_processed / n_report)
            if j_report > i_report
                i_report += 1
                t1 = Dates.Time(Dates.now())
                dt = t1 - t0
                print("\r", round(n_processed / n_all * 100, digits=1), "% ",
                    round(dt.value * 1e-9, digits=1), " s (",
                    round(dt.value * 1e-9 * n_all / n_processed, digits=0),
                    " s)          ")
            end

            left_pos = right_pos + 1
        end
    end
    println()
    t1 = Dates.Time(Dates.now())
    dt = t1 - t0
    @printf("%12d events in %8.3f s\n", n_all, dt.value * 1.0e-9)
    return matrix

end


end
