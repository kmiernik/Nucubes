"""
Module for nugates analysis
Version with conditions and process functions
"""
module Nucubes

using HDF5
using Dates
using Printf
using DataFrames
using StatsBase
using BenchmarkTools

export process
export Gate
export GGGate
export GGTate
export TGate
export select_ggg!
export select_ggt!
export select_t!


abstract type Gate
end

struct GGGate <: Gate
    D::Array{Int64, 1}
    y1::Int64
    y2::Int64
    z1::Int64
    z2::Int64
    t11::Int64
    t12::Int64
    t21::Int64
    t22::Int64
    t31::Int64
    t32::Int64
    E_unit::Int64
    t_unit::Int64
end


struct GGTGate <: Gate
    D::Array{Int64, 2}
    y1::Int64
    y2::Int64
    z1::Int64
    z2::Int64
    E_unit::Int64
    t_unit::Int64
end


struct TGate <: Gate
    D::Array{Int64, 2}
    E_unit::Int64
    t_unit::Int64
end


"""
Constructor of GGGate from ranges in keV and calculating
Int ranges from Energy units used in file
"""
function GGGate(D::Array{Int64, 1},
             gate_z::Array{Float64, 1}, 
             gate_y::Array{Float64, 1}, 
             prompt::Array{Float64, 1},
             delayed::Array{Float64, 1},
             ttype::String,
             E_unit::Int64,
             t_unit::Int64)

    conditions = [round(Int, gate_y[1] * E_unit), 
                  round(Int, gate_y[2] * E_unit),
                  round(Int, gate_z[1] * E_unit),
                  round(Int, gate_z[2] * E_unit)]
    @assert length(ttype) == 3
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
    GGGate(D, conditions..., E_unit, t_unit)
end


"""
Constructor of GGTGate from ranges in keV and calculating
Int ranges from Energy units used in file
"""
function GGTGate(D::Array{Int64, 1},
             gate_z::Array{Float64, 1}, 
             gate_y::Array{Float64, 1}, 
             E_unit::Int64,
             t_unit::Int64)

    conditions = [round(Int, gate_y[1] * E_unit), 
                  round(Int, gate_y[2] * E_unit),
                  round(Int, gate_z[1] * E_unit),
                  round(Int, gate_z[2] * E_unit)]
    GGTGate(D, conditions..., E_unit, t_unit)
end



"""
Calculate gamma-time distribution of all events
"""
function select_t!(data::Array{UInt32, 2}, 
                                    matrix::Array{Int64, 2}, 
                                    c::TGate)
    n = size(data)[2]
    for i in 1:n
        for j in 1:3
            k = trunc(Int64, data[j, i] / c.E_unit) + 1
            t = trunc(Int64, data[j+3, i] / c.t_unit) + 1
            if 1 <= k <= c.D[1] && 1 <= t <= c.D[2]
                matrix[k, t] += 1
            end
        end
    end
end


"""
Calculate gamma-gamma-gamma gate with timing conditions (prompt/delayed) 
on all three gammas
"""
function select_ggg!(data::Array{UInt32, 2}, 
                      matrix::Array{Int64, 1}, 
                      c::GGGate)
    n = size(data)[2]
    for loc in [[1, 2, 3], [1, 3, 2], [2, 3, 1], 
                                    [2, 1, 3], [3, 1, 2], [3, 2, 1]]
        for i in 1:n
            if (   (c.y1 <= data[loc[1], i] < c.y2) 
                && (c.z1 <= data[loc[2], i] < c.z2) 
                && (c.t11 <= data[loc[1]+3, i] <= c.t12)
                && (c.t21 <= data[loc[2]+3, i] <= c.t22)
                && (c.t31 <= data[loc[3]+3, i] <= c.t32)
               )
                k = trunc(Int64, data[loc[3], i] / c.E_unit) + 1
                if k <= c.D[1]
                    matrix[k] += 1
                end
            end
        end
    end
end


"""
Calculate gamma-gamma gate returns gamma-time distribution
"""
function select_ggt!(data::Array{UInt32, 2}, 
                      matrix::Array{Int64, 2}, 
                      c::GGTGate)
    n = size(data)[2]
    for loc in [[1, 2, 3], [1, 3, 2], [2, 3, 1], 
                                    [2, 1, 3], [3, 1, 2], [3, 2, 1]]
        for i in 1:n
            if (   (c.y1 <= data[loc[1], i] < c.y2) 
                && (c.z1 <= data[loc[2], i] < c.z2) 
               )
                k = trunc(Int64, data[loc[3], i] / c.E_unit) + 1
                t = trunc(Int64, data[loc[3]+3, i] / c.t_unit) + 1
                if 1 <= k <= c.D[1] && 1 <= t <= c.D[2]
                    matrix[k, t] += 1
                end
            end
        end
    end
end


"""
This functions browses through all given muliplicity datasets
and applies 'select' function with 'c' gate conditions.

* fin: HDF5 datafile with GeGeGe group and multiplicity datasets
* gate_m: if empty list, a whole range in file is used, if length is 1, it is 
          taken as lower limit, if length is 2 - low and high end, and
          any other, as a list of multiplicities
* c: Gate type struct, must be given accordingly to the select function
* select: function performing actual 

return Array of sizes c.D
"""        
function process(fin::HDF5File, gate_m::Array{Int64, 1}, 
                c::Gate, select::Function)::Array{Int64, length(c.D)}

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
   
    matrix = Array{Int64, length(c.D)}(undef, c.D...)
    matrix = zero(matrix)
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

            select(data, matrix, c)

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

end
