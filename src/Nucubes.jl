"""
Module for NuBall cubic gates analysis
Version with conditions and process functions
"""
module Nucubes

using HDF5
using Dates
using Base.Threads
using Distributed

export process
export Gate
export GGGate
export GGTate
export TGate
export select_ggg!
export select_ggt!
export select_t!

"""
This is abstract Gate type that is passed to the process function
The concrete type should define an array
    D::Array{Int64, 1}
with return array sizes (e.g. [4096] will result in 1D array of 4096 channels, 
and [4096, 4096] will result in 2D array 4096*4096).
And any other needed conditions. The Gate will be passed to the select
function which will populate the return array and use conditions.
"""
abstract type Gate
end

"""
Gate for gamma-gamma-gamma selection, the D should be single entry
(e.g. [4096] or [2048]) which defines a histogram 1 keV/bin.
Other fields are:
    * y1 and y2 - gate range in the same units as in the file y1 <= x < y
    * z1 and z2 - as above for the second gamma
    * t11 and t12 - time range for the first gamma in units as in the data file
                    t11 <= t1 <= t12
    * t21 and t22 - as above for the second gamma
    * t31 and t32 - as above for the third gamma
    * E_unit - unit value to convert from data file to keV (e.g. 100, 
                than value E(keV) = E / 100)
    * t_unit - same for the time 
"""
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

"""
Gate for creating gamma-time distribution for a given gamma-gamma condition
see GGGate for more details on fields
"""
struct GGTGate <: Gate
    D::Array{Int64, 1}
    y1::Int64
    y2::Int64
    z1::Int64
    z2::Int64
    E_unit::Int64
    t_unit::Int64
end


"""
Gate for creating total gamma-time distribution for all gammas
see GGGate for more details on fields
"""
struct TGate <: Gate
    D::Array{Int64, 1}
    E_unit::Int64
    t_unit::Int64
end


"""
    GGGate(D::Array{Int64, 1},
           gate_z::Array{Float64, 1}, 
           gate_y::Array{Float64, 1}, 
           prompt::Array{Float64, 1},
           delayed::Array{Float64, 1},
           ttype::String,
           E_unit::Int64,
           t_unit::Int64)
Constructor of GGGate from ranges in keV and calculating
Int ranges for GGGate from energy and time units as used in the data file
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
    GGTGate(D::Array{Int64, 1},
            gate_z::Array{Float64, 1}, 
            gate_y::Array{Float64, 1}, 
            E_unit::Int64,
            t_unit::Int64)

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
    select_t!(data::Array{UInt32, 2}, 
                   matrix::Array{Int64, 2}, 
                   c::TGate, 
                   matrix_lock::ReentrantLock)

Calculate gamma-time distribution of all events, populates matrix (histogram)
"""
function select_t!(data::Array{UInt32, 2}, 
                   matrix::Array{Int64, 2}, 
                   c::TGate, 
                   matrix_lock::ReentrantLock)
    n = size(data)[2]
    for i in 1:n
        for j in 1:3
            k = trunc(Int64, data[j, i] / c.E_unit) + 1
            t = trunc(Int64, data[j+3, i] / c.t_unit) + 1
            if k <= c.D[1] && t <= c.D[2]
                lock(matrix_lock) do
                    matrix[k, t] += k
                end
            end
        end
    end
end


"""
    select_ggg!(data::Array{UInt32, 2}, 
                matrix::Array{Int64, 1}, 
                c::GGGate, 
                matrix_lock::ReentrantLock)

Calculate gamma-gamma-gamma gate with timing conditions (prompt/delayed) 
on all three gammas. Populates matrix (histogram).
"""
function select_ggg!(data::Array{UInt32, 2}, 
                     matrix::Array{Int64, 1}, 
                     c::GGGate, 
                     matrix_lock::ReentrantLock)
    ml = [zeros(size(matrix)) for i in 1:nthreads()]
    n = size(data)[2]
    @inbounds for loc in [[1, 2, 3], [1, 3, 2], [2, 3, 1], 
                [2, 1, 3], [3, 1, 2], [3, 2, 1]]
        @threads for i in 1:n
            if (   (c.y1 <= data[loc[1], i] < c.y2) 
                && (c.z1 <= data[loc[2], i] < c.z2) 
                && (c.t11 <= data[loc[1]+3, i] <= c.t12)
                && (c.t21 <= data[loc[2]+3, i] <= c.t22)
                && (c.t31 <= data[loc[3]+3, i] <= c.t32)
               )
                k = trunc(Int64, data[loc[3], i] / c.E_unit) + 1
                if (k <= c.D[1])
                    ml[threadid()][k] += 1
                end
            end
        end
    end
    lock(matrix_lock) do
        for i in 1:nthreads()
            matrix .+= ml[i]
        end
    end
end


"""
    select_ggt!(data::Array{UInt32, 2}, 
                      matrix::Array{Int64, 2}, 
                      c::GGTGate,
                      matrix_lock::ReentrantLock)

Based on gamma-gamma gate calculate gamma-time distribution.
Populates matrix (histogram).
"""
function select_ggt!(data::Array{UInt32, 2}, 
                      matrix::Array{Int64, 2}, 
                      c::GGTGate,
                      matrix_lock::ReentrantLock)
    n = size(data)[2]
    for loc in [[1, 2, 3], [1, 3, 2], [2, 3, 1], 
                                    [2, 1, 3], [3, 1, 2], [3, 2, 1]]
        for i in 1:n
            if (   (c.y1 <= data[loc[1], i] < c.y2) 
                && (c.z1 <= data[loc[2], i] < c.z2) 
               )
                k = trunc(Int64, data[loc[3], i] / c.E_unit) + 1
                t = trunc(Int64, data[loc[3]+3, i] / c.t_unit) + 1
                if k <= c.D[1] && t <= c.D[2]
                    lock(matrix_lock) do
                        matrix[k, t] += 1
                    end
                end
            end
        end
    end
end


"""
    process(input_file::String, gate_m::Array{Int64, 1}, 
            c::Gate, select::Function)::Array{Int64, length(c.D)}

This functions browses through all given muliplicity datasets
and applies 'select' function with 'c' gate conditions.

* fin: HDF5 datafile with GeGeGe group and multiplicity datasets
* gate_m: if empty list, a whole range in file is used, if length is 1, it is 
          taken as lower limit, if length is 2 - low and high end, and
          any other, as a list of multiplicities
* c: Gate type struct, must be given accordingly to the select function
* select: function performing actual calculations
* returns Array as defined by c.D (e.g. if c.D = [4096] the return array
   is 1D of 4096 channels)

This function is multithreaded, so data loading and processing use separate
threads (if possible), as well as process function may be multithreade
for a given data subset.

"""        
function process(input_file::String, gate_m::Array{Int64, 1}, 
                c::Gate, select::Function)::Array{Int64, length(c.D)}

    matrix = Array{Int64, length(c.D)}(undef, c.D...)
    matrix = zero(matrix)

    HDF5.h5open(input_file, "r") do fin
        group = fin["GeGeGe"]

        M = Int64[]
        for m_set in keys(group)
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
    
        matrix_lock = Base.ReentrantLock()
        n_all = 0
        n_processed = 0
        for m in multi
            dataset = group[string(m)]
            n_all += size(dataset)[2]
        end

        t0 = Dates.Time(Dates.now())
        
        #i_report = 0
        #n_report = 10_000_000
        workers = Array{Task, 1}()
        for m in multi
            dataset = group[string(m)]
            n = size(dataset)[2]
            chunk_size = min(n, 10_000)
            try
                chunk_size = HDF5.get_chunk(dataset)[2]
            catch
            end
            left_pos = 1
            is_something_left = true
            while is_something_left
                filter!(w->!istaskdone(w), workers)
                right_pos = left_pos + chunk_size - 1
                if right_pos > n
                    right_pos = n
                    is_something_left = false
                elseif left_pos >= right_pos
                    break
                end
                data = dataset[:, left_pos:right_pos]

                push!(workers, @Threads.spawn select(data, matrix, c,
                                                     matrix_lock))

                n_processed += right_pos - left_pos

                #j_report = trunc(Int64, n_processed / n_report)
                #if j_report > i_report
                #    i_report += 1
                #    t1 = Dates.Time(Dates.now())
                #    dt = t1 - t0
                #    print("\r", round(n_processed / n_all * 100, digits=1), 
                #          "% ", round(dt.value * 1e-9, digits=1), " s (",
                #        round(dt.value * 1e-9 * n_all / n_processed, digits=0),
                #        " s)          ")
                #end

                left_pos = right_pos + 1
            end
        end
        for w in workers
            wait(w)
        end
        t1 = Dates.Time(Dates.now())
        dt = t1 - t0
        #print("\r+ gate")
        print("+ gate")
        try
            print(" z:[", round(c.z1 / c.E_unit, digits=1), ", ",
                          round(c.z2 / c.E_unit, digits=1), "]",
                  " y:[", round(c.y1 / c.E_unit, digits=1), ", ",
                          round(c.y2 / c.E_unit, digits=1), "]"
                 )
        catch
        end
        println(" m:", gate_m, " ", round(n_all / 1e6, digits=2), 
                "M hits in ", round(dt.value * 1e-9, digits=1), " s @", myid())
    end
    return matrix

end

end
