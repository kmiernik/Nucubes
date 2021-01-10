import HDF5
import Dates
using Distributed

@everywhere import Nucubes

function open_create_group(parent, group_name::String)
    try
        #HDF5.h5e_set_auto(HDF5.H5E_DEFAULT, C_NULL, C_NULL)
        global group = parent[group_name]
    catch 
        global group = HDF5.create_group(parent, group_name)
    end
    return group
end


function save_h5(g::Dict{String, Any}, matrix::AbstractArray)

    y = round(Int, (g["gate_y"][1] + g["gate_y"][2]) / 2)
    z = round(Int, (g["gate_z"][1] + g["gate_z"][2]) / 2)
    data_name = "$z" * "_" * "$y" 
    if size(g["gate_m"])[1] > 0
        data_name = data_name * "_m"
        for i in 1:length(g["gate_m"])-1
                            mi = g["gate_m"][i]
            data_name = data_name * "$mi" * "_"
        end
        mi = last(g["gate_m"])
        data_name = data_name * "$mi"
    end

    # Disables error messages caugth now silently by the try blocks
    HDF5.h5e_set_auto(HDF5.H5E_DEFAULT, C_NULL, C_NULL)

    HDF5.h5open(g["output_file"], "cw") do fout
        target_group = open_create_group(fout, g["target"])
        isotope_group = open_create_group(target_group, g["isotope"])
        detectors_group = open_create_group(isotope_group, g["detectors"])
        type_group = open_create_group(detectors_group, g["ttype"])

        dataset_created = false
        while !dataset_created
            try
                d = HDF5.create_dataset(type_group, data_name, 
                                        HDF5.datatype(Int64), 
                                    HDF5.dataspace(size(matrix)))
                dataset_created = true
                d[:] = matrix[:]
                HDF5.attributes(d)["y"] = g["gate_y"]
                HDF5.attributes(d)["z"] = g["gate_z"]
                HDF5.attributes(d)["m"] = g["gate_m"]
                HDF5.attributes(d)["prompt"] = g["prompt"]
                HDF5.attributes(d)["delayed"] = g["delayed"]
            catch
                data_name = data_name * "*"
            end
        end
    end

end

function parse_input(gate_file::String)
    
    gates = []
    for line in eachline(gate_file)
        line = strip(line)
        if startswith(line, "#")
            continue
        end
        gate = Dict{String, Any}()
        try
            words = split(line)
            gate["input_file"] = String(words[1])
            gate["output_file"] = String(words[2])
            gate["target"] = String(words[3])
            gate["isotope"] = String(words[4])
            gate["detectors"] = String(words[5])
            gate["ttype"] = String(words[6])
            gate["gate_z"] = [parse(Float64, words[7]), 
                              parse(Float64, words[8])]
            gate["gate_y"] = [parse(Float64, words[9]), 
                              parse(Float64, words[10])]
            gate["gate_m"] = [parse(Int64, words[11]), parse(Int64, words[12])]
            gate["prompt"] = [parse(Float64, words[13]), 
                              parse(Float64, words[14])]
            gate["delayed"] = [parse(Float64, words[15]), 
                            parse(Float64, words[16])]

            if gate["detectors"] == "ggg"
                gate["ggg"] = Nucubes.GGGate([4096], gate["gate_z"],
                                            gate["gate_y"], gate["prompt"],
                                            gate["delayed"], gate["ttype"],
                                            100, 1000)
                gate["select"] = Nucubes.select_ggg!
            else
                println("Gate on detectors: \"", gate["detectors"], 
                        "\" not implemented")
                continue
            end
            println("* gate ", gate["target"], " ", gate["isotope"], " ",
                    gate["ttype"], " z:", gate["gate_z"], " y:", gate["gate_y"],
                    " m:", gate["gate_m"], " added")
            push!(gates, gate)
        catch err
            if isa(err, LoadError)
                continue
            end
            throw(err)
        end
    end
    return gates
end


function main(gate_file::String)
    
    gates = parse_input(gate_file)
    t0 = Dates.Time(Dates.now())
    results = pmap(gate->Nucubes.process(gate["input_file"], gate["gate_m"],
                                        gate["ggg"], Nucubes.select_ggg!), 
                   gates)
    t1 = Dates.Time(Dates.now())
    dt = t1 - t0
    println("* ", length(gates), " gates processed in ", 
            round(dt.value * 1e-9, digits=1), " s ")
    for i in eachindex(results)
        save_h5(gates[i], results[i])
        println("-> gate ", gates[i]["target"], " ", gates[i]["isotope"], " ",
                gates[i]["ttype"], " z:", gates[i]["gate_z"], " y:",
                gates[i]["gate_y"], " m:", gates[i]["gate_m"], 
                ", counts ", sum(results[i]))
    end
end


if length(ARGS) >= 1
    main(ARGS[1])
else
    println("USAGE: nucubes.jl gate_file.txt")
    println()
    println("Where gate_file contain the following space separated columns")
    println("   (1) input HDF5 file name (string) ")
    println("   (2) output HDF5 file name (string) ")
    println("   (3) target alias (string) ")
    println("   (4) isotope of interest  (string) ")
    println("   (5) detectors (ggg, ggl, gll)  (string) ")
    println("   (6) timing gate type (aaa, ppp, ppd, pdd, ddd) (string) ")
    println("   (7) Z-axis gate begin (>=) (float) ")
    println("   (8) Z-axis gate end (<) (float) ")
    println("   (9) Y-axis gate begin (>=) (float) ")
    println("   (10) Y-axis gate end (<) (float) ")
    println("   (11) multiplicity gate begin (>=) (int) ")
    println("   (12) multiplicity gate end (<=) (int) ")
    println("   (13) prompt gate begin (>=) (float) ")
    println("   (14) prompt gate end (<=) (float) ")
    println("   (15) delayed gate begin (>=) (float) ")
    println("   (16) delayed gate end (<=) (float) ")
    println("  Comment lines are indicated by '#' sign at the beginning")
end

# Types warnings
#@InteractiveUtils.code_warntype XXX

# Profiling
#@profile XXX
#Profile.print(format=:flat)
#ProfileView.view()
#println("Press enter")
#s = readline()

