import HDF5
import Dates
using Distributed
using TOML

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
        if startswith(line, "#") || length(line) == 0
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
                gate["ggg"] = Nucubes.GG_G([8192], gate["gate_z"],
                                            gate["gate_y"], gate["prompt"],
                                            gate["delayed"], gate["ttype"],
                                            100, 1000)
                gate["select"] = Nucubes.select_gg_g!
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


function parse_toml(tomlfile)
    config = TOML.parsefile(tomlfile)
    gates = Array{Dict{String, Any}, 1}(undef, 0)

    conv = Dict("input" => "input_file", "output" => "output_file",
                "target" => "target", "detectors" => "detectors",
                "ttype" => "ttype", "m" => "gate_m", 
                "prompt" => "prompt", "delayed" => "delayed")

    default = Dict{String, Any}()
    for c in conv
        default[c.second] = config["default"][c.first]
    end

    for entry in config["gate"]
        if haskey(entry.second, "valid")
            if !entry.second["valid"]
                continue
            end
        end
        gate = copy(default)
        gate["isotope"] = entry.second["isotope"]
        
        for x in ["z", "y"]
            if isa(entry.second[x], Vector)
                gate["gate_" * x] = entry.second[x]
            elseif isa(entry.second[x], Number)
                hw = config["default"]["halfwidth"]
                if haskey(entry.second, "halfwidth")
                    hw = entry.second["halfwidth"]
                end
                gate["gate_" * x] = [entry.second[x] - hw, 
                                     entry.second[x] + hw]
            end
        end

        for par in keys(entry.second)
            if par in keys(conv)
                gate[conv[par]] = entry.second[par]
            end
        end

        gate["ggg"] = Nucubes.GG_G([8192], gate["gate_z"],
                                    gate["gate_y"], gate["prompt"],
                                    gate["delayed"], gate["ttype"],
                                    100, 1000)
        gate["select"] = Nucubes.select_gg_g!
        println("* gate ", gate["target"], " ", gate["isotope"], " ",
                gate["ttype"], " z:", gate["gate_z"], " y:", gate["gate_y"],
                " m:", gate["gate_m"], " added")
        push!(gates, gate)
    end
    gates
end


function main(gate_file)
    
    gates = parse_toml(gate_file)
    t0 = Dates.Time(Dates.now())
    results = pmap(gate->Nucubes.process(gate["input_file"], gate["gate_m"],
                                         gate["ggg"], gate["select"]), 
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
    println("USAGE: nucubes.jl gate_file.toml")
end
