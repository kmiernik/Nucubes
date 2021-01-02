import HDF5
#import Profile
#import InteractiveUtils
#import ProfileView
#import Printf
#import BenchmarkTools

include("Nucubes.jl")
import .Nucubes


function open_create_group(parent, group_name::String)
    try
        #HDF5.h5e_set_auto(HDF5.H5E_DEFAULT, C_NULL, C_NULL)
        global group = parent[group_name]
    catch 
        global group = HDF5.create_group(parent, group_name)
    end
    return group
end


function save_h5(output_file, matrix, gate_y, gate_z, gate_m, prompt, delayed, 
                 target, isotope, detectors, ttype)
    y = round(Int, (gate_y[1] + gate_y[2]) / 2)
    z = round(Int, (gate_z[1] + gate_z[2]) / 2)
    data_name = "$z" * "_" * "$y" 
    if size(gate_m)[1] > 0
        data_name = data_name * "_m"
        for mi in gate_m
            data_name = data_name * "$mi" * "_"
        end
        data_name = strip(data_name, '_')
    end

    HDF5.h5e_set_auto(HDF5.H5E_DEFAULT, C_NULL, C_NULL)
    fout = HDF5.h5open(output_file, "cw")
    target_group = open_create_group(fout, target)
    isotope_group = open_create_group(target_group, isotope)
    detectors_group = open_create_group(isotope_group, detectors)
    type_group = open_create_group(detectors_group, ttype)

    dataset_created = false
    while !dataset_created
        try
            d = HDF5.create_dataset(type_group, data_name, 
                                    HDF5.datatype(Int64), 
                                HDF5.dataspace(size(matrix)))
            dataset_created = true
            d[:] = matrix[:]
            HDF5.attributes(d)["y"] = gate_y
            HDF5.attributes(d)["z"] = gate_z
            HDF5.attributes(d)["m"] = gate_m
            HDF5.attributes(d)["prompt"] = prompt
            HDF5.attributes(d)["delayed"] = delayed
        catch
            data_name = data_name * "*"
        end
    end
    close(fout)

end


function main(gate_file::String)

    for line in eachline(gate_file)
        line = strip(line)
        if startswith(line, "#")
            continue
        end
        try
            words = split(line)
            input_file = String(words[1])
            output_file = String(words[2])
            target = String(words[3])
            isotope = String(words[4])
            detectors = String(words[5])
            ttype = String(words[6])
            gate_z = [parse(Float64, words[7]), parse(Float64, words[8])]
            gate_y = [parse(Float64, words[9]), parse(Float64, words[10])]
            gate_m = [parse(Int64, words[11]), parse(Int64, words[12])]
            prompt = [parse(Float64, words[13]), parse(Float64, words[14])]
            delayed = [parse(Float64, words[15]), parse(Float64, words[16])]

            fin = HDF5.h5open(input_file, "r")
            println("Processing gate $target $isotope $ttype z:$gate_z y:$gate_y m:$gate_m")
            matrix = Nucubes.gegege(fin, gate_z, gate_y, prompt, delayed,
                                    gate_m, ttype)
            println("Total counts ", sum(matrix))
            close(fin)

            save_h5(output_file, matrix, gate_y, gate_z, gate_m, 
                    prompt, delayed, target, isotope, detectors, ttype)
        catch err
            if isa(err, LoadError)
                continue
            end
            throw(err)
        end
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

