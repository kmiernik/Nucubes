"""
    convert(listfile, tomlfile)

    Convert text file with list of gates to toml config file 
    required by nucubes.jl

    Structure of input file is
    name z_gate y_gate

    where name is isotope/gate name and z_gate y_gate are floating
    point values
"""
function convert(listfile, tomlfile)
    fin = open(listfile, "r")

    gates = Dict{String, Any}()
    for line in readlines(fin)
        if startswith(line, "#")
            continue
        end
        data = split(line)
        isotope = data[1]
        z = parse(Float64, data[2])
        y = parse(Float64, data[3])
        name = isotope
        sub = "abcdefghijklmnoprstuwxyz"
        for i in 1:length(sub)
            if haskey(gates, name)
                name = isotope * sub[i]
            else
                break
            end
        end
        gates[name] = [z, y]
        gates[name * "bg"] = [z, y + 10]
    end

    labels = sort(collect(keys(gates)))

    fout = open(tomlfile, "w")

default = """
[default]
    input = "ggg_ppp_all.h5"
    output = "even_even.h5"
    target = "au"
    detectors = "ggg"
    ttype = "ppp"
    halfwidth = 3.0
    m = [3, 9]
    prompt = [35.0, 65.0]
    delayed = [80.0, 350.0]
"""
    println(fout, default)
    
    for entry in labels
        print(fout, "[gate.$entry]\n")
        print(fout, "    isotope = \"", entry, "\"\n")
        print(fout, "    z = ", gates[entry][1], "\n")
        print(fout, "    y = ", gates[entry][2], "\n\n")
    end

    close(fout)

    #=
    parnames = ["dE", "Emax", "dt_ge", "dt_la", "tmax", "Mmax", "chunk_size",
                "coin", "bgo_low", "ge_low", "t_prompt_low", "t_prompt_high"]
    for par in parnames
        print(fout, "\t$par = ", config["spectra"][par], "\n")
    end
    print(fout, "\n")

    for entry in labels
        print(fout, "[label.", entry.first, "]\n")
        print(fout, "\tname = \"", entry.second["name"], "\"\n")
        print(fout, "\tpos = \"", entry.second["pos"], "\"\n")
        print(fout, "\ttype = \"", entry.second["type"], "\"\n")
        print(fout, "\tdt = ", round(entry.second["dt"], digits=3), "\n")
        if haskey(entry.second, "cal")
            print(fout, "\tcal = ", entry.second["cal"], "\n")
        else
            print(fout, "\tcal = [0.0, 1.0, 0.0]\n")
        end
        if haskey(entry.second, "walk")
            print(fout, "\twalk = ", entry.second["walk"], "\n")
        else
            print(fout, "\twalk = [0.0, 1.0, 0.0]\n")
        end
        print(fout, "\tvalid = ", entry.second["valid"], "\n")
        if haskey(entry.second, "comment")
            print(fout, "\tcomment = \"", entry.second["comment"], "\"\n")
        else
            print(fout, "\tcomment = \" \"\n")
        end
        print(fout, "\n")
    end
    =#
    close(fout)
end
