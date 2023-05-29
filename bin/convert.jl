using CSV


"""
    list2toml(listfile, tomlfile)

    Convert text file with list of gates to toml config file 
    required by nucubes.jl

    Structure of input file is
    name z_gate y_gate

    where name is isotope/gate name and z_gate y_gate are floating
    point values
"""
function list2toml(listfile, tomlfile, bg=false)
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
        sub = 0
        while true
            if haskey(gates, name)
                sub += 1
                name = isotope * "_$sub"
            else
                break
            end
        end
        gates[name] = [z, y]
        if bg
            gates[name * "bg"] = [z, y + 10]
        end
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
end


"""
    fits2toml(fitsfile, tomlfile)

    Convert CSV file with fits to gates as required by nucubes.jl

    This one is for the beta-decay background check in triples

    Structure of fits file is
    isotope,Z,A,E1,I1,J1,E2,I2,J2,(E3,I3,J3,Af,Ef,sf,dAf,KS)

    E1 and E2 are actually needed for gate to be created
"""
function fits2toml(fitsfile, tomlfile)
    fits = CSV.File(fitsfile, comment="#")

    gates = Dict{String, Any}()
    for row in fits
        z = row.E1
        y = row.E2
        isotope = row.isotope
        name = isotope
        sub = 0

        already = false
        while true
            if haskey(gates, name)
                if gates[name] == [z, y]
                    already = true
                    break
                end
                sub += 1
                name = isotope * "_$sub"
            else
                break
            end
        end
        if already
            continue
        end
        gates[name] = [z, y]
    end

    labels = sort(collect(keys(gates)))

    fout = open(tomlfile, "w")

default = """
[default]
    input = "ggg_all.h5"
    output = "beta_bg_fits.h5"
    target = "au"
    detectors = "ggg"
    ttype = "ddd"
    halfwidth = 3.0
    m = [3, 9]
    prompt = [35.0, 65.0]
    delayed = [370.0, 400.0]
"""
    println(fout, default)
    
    for entry in labels
        print(fout, "[gate.$entry]\n")
        print(fout, "    isotope = \"", entry, "\"\n")
        print(fout, "    z = ", gates[entry][1], "\n")
        print(fout, "    y = ", gates[entry][2], "\n\n")
    end

    close(fout)
end
