# Nucubes
This package provides tools for nuBall HDF5 data files gating and
analysis. The main script is nucubes.jl which calculates gamma-gamma-gamma
gates based on a configuration file (see example file gates_example.toml).
This code is prepared for multi-threading and distributed computing.

## Setup
1. Copy or clone the git repository to your local machine to some `my_path` 
  (e.g. /home/user/codes/Julia/Nucubes)

2. Run julia, and in the REPL
```julia-repl
julia> using Pkg
julia> Pkg.add(url="https://github.com/kmiernik/Nucubes.git")
```

3. For distributed computing:
    * Prepare passwordless login on all your machines (e.g. using keys, 
      see [this](https://wiki.archlinux.org/index.php/SSH_keys) for
      instructions)
    * Perform step 2 on all your machines
    * Place the data files in the same location (e.g. /home/user/data/data.h5, or network drive) on all machines
    * For convenience create symbolic link to nucubes.jl script
    ```
    $ cd /home/user/data/
    $ sudo ln -s my_path/bin/nucubes.jl .
    ```
    * On your main computer create machines.txt file listing all your 
    machines by IP or name, e.g.
    ```
    server.somewhere.edu.pl
    10.13.2.1
    ```
    * Prepare gates file based on `gates_example.txt`, with the above 
      data file path

## Running
    
* For local run
    ```
    $ julia [-t n] [-p m] nucubes.jl gates.txt
    ```
    where `-t n` setups `n` cores to be used and `-p m` to use `m` processes

* For distributed run
    ```
    $ julia [-t n] [-p m] --machine-file machines.txt nucubes.jl gates.txt
    ```
    where `-t n` setups number of cores to be used per machine, machine
    file defines where the calculations are to be run and `gates.txt`
    contain gates data. If you want to make calculations on your local
    computer as well add `-p m`.

## Setup file (gates)
Gates are created in a TOML setup file. An example with explanations is given below

```
# Set default values for these parameters
[default]
    input = "data.h5"
    output = "gates.h5"
    target = "target_name"
    detectors = "ggg"
    ttype = "ppp"
    halfwidth = 1.0
    m = [3, 9]
    prompt = [20.0, 60.0]
    delayed = [80.0 400.0]

# Each gate always needs "z", "y" and "isotope"
# Each gate has some unique identifier (anything)
# comment is optional and is skipped during parsing
# valid = false can be used to skip the gate (valid = true is default,
and does not need to be placed)
[gate.1]
    z = 606.6
    y = 423.4
    isotope = "te136"
    comment = "4+->2+->0+"
    valid = false

# Default parameters can be overwritten for a given gate
[gate.2]
    input = "data2.h5"
    target = "some_other_target"
    z = 423.4
    y = 352.3
    isotope = "te136"
    m = [4, 9]

# If z(y) is number, the gate is +/- halfwidth, 
# if it is a vector, it will be used as lower and upper gate limit
# (so it can be non-symmetric if needed)
[gate.3]
    z = [422.4, 425.4]
    y = 352.3
    isotope = "te136"
    ttype = "ppd"
    comment = "4+->2+->0+"
```
