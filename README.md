# Nucubes
This package provides tools for nuBall HDF5 data files gating and
analysis. The main script is nucubes.jl which calculates gamma-gamma-gamma
gates based on a configuration file (see example file gates_example.txt).
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
    * Place the data files in the same location (e.g. /home/user/data/data.h5)
      on all machines
    * For convencince create symbolic link to nucubes.jl script
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
    * Prapare gates file based on `gates_example.txt`, with the above 
      data file path

## Running
    
* For local run
    ```
    $ julia [-t n] [-p m] nucubes.jl gates.txt
    ```
    where `-t n` setups `n` cores to be used and `-p m` to use `m` processes

* For distrubuted run
    ```
    $ julia [-t n] [-p m] --machine-file machines.txt nucubes.jl gates.txt
    ```
    where `-t n` setups number of cores to be used per machine, machine
    file defines where the calculations are to be run and `gates.txt`
    contain gates data. If you want to make calculations on your local
    computer as well add `-p m`.
