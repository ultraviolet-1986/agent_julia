# Agent Julia

Project code for Newcastle University MSc Bioinformatics project submission for
academic year 2020/2021.

## Table of Contents

- [Description](#description)
- [Requirements](#requirements)
- [Execution Instructions](#execution-instructions)
  - [Gillespie SSA Model](#gillespie-ssa-model)
  - [Agent Julia Model](#agent-julia-model)
    - [A Note on Video Rendering](#a-note-on-video-rendering)
- [Development Environment](#development-environment)
- [Resources](#resources)
- [References](#references)

## Description

This repository contains the code for the `Agent Julia` program. This is a model
for measuring the dynamics of mtDNA within a 3-dimensional space over the course
of a human life-span of approximately 90 years. It is written in the Julia
programming language (v1.6.x) and makes use of the `Agents.jl` package.

The timeline of this project is: `26/04/2021` to `27/08/2021` and package
versions will be contemporary to this time and no changes will be applied once
the code is ready for submission. Updates past this time represent a
continuation of the project outside of an academic context.

## Requirements

At the time of writing, this software requires `Julia` versions greater than or
equal to `v1.6.0` and previous versions have not been accounted for.

Each script contained within this repository will install any dependencies they
require before execution so no extra software is required except a `Julia REPL`
and a `BASH` terminal for video rendering on integrated graphics scenarios (see
below).

## Execution Instructions

This program was written to be executed from the command line but may also be
executed from within the Julia REPL.

### Gillespie SSA Model

The code for the Gillespie stochastic simulation algorithm is stored within the
`gillespie_model.jl` file within the `Gillespie SSA Model` directory. It
contains the functions required to perform a stochastic simulation given
user-defined parameters which may be written to another `Julia` script. These
scripts can contain parameters for simulating different situations and lengths
of time.

It is possible to execute a use-case script from the Terminal (BASH or
compatible) by using the following command from the appropriate working
directory:

```bash
julia name_of_use_case.jl

# OR

./name_of_use_case.jl
```

It is also possible to execute these scripts directly from within the Julia REPL
by using the following command:

```julia
include("name_of_use_case.jl")
```

In both events, the use-case script will define the simulation parameters and
call the model script to perform the simulation, finally providing a report and
a plot showing the results of the simulation.

Note that *Microsoft Visual Studio Code* users may simply open the simulation
file of their choice and run the file using the built-in debugger, or the `Run`
button toward the right of the tab bar.

### Agent Julia Model

This model also contains use-case simulation files in the same manner as the
`Gillespie SSA` model, these work in the same way. Plots will be placed within
the same directory as the `agent_julia_model.jl` script. This model will export
two types of result by default: **1)** a CSV file containing agent counts at
each given time point (duration of 1 month by default), and **2)** a graphical
plot of these data. It is also possible to render a video of a simulation (see
below).

#### A Note on Video Rendering

In the `agent_julia_model.jl` file, it is possible to generate a video of the
simulation by using the command `simulation_to_video()` function within the
Julia REPL once the simulation has completed. Because of the use of a set seed,
results should be the same as results shown within the plot, but this cannot be
verified and there is currently no way to export the simulation data directly to
video and so the same simulation will be executed again but results are written
to a video file as the simulation occurs.

Note that these scripts have a companion called `integrated_gpu_support.sh`,
which is intended to be ran in the event that the user has a machine with an
integrated GPU. This is due to a change within the `CairoMakie` package, this
package is used to render the output. Users of NVIDIA or AMD graphics chipsets
need not use this BASH script.

It is possible to execute the companion script using the Terminal with the
following command (This script may also be executed using the Julia REPL by
switching into shell mode by using the `;` key and using the same commands):

```bash
bash integrated_gpu_support.sh

# OR

./integrated_gpu_support.sh
```

This script uses `sudo` privileges to execute the model and thus create the
video output. Once the simulation has completed, the script will use the same
permissions to change the owner of the video file to the current user, otherwise
the owner of the file will be `root`. In a `Toolbox` environment, you will not
be asked for your password and the results should be the same.

## Development Environment

The following software and tools form the development environment (at
time of writing):

- [Fedora Silverblue Linux 34](https://silverblue.fedoraproject.org/)
- [Toolbox](https://github.com/containers/toolbox)
- [Microsoft Visual Studio Code](https://code.visualstudio.com/)
- [Julia Programming Language v1.6.x](https://julialang.org/)
- [Agents.jl](https://juliadynamics.github.io/Agents.jl/stable/)

## Resources

This section will contain a list of resources including code and sources of data
which were used to create this software.

- [Gillespie Model in Julia](https://nextjournal.com/bebi5009/gillespie-julia)
- [Continuous space social distancing for COVID-19](https://git.io/Jc1w6)

## References

Datseris, G., Vahdati, A.R., and DuBois, T.C. (2021) 'Agents.jl: A performant
and feature-full agent based modelling software of minimal code complexity'.
[Preprint]. Available at: <https://arxiv.org/abs/2101.10072>
(Accessed: 26/04/2021).
