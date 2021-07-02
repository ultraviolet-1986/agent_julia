# Agent Julia

Project code for Newcastle University MSc Bioinformatics project submission for
academic year 2020/2021.

## Table of Contents

- [Description](#description)
- [Requirements](#requirements)
- [Execution Instructions](#execution-instructions)
  - [Agent Julia Model](#agent-julia-model)
  - [Gillespie SSA Model](#gillespie-ssa-model)
- [Development Environment](#development-environment)
- [Resources](#resources)
- [References](#references)

## Description

This repository contains the code for the `Agent Julia` program. This is a model
for measuring the dynamics of mtDNA within a 3-dimensional space over the course
of a human life-span of approximately 90 years. It is written in the Julia
programming language (v1.6.0) and makes use of the `Agents.jl` package.

The timeline of this project is: `26/04/2021` to `27/08/2021` and package
versions will be contemporary to this time and no changes will be applied once
the code is ready for submission. Updates past this time represent a
continuation of the project outside of an academic context.

## Requirements

At the time of writing, this software requires `Julia` versions greater than or
equal to `v1.6.0` and previous versions have not been accounted for.

## Execution Instructions

This program was written to be executed from the command line but may also be
executed from within the Julia REPL.

### Agent Julia Model

Note that this script has a companion called `integrated_gpu_support.sh`, which
is intended to be ran in the event that the user has a machine with an
integrated GPU. This is due to a change within the `CairoMakie` package, this
package is used to render the output.

It is possible to execute the companion script using the Terminal with the
following command (This script may also be executed using the Julia REPL by
switching into shell mode by using the `;` key and using the same commands):

```bash
bash integrated_gpu_support.sh

# OR

./integrated_gpu_support.sh
```

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

## References

This section will contain a list of academic papers which were used to create
this software.
