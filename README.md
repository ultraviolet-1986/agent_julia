# Agent Julia

Project code for Newcastle University MSc Bioinformatics project
submission for academic year 2020/2021.

## Table of Contents

- [Description](#description)
- [Execution Instructions](#execution-instructions)
- [Development Environment](#development-environment)
- [Resources](#resources)
- [References](#references)

## Description

This repository contains the code for the `Agent Julia` program. This is
a model for measuring the dynamics of mtDNA within a 3-dimensional space
over the course of a human life-span of approximately 90 years. It is
written in the Julia programming language (v1.6.0) and makes use of the
`Agents.jl` package.

The timeline of this project is: `26/04/2021` to `20/08/2021` and
package versions will be contemporary to this time and no changes will
be applied once the code is ready for submission. Updates past this
time represent a continuation of the project outside of an academic
context.

## Execution Instructions

This program was written to be executed from the command line. On the
command-line, use the following command(s):

```bash
# Execute Agent-based simulation.
julia agent_julia.jl

# Execute stochastic Gillespie simulation.
julia gillespie_model.jl
```

Note that either of these scripts may be executed from the Julia REPL by using
the following command(s):

```julia
# Execute Agent-based simulation.
include("agent_julia.jl")

# Execute stochastic Gillespie simulation.
include("gillespie_model.jl")
```

## Development Environment

The following software and tools form the development environment (at
time of writing):

- [Fedora Silverblue Linux 34](https://silverblue.fedoraproject.org/)
- [Toolbox](https://github.com/containers/toolbox)
- [Microsoft Visual Studio Code](https://code.visualstudio.com/)
- [Julia Programming Language](https://julialang.org/)
- [Agents.jl](https://juliadynamics.github.io/Agents.jl/stable/)

## Resources

This section will contain a list of resources including code and
sources of data which were used to create this software.

## References

This section will contain a list of academic papers which were used to
create this software.
