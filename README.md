# Agent Julia

Project code for Newcastle University MSc Bioinformatics project submission for
academic year 2020/2021.

## Table of Contents

- [Description](#description)
- [Requirements](#requirements)
- [Execution Instructions](#execution-instructions)
  - [Prerequisites](#prerequisites)
  - [Gillespie SSA Model](#gillespie-ssa-model)
  - [Agent Julia Model](#agent-julia-model)
    - [A Note on Video Rendering](#a-note-on-video-rendering)
- [Development Environment](#development-environment)
- [Resources](#resources)
- [References](#references)

## Description

This repository contains the code for the `Agent Julia` program. This is a model
for measuring the dynamics of mtDNA within a 2-dimensional space over the course
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

### Prerequisites

The models contained within these scripts require many packages. To reduce the
amount of downloads and updates required, a script named `dependencies.jl` is
included. This script must be ran before using **any** of the model files. This
script can be launched from the Julia REPL by using the following command:

```julia
include("dependencies.jl")
```

This script may take time to complete, but will reduce the amount of time
required to run a simulation, particularly if you are using many packages
already. Please note this script may also perform updates on installed packages.

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
`Gillespie SSA` model, and work in the same way. This model will export three
types of result by default:

1. A CSV file in the format of `patient_type.csv` within the `data` directory.
2. Two plots in the format of `0X plot_type.png` within the `plots/patient_type`
   directory.
3. A video file in the format of `patient_type.mp4` within the `videos`
   directory.

This approach allows the user to keep all output in a consistent manner and all
examples of these data and outputs may be found within the project's `examples`
directory. These examples are updated along with the code regularly to ensure
that they are current and relevant.

#### A Note on Video Rendering

By default, Agent Julia simulation files will generate a video simulation using
the given parameters will be created. These simulation videos are not
representative of the values exported in CSV format or shown in the plots. They
are separate, but they do make use of the same parameters as the current
simulation and also the set random seed - this means that although they are a
different simulation, they will be very similar to the results shown within the
plots and CSV file.

Note that these scripts have a companion called `integrated_gpu_support.sh`,
which is intended to be ran in the event that the user has a machine with an
integrated GPU. This is due to a change within the `CairoMakie` package, this
package is used to render the output. Users of **NVIDIA** or **AMD** graphics
chipsets need not use this BASH script.

It is possible to execute the companion script using the Terminal with the
following command (This script may also be executed using the Julia REPL by
switching into shell mode by using the `;` key and using the same commands):

```bash
bash "integrated_gpu_support.sh" "name_of_use_case.jl"

# OR

./"integrated_gpu_support.sh" "name_of_use_case.jl"
```

This script uses `sudo` privileges to execute the model and thus create the
video output. This means that a user may need to manually change the permissions
of the video file once it has been created, this is possible using the following
command as your current user (not `root`):

```bash
sudo chown "$USER":"$USER" "name_of_file.mp4"
```

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
