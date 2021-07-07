#!/usr/bin/env julia

###########
# License #
###########

# Agent Julia
# Project code for Newcastle University MSc Bioinformatics project
# submission.
# Copyright (C) 2021 William Willis Whinn

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

##############
# References #
##############

# - [Agents.jl Tutorial](https://git.io/Jc1w6)

###########
# Imports #
###########

import Pkg

Pkg.add("Agents")
Pkg.add("CairoMakie")
Pkg.add("InteractiveDynamics")
Pkg.add("Random")

using Agents
using CairoMakie
using InteractiveDynamics
using Random

#################
# Prerequisites #
#################

# Ensure results are reproducible.
Random.seed!(41269)

#############
# Functions #
#############

mutable struct Wild_mtDNA <: AbstractAgent
    id::Int
    pos::NTuple{2,Float64}
    vel::NTuple{2,Float64}
    mass::Float64
end


# Attempt to inherit the 'Wild' struct.
mutable struct Mutant_mtDNA <: AbstractAgent
    id::Int
    pos::NTuple{2,Float64}
    vel::NTuple{2,Float64}
    mass::Float64
    days_mutated::Int  # number of days since mutated, may be useless
    status::Symbol     # :S, :I or :R, may be useless
    β::Float64         # May also be useless
end


function ball_model(; speed = 0.002)
    space2d = ContinuousSpace((1, 1), 0.02)

    model = ABM(
        Wild_mtDNA,
        space2d,
        properties = Dict(:dt => 1.0),
        rng = MersenneTwister(41269)
    )

    # model2 = ABM(
    #     Mutant_mtDNA,
    #     space2d,
    #     properties = Dict(:dt => 1.0),
    #     rng = MersenneTwister(41269)
    # )

    # And add some agents to the model.
    # NOTE Change this to include the wild-type and mutant mtDNA.
    for ind in 1:200
        pos = Tuple(rand(model.rng, 2))
        vel = sincos(2π * rand(model.rng)) .* speed
        add_agent!(pos, model, vel, 1.0)
    end
    return model
end

#############
# Kickstart #
#############

model = ball_model()
agent_step!(agent, model) = move_agent!(agent, model, model.dt)

try
    abm_video(
        "bacteria.mp4",
        resolution=(1920, 1080),
        model,
        agent_step!;
        title = "Ball Model",
        frames = 50,
        spf = 2,
        framerate = 25,
    )
catch
    println("\nERROR: Failed to create video file.")
    println("NOTE:  Please run the 'integrated_gpu_support.sh' script.")
end

# End of File.
