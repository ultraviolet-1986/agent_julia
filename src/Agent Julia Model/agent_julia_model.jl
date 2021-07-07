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
Pkg.add("DrWatson")
Pkg.add("CairoMakie")
Pkg.add("InteractiveDynamics")
Pkg.add("Random")

using Agents
using DrWatson: @dict
using CairoMakie
using InteractiveDynamics
using Random

#################
# Prerequisites #
#################

# Share seed throughout program.
const seed = 41269

# Ensure results are reproducible.
Random.seed!(seed)

#############
# Variables #
#############

# Set temporal units.
const day = 24.0
const year = day * 365.0
const tend = year * 80.0

###########
# Structs #
###########

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
    days_mutated::Int              # number of days since mutated
    status::Symbol                 # :S, :I or :R
    mutation_probability::Float64  # May also be useless
end


#############
# Functions #
#############

# EDIT Comment to import S.I.R. model.
# function ball_model(; speed = 0.002)
#     space2d = ContinuousSpace((1, 1), 0.02)

#     model = ABM(
#         Wild_mtDNA,
#         space2d,
#         properties = Dict(:dt => 1.0),
#         rng = MersenneTwister(seed)
#     )

#     # And add some agents to the model.
#     # NOTE Change this to include the wild-type and mutant mtDNA.
#     for i in 1:200
#         pos = Tuple(rand(model.rng, 2))
#         vel = sincos(2π * rand(model.rng)) .* speed
#         add_agent!(pos, model, vel, 1.0)
#     end

#     return model
# end

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
        frames = 1000,
        spf = 2,
        framerate = 60
    )
catch
    println("\nERROR: Failed to create video file.")
    println("NOTE:  Please run the 'integrated_gpu_support.sh' script.")
end

# End of File.
