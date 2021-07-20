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

#########
# Notes #
#########

# TODO Write event-trigger timing mechanism.
# TODO Sort agents and choose the top to pass to said mechanism, and
#      randomly select the event to occur: replicate, degrade, mutate.
#      Note that agents are to be sorted by their 'time of interaction'
#      being the shortest.
# TODO Perform the simulation, then optionally output to video later.
#      The important thing is to be able to access the data through the
#      REPL.
# NOTE Wild-type and mutants have a slightly different mass,
#      investigate this.

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

random_seed = 41269
Random.seed!(random_seed)

#############
# Variables #
#############

# TEMPORAL UNITS

# NOTE Stepping functions require type 'Int'.

hour = Int(1)
day = Int(hour * 24)
year = Int(day * 365)
month = Int(year / 12)

tend = Int(year * 80)

# COLOURS

green = "\e[32m"
red = "\e[31m"
yellow = "\e[33m"
reset = "\e[0m"

green_hex = "#338c54"  # Wild-type mtDNA
red_hex = "#bf2642"    # Mutant mtDNA

# AGENT PROPERTIES

agent_max = 200

###########
# Structs #
###########

mutable struct mtDNA <: AbstractAgent
    id::Int
    pos::NTuple{2,Float64}
    vel::NTuple{2,Float64}
    mass::Float64
    days_mutated::Int  # Days since agent is mutated.
    status::Symbol     # :W (Wild-type mtDNA), :M (Mutant mtDNA)
    β::Float64         # Mutation probability.
end

#############
# Functions #
#############

function mutation_initiation(;
    mutation_probability = 0.1, # 0.05,
    isolated = 0.0,
    interaction_radius = 0.012,
    dt = 1.0,
    speed = hour, # 0.002,
    death_rate = day * 260.0, # 0.044,
    N = agent_max,        # Set N to not go above 'agent_max'.
    initial_mutated = 10,  # Tied to initial conditions.
    seed = random_seed,
    βmin = 0.0,
    βmax = 0.2,
)

    properties = @dict(
        mutation_probability,
        death_rate,
        interaction_radius,
        dt,
    )

    # space = ContinuousSpace((10, 10), 0.02)
    space = ContinuousSpace((24, 12), interaction_radius)

    model = ABM(
        mtDNA,
        space,
        properties=properties,
        rng=MersenneTwister(seed)
    )

    # Add initial individuals
    for ind in 1:N
        pos = Tuple(rand(model.rng, 2))
        status = ind ≤ N - initial_mutated ? :W : :M
        isisolated = ind ≤ isolated * N
        mass = isisolated ? Inf : 1.0
        vel = isisolated ? (0.0, 0.0) : sincos(2π * rand(model.rng)) .* speed

        β = (βmax - βmin) * rand(model.rng) + βmin
        add_agent!(pos, model, vel, mass, 0, status, β)
    end

    return model
end


function model_step!(model)
    r = model.interaction_radius

    for (a1, a2) in interacting_pairs(model, r, :nearest;)  # Errors here on n > 1.
        elastic_collision!(a1, a2, :mass)
        println("Bounce!")
        println("Agents remaining: $(length(model.agents)), $(a1.status), $(a2.status)")
    end
end


function agent_step!(agent, model)
    move_agent!(agent, model, model.dt)
    random_action!(agent, model)
end


function random_action!(agent, model)
    roll = rand()
    mother_position = agent.pos

    # Degrade mtDNA
    if roll <= (1 / 3)
        # Prevent population extinction
        if Int(length(model.agents)) > 50
            kill_agent!(agent, model)
            println("Degrade mtDNA")
        end

    elseif roll <= (2 / 3)
        # Prevent population explosion
        if Int(length(model.agents)) < agent_max

            # Replicate wild-type mtDNA
            if agent.status == :W
                add_agent!(agent, mother_position, model)
                println("Replicate wild-type mtDNA")

            # Replicate mutant mtDNA
            else
                add_agent!(agent, mother_position, model)
                agent.status = :M
                agent.days_mutated = 0
                println("Replicate mutant mtDNA")
            end
        end

    # Mutate wild-type mtDNA
    else
        if agent.status == :W
            agent.status = :M
            agent.days_mutated = 0
            println("Mutate wild-type mtDNA")
        end
    end

    # Update agent
    if agent.status == :M
        agent.days_mutated += 1
    end

    println("Agents remaining: $(length(model.agents)), $(agent.status)")
    return
end


function render_video()
    # RENDER ONLY
    simulation = abm_video(
        "agent_julia_simulation.mp4",
        agent_julia_model,
        agent_step!,
        model_step!;
        title = "mtDNA population dynamics",
        frames = 100, # tend, # frames = 1000,
        ac = model_colours,
        as = 10,
        spf = 1, #month, # spf = 1,
        framerate = 60,
    )

    # # RUN ONLY
    # simulation = run!(
    #     agent_julia_model,
    #     agent_step!,
    #     model_step!,
    #     30;
    # )

    # Can be accessed from REPL by using below:
    # model[1].pos  # Get first agent's position
    # sir_model[1].mass  # Get SIR mass

    return simulation
end

#############
# Kickstart #
#############

println("\n")

print("Defining simulation colour palette... ")
model_colours(a) = a.status == :W ? green_hex : red_hex
println("$(green)Done$(reset)")

print("Creating mtDNA population dynamics model... ")
agent_julia_model = mutation_initiation()
println("$(green)Done$(reset)")

# print("Rendering simulation output as 'agent_julia_simulation.mp4'... ")
# try
#     simulation = render_video()
#     println("$(green)Done$(reset)")
# catch
#     println("$(red)Error$(reset)")
#     println("$(yellow)Please run the 'integrated_gpu_support.sh' script.$(reset)")
# end
simulation = render_video()

# End of File.
