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

hour = 1
day = hour * 24
year = day * 365

tend = year * 80

# COLOURS

green = "\e[32m"
red = "\e[31m"
yellow = "\e[33m"
reset = "\e[0m"

grey_hex = "#2b2b33"   # Wild-type mtDNA
green_hex = "#338c54"  # May remove this, mutants cannot recover.
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
    status::Symbol     # :S, :I or :R, This will change.
    β::Float64         # Mutation probability.
end

#############
# Functions #
#############

function mutation_initiation(;
    infection_period = 30 * day,     # Irrelevant as above for this model.
    detection_time = 14 * day,       # Irrelevant as above for this model.
    reinfection_probability = 0.05,  # Irrelevant as above for this model.
    isolated = 0.0,                  # Irrelevant as above for this model.
    interaction_radius = 0.012,      # Irrelevant as above for this model.
    dt = 1.0,
    speed = 0.002,
    death_rate = 0.044,   # Change to factor in degradation.
    N = agent_max,        # Set N to not go above 'agent_max'.
    initial_mutated = 1,  # Tied to initial conditions.
    seed = random_seed,
    βmin = 0.4,
    βmax = 0.8,
)

    properties = @dict(
        infection_period,
        reinfection_probability,
        detection_time,
        death_rate,
        interaction_radius,
        dt,
    )

    space = ContinuousSpace((16, 9), 0.02)

    model = ABM(
        mtDNA,
        space,
        properties=properties,
        rng=MersenneTwister(random_seed)
    )

    # Add initial individuals
    for ind in 1:N
        pos = Tuple(rand(model.rng, 2))
        status = ind ≤ N - initial_mutated ? :S : :I
        isisolated = ind ≤ isolated * N
        mass = isisolated ? Inf : 1.0
        vel = isisolated ? (0.0, 0.0) : sincos(2π * rand(model.rng)) .* speed

        # TODO Trigger event here
        # very high transmission probability
        # we are modelling close encounters after all
        β = (βmax - βmin) * rand(model.rng) + βmin
        add_agent!(pos, model, vel, mass, 0, status, β)
    end

    return model
end


function transmit!(a1, a2, rp)
    # for transmission, only 1 can have the disease (otherwise nothing happens)
    count(a.status == :I for a in (a1, a2)) ≠ 1 && return
    infected, healthy = a1.status == :I ? (a1, a2) : (a2, a1)

    rand(sir_model.rng) > infected.β && return

    if healthy.status == :R
        rand(sir_model.rng) > rp && return
    end
    healthy.status = :I
end


function random_action!(agent, model)
    roll = rand()

    # Replicate
    if roll <= (1 / 3)
        mother_position = agent.pos
        add_agent!(agent, mother_position, model)

    # Degrade
    elseif roll <= (2 / 3)
        kill_agent!(agent, model)

    # Mutate
    else
        agent.status = :M
    end
end


function sir_model_step!(model)
    r = model.interaction_radius

    for (a1, a2) in interacting_pairs(model, r, :nearest)
        transmit!(a1, a2, model.reinfection_probability)
        elastic_collision!(a1, a2, :mass)
    end
end


function sir_agent_step!(agent, model)
    move_agent!(agent, model, model.dt)
    update!(agent)
    replicate_or_degrade!(agent, model)
end


update!(agent) = agent.status == :I && (agent.days_mutated += 1)


function replicate_or_degrade!(agent, model)
    if agent.days_mutated ≥ model.infection_period

        # Degrade agent
        if rand(model.rng) ≤ model.death_rate
            kill_agent!(agent, model)  # Degrade.

        # Replicate agent
        else
            # Update original agent.
            # NOTE Original code block.
            agent.status = :R       # Change status of agent.
            agent.days_mutated = 0

            # Get mother's position for placing daughter instance.
            mother_position = agent.pos

            # Create replicated agent.
            # NOTE Additional code block.
            # NOTE Replicated agent will be placed randomly adjacent as
            #      overlapping is not possible with this model's
            #      physics.
            new_agent = add_agent!(agent, mother_position, model)
            new_agent.status = :R  # Set to recovered colour (green).
            new_agent.days_mutated = 0
        end
    end
end


function render_video()
    # RENDER ONLY
    abm_video(
        "agent_julia_simulation.mp4",
        sir_model,
        sir_agent_step!,
        sir_model_step!;
        title = "mtDNA population dynamics",
        frames = tend,
        ac = sir_colors,
        as = 10,
        spf = 1,
        framerate = 60,
    )

    # RUN ONLY
    # run!(
    #     sir_model,
    #     sir_agent_step!,
    #     sir_model_step!,
    #     1000;
    # )

    # Can be accessed from REPL by using below:
    # model[1].pos  # Get first agent's position
    # sir_model[1].mass  # Get SIR mass
end

#############
# Kickstart #
#############

println("\n")

print("Defining simulation colour palette... ")
sir_colors(a) = a.status == :S ? grey_hex : a.status == :I ? red_hex : green_hex
println("$(green)Done$(reset)")

print("Creating mtDNA population dynamics model... ")
sir_model = mutation_initiation()
println("$(green)Done$(reset)")

print("Rendering simulation output as 'agent_julia_simulation.mp4'... ")
try
    simulation = render_video()
    println("$(green)Done$(reset)")
catch
    println("$(red)Error$(reset)")
    println("$(yellow)Please run the 'integrated_gpu_support.sh' script.$(reset)")
end

# End of File.
