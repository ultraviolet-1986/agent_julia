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

# TODO Slow the kinetic rates by a factor of 1000, or increase the
#      speed of the model to account for the higher values. Temporal
#      rates MUST be integer for the stepping function.
# TODO Perform the simulation, then optionally output to video later.
#      The important thing is to be able to access the data through the
#      REPL.
# NOTE Wild-type and mutants have a slightly different mass,
#      investigate this.
# TODO Take kinetic rates from the Gillespie model and translate them
#      for this model. Reduce the mutation rate.

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
Pkg.add("CSV")
Pkg.add("Distributions")
Pkg.add("DrWatson")
Pkg.add("InteractiveDynamics")
Pkg.add("Random")

using Agents
using CairoMakie
using CSV
using Distributions
using DrWatson: @dict
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

# NOTE Stepping functions require type 'Int'. Will account for this.

hour = 1.0 / 100.0
day = hour * 24.0
year = day * 365.0
month = year / 12.0

tend = Int(round(year * 80.0))

δ = month

# COLOURS (USED FOR VIDEO RENDER ONLY)

green = "\e[32m"
red = "\e[31m"
yellow = "\e[33m"
reset = "\e[0m"

green_hex = "#338c54"  # Wild-type mtDNA
red_hex = "#bf2642"    # Mutant mtDNA

# AGENT PROPERTIES

agent_max = rand(Poisson(200))
initial_mutants = 10

βmin = 0.0
βmax = 0.0001

λ = day * 260.0

###########
# Structs #
###########

mutable struct mtDNA <: AbstractAgent
    id::Int
    pos::NTuple{2,Float64}
    vel::NTuple{2,Float64}
    mass::Float64
    age_in_days::Int
    status::Symbol  # :W (Wild-type mtDNA), :M (Mutant mtDNA)
    β::Float64      # Mutation probability.
end

#############
# Functions #
#############

function mutation_initiation(;
    mutation_probability = 0.05,
    isolated = 0.0,
    interaction_radius = 0.012,
    dt = δ,
    speed = δ / 500,
    death_rate = λ,  # Terminates simulation itself at step 6240 (λ).
    N = agent_max,
    initial_mutated = initial_mutants,
    seed = random_seed,
    βmin = 0.4,
    βmax = 0.8,
)

    properties = @dict(
        mutation_probability,
        death_rate,
        interaction_radius,
        dt,
    )

    space = ContinuousSpace((24, 12), interaction_radius)

    model = ABM(
        mtDNA,
        space,
        properties=properties,
        rng=MersenneTwister(seed)
    )

    # Add initial agents
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

    for (a1, a2) in interacting_pairs(model, r, :nearest;)
        elastic_collision!(a1, a2, :mass)
    end
end


function agent_step!(agent, model)
    move_agent!(agent, model, model.dt)
    agent.age_in_days += 1  # Update agent's age
    random_action!(agent, model)
end


function random_action!(agent, model)
    roll = rand()
    mother_position = agent.pos

    # Degrade mtDNA
    if roll <= (1 / 3)
        if agent.age_in_days ≥ λ
            kill_agent!(agent, model)
        end

    # Replicate mtDNA
    elseif roll <= (2 / 3)

        # Wild-type mtDNA
        if agent.status == :W
            add_agent!(agent, mother_position, model)
            agent.status = :W
            agent.age_in_days = 0

        # Mutant mtDNA
        else
            add_agent!(agent, mother_position, model)
            agent.status = :M
            agent.age_in_days = 0
        end

    # Mutate wild-type mtDNA
    else
        if agent.status == :W
            # Calculate probability of mutation
            β = rand(agent_julia_model.rng)

            # Mutate if calculated probability falls within range
            if β > βmin && β < βmax
                agent.status = :M
            end
        end
    end
end


function render_plot(data)
    figure = Figure()

    ax = figure[1, 1] = Axis(
        figure;
        title = "mtDNA population dynamics",
        ylabel = "Mutation level",
        xlabel = "Time",
    )

    mutation_line = lines!(
        ax,
        data[:, dataname((:status, mutant))],
        color = :red
    )

    save("agent_julia_plot.svg", figure)

    return
end


function render_video()
    abm_video(
        "agent_julia_simulation.mp4",
        agent_julia_model,
        agent_step!,
        model_step!;
        title = "mtDNA population dynamics",
        frames = tend, # frames = 1000,
        ac = model_colours,
        as = 10,
        spf = 1, # month, # spf = 1,
        framerate = 60,
    )
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

# Get counts of elements
wild(x) = count(i == :W for i in x)
mutant(x) = count(i == :M for i in x)
adata = [(:status, wild), (:status, mutant)]

# Save simulation data
# data, _ = run!(agent_julia_model, agent_step!, model_step!, tend; adata)

# total_wild = eachcol(data)[2]
# total_mutants = eachcol(data)[3]

# CSV.write("agent_julia_data.csv", data)

# render_plot(data)

# End of File.
