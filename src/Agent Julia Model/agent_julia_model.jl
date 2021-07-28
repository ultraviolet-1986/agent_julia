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

# TODO Discover how to include a hard boundary.
# TODO Adjust space to more closely match a muscle fibre cell.
# TODO Include measure for discretising time so we have 960 states in
#      months rather than recording every step. This must be tied to δ.
# TODO Have the model populated with evenly dispersed agents.
# TODO Include positional data as part of the CSV export.
# TODO Export each frame as PNG and use FFMPEG to stitch together a
#      video and remove the 'abm_video' function call. This helps to
#      ensure plots and states match that simulation.

##############
# References #
##############

# - [Agents.jl Tutorial](https://git.io/Jc1w6)

###########
# Imports #
###########

import Pkg
import Plots as plt

Pkg.add("Agents")
Pkg.add("CairoMakie")
Pkg.add("CSV")
Pkg.add("Distributions")
Pkg.add("DrWatson")
Pkg.add("InteractiveDynamics")
Pkg.add("Plots")
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

# NOTE Stepping functions require type 'Int'.

time_factor = 100.0

hour = 1.0 / time_factor
day = hour * 24.0
week = day * 7.0
fortnight = week * 2.0
year = day * 365.0
month = year / 12.0

tend = Int(round(year * 80.0))

δ = fortnight

# COLOURS (TEXT OUTPUT ONLY)

green = "\e[32m"
red = "\e[31m"
yellow = "\e[33m"
reset = "\e[0m"

# COLOURS (USED FOR VIDEO RENDER ONLY)

green_hex = "#338c54"  # Wild-type mtDNA
red_hex = "#bf2642"    # Mutant mtDNA

# AGENT PROPERTIES

agent_max = rand(Poisson(200))
initial_mutants = 10

# Mutant replication boundaries.
βmin = 0.0
βmax = 0.0001

# mtDNA half-life (death rate).
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
    death_rate = day,
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
    agent.age_in_days += 1
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
            β = rand(agent_julia_model.rng)

            if β > βmin && β < βmax
                agent.status = :M
            end
        end
    end
end


function perform_simulation()
    wild(x) = count(i == :W for i in x)
    mutant(x) = count(i == :M for i in x)

    adata = [(:status, wild), (:status, mutant)]

    print("Performing simulation. Please wait... ")
    data, _ = run!(agent_julia_model, agent_step!, model_step!, tend; adata)
    println("$(green)Done$(reset)")

    CSV.write("agent_julia_data.csv", data)

    render_plot(data)

    return data
end


function render_plot(data)
    x = eachcol(data)[1] / time_factor  # Times (years)
    y = eachcol(data)[3]                # Mutation level (n)

    print("Building plot. Please wait... ")
    fig = plt.plot(
        x,
        y,
        xlims=(0, 80),
        title="mtDNA population dynamics (agent)",
        xlabel="Time (years)",
        ylabel="Mutation level (n)",
        legend=false,
        dpi=1200,
    )
    println("$(green)Done$(reset)")

    print("Rendering plot to $(yellow)agent_julia_mutation_plot.png$(reset)... ")
    plt.savefig(fig, "agent_julia_mutation_plot.png")
    println("$(green)Done$(reset)\n")
end


function simulation_to_video()
    print("Defining simulation colour palette... ")
    model_colours(a) = a.status == :W ? green_hex : red_hex
    println("$(green)Done$(reset)")

    print("Rendering simulation output as $(yellow)agent_julia_simulation.mp4$(reset)... ")

    try
        abm_video(
            "agent_julia_simulation.mp4",
            agent_julia_model,
            agent_step!,
            model_step!;
            title = "mtDNA population dynamics",
            frames = tend,
            ac = model_colours,
            as = 10,
            spf = 1,
            framerate = 60,
        )
        println("$(green)Done$(reset)")
    catch
        println("$(red)Error$(reset)")
        println("Please run the $(yellow)integrated_gpu_support.sh$(reset) script.")
    end
end


#############
# Kickstart #
#############

println("\n")

print("Creating mtDNA population dynamics model... ")
agent_julia_model = mutation_initiation()
println("$(green)Done$(reset)")

# NOTE Use only one method below.
simulation_to_video()
# data = perform_simulation()

# End of File.
