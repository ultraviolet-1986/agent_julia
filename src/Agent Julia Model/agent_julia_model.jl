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
# TODO Have the model populated with evenly dispersed agents.
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
Pkg.add("DataFrames")
Pkg.add("Distributions")
Pkg.add("DrWatson")
Pkg.add("InteractiveDynamics")
Pkg.add("Plots")
Pkg.add("Random")

using Agents
using CairoMakie
using CSV
using DataFrames
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

month = 1.0
year = month * 12.0
day = year / 365.0
hour = day / 24.0
week = day * 7.0

tend = Int(round(year * 80.0))

δ = month

# COLOURS > TEXT OUTPUT

green = "\e[32m"
red = "\e[31m"
yellow = "\e[33m"
reset = "\e[0m"

# COLOURS > SIMULATION VIDEO

green_hex = "#338c54"  # Wild-type mtDNA
red_hex = "#bf2642"    # Mutant mtDNA

# GRAPH PROPERTIES

graph_width = 24.0
graph_height = 12.0

# AGENT PROPERTIES > INITIAL CONDITIONS

agent_min = 150
agent_max = rand(Poisson(200))

initial_mutants = 25

# SIMULATION PARAMETERS

loops = 10

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
end

#############
# Functions #
#############

function mutation_initiation(;
    interaction_radius = 0.012,
    dt = δ,
    reaction_rate = hour,
    N = agent_max,
    initial_mutated = initial_mutants,
    seed = random_seed,
)

    properties = @dict(
        reaction_rate,
        interaction_radius,
        dt,
    )

    space = ContinuousSpace((graph_width, graph_height), interaction_radius)

    model = ABM(
        mtDNA,
        space,
        properties=properties,
        rng=MersenneTwister(seed),
    )

    for ind in 1:N
        pos = (rand(Uniform(0.0, graph_width)), rand(Uniform(0.0, graph_height)))
        status = ind ≤ N - initial_mutated ? :W : :M
        mass = 1.0
        vel = sincos(2π * rand(model.rng)) .* δ

        add_agent!(pos, model, vel, mass, 0, status)
    end

    return model
end


function model_step!(model)
    # No action required.
end


function agent_step!(agent, model)
    move_agent!(agent, model, model.dt)
    random_action!(agent, model)
end


function random_action!(agent, model)
    roll = rand(Uniform(0.0, 1.0))

    copy_number = length(model.agents)

    # Replicate mtDNA
    if roll <= (1 / 2)
        if copy_number <= agent_max
            pos = agent.pos
            status = agent.status
            mass = 1.0
            vel = sincos(2π * rand(model.rng)) .* δ

            add_agent!(pos, model, vel, mass, 0, status)
        end

    # Degrade mtDNA
    else
        if copy_number >= agent_min
            kill_agent!(agent, model)
        end
    end

    # Mutate mtDNA
    if agent.status == :W
        if roll <= (1 / 10000)
            agent.status = :M
        end
    end
end


function perform_simulation()
    print("Performing simulation $(loops) time(s). Please wait... ")
    data = loop_simulation(loops)
    println("$(green)Done$(reset)")

    print("Writing data to $(yellow)agent_julia_data.csv$(reset) Please wait... ")
    CSV.write("agent_julia_data.csv", data)
    println("$(green)Done$(reset)")

    render_plot(data)

    return data
end


function loop_simulation(n=1::Int64)
    temp_results = []

    for i in 1:1:n
        wild(x) = count(i == :W for i in x)
        mutant(x) = count(i == :M for i in x)
        adata = [(:status, wild), (:status, mutant)]
        data, _ = run!(agent_julia_model, agent_step!, model_step!, tend; adata)
        push!(temp_results, data)
    end

    # Create final composite results.
    steps = eachcol(temp_results[1])[1]
    wild_counts = []
    mutant_counts = []

    for j in 1:1:length(temp_results)
        wild_count = eachcol(temp_results[j])[2]
        push!(wild_counts, wild_count)
    end

    for k in 1:1:length(temp_results)
        mutant_count = eachcol(temp_results[k])[3]
        push!(mutant_counts, mutant_count)
    end

    wild_mean = mean(wild_counts)
    mutant_mean = mean(mutant_counts)

    # Place final results in composite DataFrame.
    composite_results = DataFrame(
        steps = steps,
        wild_mean = wild_mean,
        mutant_mean = mutant_mean
    )

    return composite_results
end


function render_plot(data)
    times = eachcol(data)[1] / year    # Years
    mutation_level = eachcol(data)[3]  # N

    print("Building plot. Please wait... ")
    fig = plt.plot(
        times,
        mutation_level,
        xlims=(0, 80),
        ylims=(0, agent_max),
        title="mtDNA population dynamics (agent)",
        xlabel="Time (years)",
        ylabel="Mutation level (n)",
        legend=false,
        # smooth=true,
        dpi=1200,
    )
    println("$(green)Done$(reset)")

    print("Rendering plot to $(yellow)agent_julia_mutation_plot.png$(reset)... ")
    plt.savefig(fig, "agent_julia_mutation_plot.png")
    println("$(green)Done$(reset)")
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
            spf = Int(round(δ)),
            framerate = 60,
        )
        println("$(green)Done$(reset)\n")
    catch
        println("$(red)Error$(reset)")
        println("Please run the $(yellow)integrated_gpu_support.sh$(reset) script.\n")
    end
end


#############
# Kickstart #
#############

println("\n")

print("Creating mtDNA population dynamics model... ")
agent_julia_model = mutation_initiation()
println("$(green)Done$(reset)")

# NOTE Each method creates a separate simulation.
data = perform_simulation()
simulation_to_video()

# End of File.
