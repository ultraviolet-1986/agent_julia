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

using Agents
using CairoMakie
using CSV
using DataFrames
using Distributions
using DrWatson: @dict
using InteractiveDynamics
using Random

using IterTools
using Statistics

#################
# Prerequisites #
#################

random_seed = 41269
Random.seed!(random_seed)

#############
# Variables #
#############

# TEMPORAL UNITS

day = Float64(24.0 * 3600.0)  # Seconds per day.
week = Float64(day * 7.0)         # Seconds per week.
year = Float64(day * 365.0)       # Seconds per year.
month = Float64(year / 12.0)      # Seconds per month.
hour = Float64(day / 24.0)        # Seconds per hour.

# Reaction rate
λ = 260.0
reaction_rate = log(2) / λ
reaction_rate = reaction_rate * 365 / 12

# TODO Define replication, degradation, and mutation rates seperately.

δ = month

# NOTE Stepping functions require type 'Int'.

tend = Int(round(year * 80.0 / δ))  # 1 step is equal to 1 month.

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
graph_height = 6.0

# INITIAL CONDITIONS

agent_min = 150
agent_max = 250

initial_wild = 150
initial_mutants = 50

agent_total = initial_mutants + initial_wild

# SIMULATION PARAMETERS

loops = 1000

# PATHS

data_path = "$(Base.source_dir())/data/agent_csv_files"

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


"""
`mutation_initiation(;
    interaction_radius = 0.012,
    dt = δ,
    reaction_rate = hour,
    N = agent_max,
    initial_mutated = initial_mutants,
    seed = random_seed,
)`

Create and populate the Agent Julia model with specified agents.
"""
function mutation_initiation(;
    interaction_radius = 0.012,
    dt = δ,
    reaction_rate = reaction_rate,
    N = agent_total,
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
        status = ind <= N - initial_mutated ? :W : :M
        mass = 1.0
        vel = sincos(2π * rand(model.rng)) .* δ

        add_agent!(pos, model, vel, mass, 0, status)
    end

    return model
end


"""
`model_step!(model)`

Perform all model actions at currrent step.

NOTE: Agent Julia does not require model actions.
"""
function model_step!(model); end


"""
`agent_step!(agent, model)`

Perform all agent actions at currrent step.
"""
function agent_step!(agent, model)
    move_agent!(agent, model, model.dt)
    random_action!(agent, model)
end


"""
`random_action!(agent, model)`

Randomly select an action: Replicate, Degrade, or Mutate an agent.
"""
function random_action!(agent, model)
    roll = rand(Uniform(0.0, 1.0))

    copy_number = nothing
    copy_number = length(model.agents)

    # Replicate mtDNA
    if roll <= reaction_rate
        if copy_number <= agent_max
            pos = agent.pos
            status = agent.status
            mass = 1.0
            vel = sincos(2π * rand(model.rng)) .* δ

            add_agent!(pos, model, vel, mass, 0, status)
        end

    # Degrade mtDNA
    elseif (roll > reaction_rate) && (roll <= (reaction_rate * 2))  # Change var names
        if copy_number >= agent_min
            kill_agent!(agent, model)
        end
    end

    # # Mutate mtDNA
    # if agent.status == :W
    #     if roll <= (1 / 10000)
    #         agent.status = :M
    #     end
    # end
end


"""
`perform_simulation()`

Run the simulation through the `loop_simulation(n)` function `n`
time(s).
"""
function perform_simulation()
    print("Performing simulation $(yellow)$(loops)$(reset) time(s). Please wait... ")
    data = loop_simulation(loops)
    println("$(green)Done$(reset)")

    return data
end


"""
`loop_simulation(n=1::Int64)`

Perform the Agent Julia simulation `n` number of times.
"""
function loop_simulation(n=1::Int64)
    temp_results = []

    for i in 1:1:n
        # Purge and redefine model for new simulation.
        agent_julia_model = nothing
        agent_julia_model = mutation_initiation()

        wild(x) = count(i == :W for i in x)
        mutant(x) = count(i == :M for i in x)

        # Purge and redefine 'adata' for new simulation.
        adata = nothing
        adata = [(:status, wild), (:status, mutant)]

        data, _ = run!(agent_julia_model, agent_step!, model_step!, tend; adata)
        push!(temp_results, data)
    end

    # Create final composite results.
    steps = eachcol(temp_results[1])[1]
    wild_counts = []
    mutant_counts = []
    total_counts = []

    for j in 1:1:length(temp_results)
        wild_count = eachcol(temp_results[j])[2]
        mutant_count = eachcol(temp_results[j])[3]
        total_cn = eachcol(temp_results[j])[2] + eachcol(temp_results[j])[3]

        push!(wild_counts, wild_count)
        push!(mutant_counts, mutant_count)
        push!(total_counts, total_cn)
    end

    wild_mean   = mean(wild_counts)
    mutant_mean = mean(mutant_counts)
    total_mean  = mean(total_counts)

    composite_results = DataFrame(
        steps = steps,
        wild_mean = wild_mean,
        mutant_mean = mutant_mean,
        total_mean = total_mean
    )

    complete_results = temp_results

    return composite_results, temp_results
end


"""
`output_simulation_to_video()`

Run a single simulation and output the results to video.
"""
function output_simulation_to_video()
    video_path = "$(Base.source_dir())/videos"
    mp4_path = "$(video_path)/agent_julia_simulation.mp4"

    rm(video_path; force=true, recursive=true)
    mkpath(video_path)

    print("\nDefining simulation colour palette... ")
    model_colours(a) = a.status == :W ? green_hex : red_hex
    println("$(green)Done$(reset)")

    agent_julia_model = nothing
    agent_julia_model = mutation_initiation()

    print("Rendering simulation output as $(yellow)$(mp4_path)$(reset)... ")
    try
        abm_video(
            "$(mp4_path)",
            agent_julia_model,
            agent_step!,
            model_step!;
            title = "Patient with mutant inheritance",
            frames = tend,
            ac = model_colours,
            as = 10,
            spf = Int(round(δ)),
            framerate = 60,
        )
        println("$(green)Done\n$(reset)")
    catch
        println("$(red)Error$(reset)")
        println("Please run the $(yellow)integrated_gpu_support.sh$(reset) script.\n")
    end
end


#############
# Kickstart #
#############

println("\nAGENT JULIA SIMULATION")
println("======================\n")

println("‣ Maximum possible agents....: $(green)$(agent_max)$(reset)")
println("‣ Total agents...............: $(green)$(agent_total)$(reset)")
println("‣ Wild-type agents...........: $(green)$(initial_wild)$(reset)")
println("‣ Mutant elements............: $(green)$(initial_mutants)$(reset)")
println("‣ Simulation loops...........: $(green)$(loops)$(reset)\n")

print("Creating mtDNA population dynamics model... ")
agent_julia_model = nothing
agent_julia_model = mutation_initiation()
println("$(green)Done$(reset)")

# Execute Agent Julia simulation.
data, results = perform_simulation()


# KICKSTART > DATA OUTPUT

rm(data_path; force=true, recursive=true)
mkpath(data_path)

print("Writing all data to $(yellow)$(loops)$(reset) CSV file(s)... ")
results = DataFrame(results)
for i in 1:1:loops
    wild = results[i, 2]
    mutant = results[i, 3]

    temp = DataFrame(wild_count=wild, mutant_count=mutant)

    padding = Int(length(string(length(results.step))))

    fname = "$(lpad(i, padding, '0')).csv"
    CSV.write("$(data_path)/$(fname)", temp)
end
println("$(green)Done$(reset)")


# KICKSTART > PROMPT ADDITIONAL STEPS

println("""\nOPTIONAL: Use $(yellow)output_simulation_to_video()$(reset) function.""")
println("""Plots may be created using the R script $(yellow)create_plots.R$(reset).\n""")

# End of File.
