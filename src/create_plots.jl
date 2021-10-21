#!/usr/bin/env julia

###########
# Imports #
###########

using CSV
using DataFrames
using Plots
using Statistics
using StatsBase

#################
# Prerequisites #
#################

# plotly()

#############
# Variables #
#############

# COLOURS > TEXT OUTPUT

green = "\e[32m"
red = "\e[31m"
yellow = "\e[33m"
reset = "\e[0m"

#############
# Functions #
#############


"""
# `load_csv_data(model="agent"::String)`

Load a model's CSV files to memory and return a data frame. Note that
by default, this function will attempt to load the Agent Julia model's
data.

## Examples

- `load_csv_data()`
- `load_csv_data("agent")`
- `load_csv_data("gillespie")`
"""
function load_csv_data(model="agent"::String)
    path = "$(Base.source_dir())/data/$(model)_csv_files"
    data = []

    if isdir("$(path)")
        print("Reading CSV files in directory $(yellow)$(model)_csv_files$(reset)... ")
        try
            for file in readdir(path; join=true)
                csv_data = CSV.File("$(file)")
                push!(data, csv_data)
            end
            println("$(green)Done$(reset)")
            return DataFrame(data)
        catch
            println("$(red)Error$(reset)")
            error("No files detected or invalid data structure.")
            return
        end
    else
        error("Directory $(yellow)$(model)_csv_files$(reset) does not exist.")
        return
    end
end


"""
# `create_model_plots(data::DataFrame; model="agent"::String)`

Use the loaded model's CSV data to create plots. Note that the model
data should have the same shape regardless of which model is being
utilised.

## Examples

- `create_model_plots(agent_data)`
- `create_model_plots(gillespie_data; model="gillespie")`
"""
function create_model_plots(data::DataFrame; model="agent"::String)
    println("$(titlecase(model))")
    println("$(typeof(data))\n")

    tend = length(data[:, 1][1])
    year = tend / 80.0
    decade = year * 10.0

    agent_max = last(sort(maximum(data.wild_count + data.mutant_count)))

    totals = data.wild_count + data.mutant_count
    # results = mean.(normalize(totals) * agent_max)  # Normalized totals %

    upper_percentile = percentile.(totals, 95)
    middle_percentile = percentile.(totals, 50)
    lower_percentile = percentile.(totals, 5)

    percentiles = [upper_percentile, middle_percentile, lower_percentile]

    total_mean = mean.(totals)

    fig1 = plot(
        (data.wild_count + data.mutant_count),
        title="$(titlecase(model)) model",
        xlabel="Time (years)",
        ylabel="mtDNA copy numbers",
        xticks=(0:decade:tend, 0:10:80),
        ylims=(0, agent_max),
        legend=false,
        label="",
        linecolor=:black,
        linealpha=0.01,
        dpi=1200
    )

    fig1 = plot!(
        results,
        xticks=(0:decade:tend, 0:10:80),
        ylims=(0, agent_max),
        xlims=(0, tend),
        palette = :Dark2_5,
        legend=:bottomright,
        label=["Total mean" "Wild mean" "Mutant mean"],
        # smooth=true,
    )

    savefig(fig1, "$(model)_figure_1.png")

    # mutation_load =

    # fig2 = plot(
    #     (data.wild_count + data.mutant_count),
    #     title="$(titlecase(model)) model",
    #     xlabel="Time (years)",
    #     ylabel="mtDNA counts (n)",
    #     xticks=(0:decade:tend, 0:10:80),
    #     ylims=(0, agent_max),
    #     legend=false,
    #     label="",
    #     linecolor=:black,
    #     linealpha=0.01,
    #     dpi=1200
    # )
end

#############
# Kickstart #
#############

# KICKSTART > OPTIONAL

# using JLD2
# @save "full_data.jld2" agent_data gillespie_data
# @load "full_data.jld2" agent_data gillespie_data

# KICKSTART > LOAD DATA

agent_data = load_csv_data()
gillespie_data = load_csv_data("gillespie")

# KICKSTART > PLOT DATA

create_model_plots(agent_data)
create_model_plots(gillespie_data; model="gillespie")

# End of File.
