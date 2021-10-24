#!/usr/bin/env julia

###########
# Imports #
###########

using CSV
using DataFrames
using LinearAlgebra
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
# `create_full_copy_level_plot(data::DataFrame; model="agent"::String)`

Use the loaded model's CSV data to create plots. Note that the model
data should have the same shape regardless of which model is being
utilised.

## Examples

- `create_full_copy_level_plot(agent_data)`
- `create_full_copy_level_plot(gillespie_data; model="gillespie")`
"""
function create_full_copy_level_plot(data::DataFrame; model="agent"::String)
    filename = "$(model)_full_copy_level_plot.png"

    tend = length(data[:, 1][1])
    year = tend / 80.0
    decade = year * 10.0

    agent_max = last(sort(maximum(data.wild_count + data.mutant_count)))

    results = [mean((data.wild_count + data.mutant_count))]

    print("Rendering plot to $(yellow)$(filename)$(reset)... ")

    figure_1 = plot(
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

    figure_1 = plot!(
        results,
        xticks=(0:decade:tend, 0:10:80),
        ylims=(0, agent_max),
        xlims=(0, tend),
        linecolor=:red,
        legend=:bottomright,
        label=["Total mean" "Wild mean" "Mutant mean"],
    )

    savefig(figure_1, "$(filename)")

    println("$(green)Done$(reset)")
end


"""
# `create_mean_mtdna_level_plots(data::DataFrame; model="agent"::String)`

Use the loaded model's CSV data to create plots. Note that the model
data should have the same shape regardless of which model is being
utilised.

## Examples

- `create_mean_mtdna_level_plots(agent_data)`
- `create_mean_mtdna_level_plots(gillespie_data; model="gillespie")`
"""
function create_mean_mtdna_level_plots(data::DataFrame; model="agent"::String)
    filename = "$(model)_mean_mtdna_levels_plot.png"

    tend = length(data[:, 1][1])
    year = tend / 80.0
    decade = year * 10.0

    agent_max = last(sort(maximum(data.wild_count + data.mutant_count)))

    results = [
        mean((data.wild_count + data.mutant_count)),
        mean(data.wild_count),
        mean(data.mutant_count)
    ]

    print("Rendering plot to $(yellow)$(filename)$(reset)... ")

    figure_1 = plot(
        results,
        title="$(titlecase(model)) model",
        xlabel="Patient age (y)",
        ylabel="mtDNA copy numbers",
        xticks=(0:decade:tend, 0:10:80),
        ylims=(0, agent_max),
        label=["Total mean" "Wild mean" "Mutant mean"],
        dpi=1200
    )

    savefig(figure_1, "$(filename)")

    println("$(green)Done$(reset)")
end


"""
# `create_cross_section_plot(data::DataFrame; model="agent"::String)`

Use the loaded model's CSV data to create plots. Note that the model
data should have the same shape regardless of which model is being
utilised.

## Examples

- `create_cross_section_plot(agent_data)`
- `create_cross_section_plot(gillespie_data; model="gillespie")`
"""
function create_cross_section_plot(data::DataFrame; model="agent"::String)
    filename = "$(model)_mtdna_cross_section_plot.png"

    totals = sum(data.wild_count + data.mutant_count)
    mutation_load = sum(data.mutant_count) ./ totals

    print("Rendering plot to $(yellow)$(filename)$(reset)... ")

    figure_1 = plot(
        (mutation_load * 100),
        title="$(titlecase(model)) model",
        xlabel="Patient age (y)",
        ylabel="Mutation load (%)",
        xlims=(30, 40),
        xticks=(30:2:40),
        ylims=(0, 100),
        yticks=(0:25:100),
        label=false,
        dpi=1200
    )

    savefig(figure_1, "$(filename)")

    println("$(green)Done$(reset)")
end


"""
# `create_percentile_plot(data::DataFrame; model="agent"::String)`

Use the loaded model's CSV data to create plots. Note that the model
data should have the same shape regardless of which model is being
utilised.

## Examples

- `create_percentile_plot(agent_data)`
- `create_percentile_plot(gillespie_data; model="gillespie")`
"""
function create_percentile_plot(data::DataFrame; model="agent"::String)
    filename = "$(model)_mtdna_percentile_plot.png"

    tend = length(data[:, 1][1])
    year = tend / 80.0
    decade = year * 10.0

    agent_max = last(sort(maximum(data.wild_count + data.mutant_count)))

    total_counts = (data.wild_count + data.mutant_count)

    upper_quantile = quantile.(total_counts, 0.95)
    middle_quantile = quantile.(total_counts, 0.5)
    lower_quantile = quantile.(total_counts, 0.05)

    results = [
        upper_quantile,
        middle_quantile,
        lower_quantile
    ]

    print("Rendering plot to $(yellow)$(filename)$(reset)... ")

    figure_1 = plot(
        results,
        title="$(titlecase(model)) model",
        xlabel="Patient age (y)",
        ylabel="mtDNA copy levels",
        xticks=(0:decade:tend, 0:10:80),
        ylims=(0, agent_max),
        legend=:bottomright,
        label=["95th Percentile" "50th Percentile" " 5th Percentile"],
        dpi=1200
    )

    savefig(figure_1, "$(filename)")

    println("$(green)Done$(reset)")
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

create_full_copy_level_plot(agent_data)
create_full_copy_level_plot(gillespie_data; model="gillespie")

create_mean_mtdna_level_plots(agent_data)
create_mean_mtdna_level_plots(gillespie_data; model="gillespie")

create_cross_section_plot(agent_data)
create_cross_section_plot(gillespie_data; model="gillespie")

create_percentile_plot(agent_data)
create_percentile_plot(gillespie_data; model="gillespie")

# End of File.
