#!/usr/bin/env julia

import Plots as plt

using DataFrames
using Distributions
using IterTools
using JLD2
using LinearAlgebra
using CSV

# Load variables from file and re-cast types.
@load "agent_data.jld2" data results
data = DataFrame(data)
results = DataFrame(results)

# agent_max = results."wild_status"[1][1] + results."mutant_status"[1][1]
agent_max = 250

# QUANTILE PLOT VARIABLES

steps = length(results."step"[1])

upper_quantile = []
middle_quantile = []
lower_quantile = []

for f in 1:1:steps
    push!(upper_quantile, quantile(Vector(data[f, 2:3]), 0.95)) # 95th percentile
    push!(middle_quantile, quantile(Vector(data[f, 2:3]), 0.5)) # 50th percentile
    push!(lower_quantile, quantile(Vector(data[f, 2:3]), 0.05)) #  5th percentile
end

quantiles = [upper_quantile, middle_quantile, lower_quantile] # N
# quantiles = (quantiles / 200) * 100                           #  %
# quantiles = (quantiles / 200)                                 # 0-1

#############
# Functions #
#############

function agent_mutation_load_plot()
    total = (sum(eachrow(results.wild_status)) + sum(eachrow(results.mutant_status)))
    total = collect(Iterators.flatten(total))

    mutation_load = sum(eachrow(results.mutant_status)) ./ total

    agents_mutation_load_plot = plt.plot(
        (mean(mutation_load) * 100),
        ylims=(0, 100),
        xlabel="Time (epochs) [80 years]",
        ylabel="Mean mutation load (%)",
        title="Agent Julia: 50/200 mutant mtDNA",
        legend=false,
        dpi=1200
    )

    plt.savefig(agents_mutation_load_plot, "agent_mutation_load_plot.png")
end


function agent_copy_levels_plot()
    mutation_load = sum(results."mutant_status") ./ length(results."mutant_status")

    total_copy_number = sum(results."wild_status" + results."mutant_status") /
        length(results."step")

    figure_2 = plt.plot(
        [total_copy_number, mutation_load],
        ylims=(0, agent_max),
        title="Patient with mutant mtDNA inheritance",
        xlabel="Time (epochs) [80 years]",
        ylabel="mtDNA levels (n)",
        label=["Total copy number" "Mutant copy number"],
        legend=:bottomright,
        dpi=1200
    )

    plt.savefig(figure_2, "agent_copy_levels_plot.png")
end


function agent_quantile_plot()
    figure_3 = plt.plot(
        quantiles,
        ylims=(0, agent_max),
        title="Patient with mutant mtDNA inheritance",
        xlabel="Time (epochs) [80 years]",
        ylabel="Mean agent levels (n)",
        label=["95th percentile" "50th percentile" "5th percentile"],
        dpi=1200,
    )

    plt.savefig(figure_3, "agent_quantile_plot.png")
end


function output_all_simulations_to_csv()
    for i in 1:length(results.step)
        wild = results[i, 2]
        mutant = results[i, 3]

        padding = Int(length(string(length(results.step))))

        mkpath("agent_csv_files")

        temp = DataFrame(wild_count=wild, mutant_count=mutant)

        fname = "$(lpad(i, padding, '0')).csv"
        CSV.write("agent_csv_files/$(fname)", temp)
    end
end

#############
# Kickstart #
#############

# agent_copy_levels_plot(); agent_mutation_load_plot(); agent_quantile_plot()

# End of File.
