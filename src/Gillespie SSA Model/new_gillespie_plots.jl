#!/usr/bin/env julia

###########
# Imports #
###########

import Plots as plt

using DataFrames
using Distributions
using IterTools
using JLD2
using StatsPlots

#############
# Functions #
#############

function gillespie_mutation_load_plot()
    mutation_load_plot = plt.plot(
        (mean(mutation_vector) * 100),
        ylims=(0, 100),
        xlabel="Time (epochs) [80 years]",
        ylabel="Mean mutation load (%)",
        title="Gillespie SSA: 50/200 mutant mtDNA",
        legend=false,
        dpi=1200
    )

    plt.savefig(mutation_load_plot, "gillespie_mutation_load_plot.png")
end

#############
# Kickstart #
#############

# Reshape mutation_load data for iteration.
mutation_load = reshape(mutation_load, 1000, 960)

# Define vector to store previous matrix data.
mutation_vector = []

# Place each 'row' into the new vector as an element.
for i in eachrow(mutation_load)
    push!(mutation_vector, i)
end

gillespie_mutation_load_plot()

# End of File.
