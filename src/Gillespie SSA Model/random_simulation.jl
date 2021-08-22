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

###########
# Imports #
###########

using CSV
using DataFrames
using Distributions
using Plots
using Random
using StatsPlots

#################
# Prerequisites #
#################

# Ensure results are not reproducible.
Random.seed!()

#############
# Variables #
#############

# VARIABLES > INITIAL CONDITIONS

# Randomly select total number of mtDNA elements using the Poisson
# distribution, then select the number of wild elements and assign the
# remainder as mutant.

max_elements = rand(Poisson(200))

# Randomly assign amounts and prevent mutation levels greater than
# approximately 50%.

random_wild = rand(120:max_elements)
random_mutant = max_elements - random_wild

# Initial concentrations of wild-type and mutant mtDNA.
u0 = [random_wild, random_mutant]

# VARIABLES > PATHS

data_path = "$(Base.source_dir())/data"
plot_path = "$(Base.source_dir())/plots/random_simulation"

plot_1_path = "$(plot_path)/01_timeline.png"
plot_2_path = "$(plot_path)/02_quantiles.png"
plot_3_path = "$(plot_path)/03_density.png"

csv_path = "$(data_path)/random_simulation.csv"

#############
# Kickstart #
#############

mkpath(data_path)
mkpath(plot_path)

# EXECUTE MODEL

# Run Gillespie SSA model with above parameters.
# NOTE Static seed will be assigned from this file.
include("$(Base.source_dir())/gillespie_model.jl")

# DEFINE MUTATION TIME-LINE PLOT (PLOT 1)

years = times / year
counts = [wild_copy_median, mutant_copy_median]  # n

print("\nCreating mutation/time plot $(yellow)plots/$(basename(plot_1_path))$(reset)... ")
fig = plot(
    years,
    counts,
    xlabel="Time (years)",
    ylabel="mtDNA count (n)",
    ylims=(-5, (agent_max + 10)),
    title="Patient with random mutant mtDNA inheritance",
    label=["Wild count" "Mutant count"],
    legend=:outertopright,
    dpi=1200
)

savefig(fig, "$(plot_1_path)")
println("$(green)Done$(reset)")

# DEFINE PERCENTILE PLOT (PLOT 2)

quantiles = [upper_quantile, middle_quantile, lower_quantile] * agent_max

print("Creating quantile plot $(yellow)plots/$(basename(plot_2_path))$(reset)... ")
fig2 = plot(
    years,
    quantiles,
    xlabel="Time (years)",
    ylabel="mtDNA count (n)",
    ylims=(-5, (agent_max + 10)),
    title="Patient with random mutant mtDNA inheritance",
    label=["95th percentile" "50th percentile" "5th percentile"],
    legend=:outertopright,
    dpi=1200
)

savefig(fig2, "$(plot_2_path)")
println("$(green)Done$(reset)")

# DEFINE CSV EXPORT FILE

data = DataFrame(
    steps = time_list,
    wild_median = wild_copy_median,
    mutant_median = mutant_copy_median,
    total_counts = total_copy_median
)

print("Writing data to $(yellow)data/$(basename(csv_path))$(reset) Please wait... ")
CSV.write("$(csv_path)", data)
println("$(green)Done$(reset)")

# End of File.
