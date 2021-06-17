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

# - Requires Julia >= v1.6.

###########
# Imports #
###########

using Random

#############
# Variables #
#############

# Define elements for random simulaiton.

max_elements = 200

random_wild = rand(1:max_elements)
random_mutant = max_elements - random_wild

# Initial concentrations of wild-type and mutant mtDNA.
u0 = [random_wild, random_mutant]

# 80 Years = 960 Months.
tend = 960.0

# Kinetic rates of reactions.
# parameters = (r=0.01, m=0.001, d=0.01)

# 1 Month / 28 Days / 24 Hours = Hourly Rate
# 1 / 28 / 24 = 0.001488095
parameters = (r=0.001488095, m=0.001488095, d=0.001488095)

# Number of simulation repeats.
loops = 1000

#############
# Kickstart #
#############

# EXECUTE MODEL

# Run Gillespie SSA model with above parameters.
include("$(pwd())/gillespie_model.jl")

# Ensure results are random.
Random.seed!()


# DEFINE PLOT

# Define plot axis elements.
x = times / 12.0  # Convert months to years.
y = [mean_wild, mean_mutant, mean_trend] * 100

print("\nWriting plot to '$(pwd())/random_simulation.png'... ")
fig = plot(
    x,  # Temporal States
    y,  # Molecule Concentration [Wild-type, Mutant]
    xlabel="Time (Years)",
    ylabel="Molecule Concentration (%)",
    ylims=(0, 100),
    title="Random Mutant mtDNA Inheritance (SSA Model)",
    label=["Wild-type" "Mutant" "Trend"],
    dpi=1200
)

# Save plot in current working directory.
savefig(fig, "$(pwd())/random_simulation.png")
println("Done")


# DEFINE QUANTILE PLOT

# Define axis elements.
y2 = [lower_quantile, upper_quantile, quantile_trend] * 100

print("Writing plot to '$(pwd())/random_simulation_quantile.png'... ")
fig2 = plot(
    x,   # Temporal States
    y2,  # Certainty [Lower Quantile, Upper Quantile]
    xlabel="Time (Years)",
    ylabel="Certainty (%)",
    ylims=(0, 100),
    title="Random Mutant mtDNA Inheritance (Quantiles)",
    label=["Lower Quantile" "Upper Quantile" "Trend"],
    dpi=1200
)

# Save plot in current working directory.
savefig(fig2, "$(pwd())/random_simulation_quantile.png")
println("Done")

# End of File.
