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
# - Time-scale for this model is that integer '1' equals 1 month.
#   - Scales are adjusted to years when defining the plot only.
#   - With 'delta' being '1', we are recording results monthly.

###########
# Imports #
###########

using Plots,
      StatsPlots

#############
# Variables #
#############

# Initial concentrations of wild-type and mutant mtDNA.
u0 = [175, 25]

# 80 Years = 960 Months.
tend = 960.0

# Number of simulation repeats.
loops = 1000

# VARIABLES > KINETIC RATES

# Original rates: tuned to favour longer simulation runtime.
# parameters = (r=0.01, m=0.001, d=0.01)

# 1 Month / 28 Days = Daily Rate (0.035714286)
# parameters = (r=0.035714286, m=0.035714286, d=0.035714286)

# 1 Month / 28 Days / 24 Hours = Hourly Rate (0.001488095)
parameters = (r=0.001488095, m=0.001488095, d=0.001488095)

#############
# Kickstart #
#############

# EXECUTE MODEL

# Run Gillespie SSA model with above parameters.
include("$(pwd())/gillespie_model.jl")


# DEFINE MUTATION TIME-LINE PLOT

# Define plot axis elements.
x = times / 12.0       # Convert months to years.
y = mean_mutant * 100  # Convert mutation level to percentage.

print("\nCreating mutation time-line plot 'patient_with_inheritance.png'... ")
fig = plot(
    x,  # Temporal States
    y,  # Mutation Level
    xlabel="Time (years)",
    ylabel="Mutation level (%)",
    ylims=(0, 100),
    title="Patient with mutant mtDNA inheritance",
    legend=false,
    dpi=1200
)

# Save plot in current working directory.
savefig(fig, "$(pwd())/patient_with_inheritance.png")
println("Done")


# DEFINE QUANTILE PLOT

# Define axis elements.
y2 = [upper_quantile, lower_quantile] * 100

print("Creating quantile plot 'patient_with_inheritance_quantile.png'... ")
fig2 = plot(
    x,   # Temporal States
    y2,  # Certainty [Lower Quantile, Upper Quantile]
    xlabel="Time (Years)",
    ylabel="Certainty (%)",
    ylims=(0, 100),
    title="Positive Mutant mtDNA Inheritance (Quantiles)",
    label=["97.5 Percentile" "2.5 Percentile"],
    dpi=1200
)

# Save plot in current working directory.
savefig(fig2, "$(pwd())/patient_with_inheritance_quantile.png")
println("Done")

# DEFINE DENSITY PLOT (PLOT 3)

# Define axis elements.
x3 = mean_mutant

print("Creating density plot 'patient_with_inheritance_density.png'... ")
fig3 = density(
    vec(x3),  # Mean of mutant levels
    title="Density of mutation mean",
    legend=false,
    dpi=1200
)

savefig(fig3, "$(pwd())/patient_with_inheritance_density.png")
println("Done")

# End of File.
