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

import Pkg

Pkg.add("Distributions")
Pkg.add("Plots")
Pkg.add("StatsPlots")

using Distributions
using Plots
using StatsPlots

#############
# Variables #
#############

# VARIABLES > INITIAL CONDITIONS

# Initial concentrations of wild-type and mutant mtDNA.
u0 = [200, 0]

# VARIABLES > PATHS

plot_path = "$(pwd())/plots/patient_without_inheritance"

plot_1_path = "$(plot_path)/01_timeline.png"
plot_2_path = "$(plot_path)/02_quantiles.png"
plot_3_path = "$(plot_path)/03_density.png"

#############
# Kickstart #
#############

# Create plot directory.
mkpath(plot_path)

# EXECUTE MODEL

# Run Gillespie SSA model with above parameters.
# NOTE Static seed will be assigned from this file.
include("$(pwd())/gillespie_model.jl")

# DEFINE MUTATION TIME-LINE PLOT (PLOT 1)

# Define plot axis elements.
x = times / 12.0         # Convert months to years.
y = mean_mutant * 100.0  # Convert mutation level to percentage.

print("\nCreating mutation/time plot '$(plot_1_path)'... ")
fig = plot(
    x,  # Temporal States
    y,  # Mutation Level
    xlabel="Time (years)",
    ylabel="Mutation level (%)",
    ylims=(-10, 100),  # Comment to automatically zoom.
    title="Patient without mutant mtDNA inheritance",
    legend=false,
    dpi=1200
)

# Save plot in current working directory.
savefig(fig, "$(plot_1_path)")
println("Done")

# DEFINE PERCENTILE PLOT (PLOT 2)

# Define axis elements.
y2 = [upper_quantile, middle_quantile, lower_quantile] * 100

print("Creating quantile plot '$(plot_2_path)'... ")
fig2 = plot(
    x,   # Temporal States
    y2,  # Mutation level [2.5th percentile, 50th percentile, 97.5th percentile]
    xlabel="Time (years)",
    ylabel="Mutation level (%)",
    ylims=(-10, 100),  # Comment to automatically zoom.
    title="Patient without mutant mtDNA inheritance",
    label=["97.5th percentile" "50th percentile" "2.5th percentile"],
    linealpha=0.5,
    dpi=1200
)

# Save plot in current working directory.
savefig(fig2, "$(plot_2_path)")
println("Done")

# DEFINE DENSITY PLOT (PLOT 3)

# Define axis elements.
x3 = vec(mean_mutant)

print("Creating density plot '$(plot_3_path)'... ")
fig3 = density(
    x3,  # Mean of mutant levels
    title="Patient without mutant mtDNA inheritance",
    xlabel="Density (mutation mean)",
    legend=false,
    dpi=1200
)

savefig(fig3, "$(plot_3_path)")
println("Done")

# End of File.
