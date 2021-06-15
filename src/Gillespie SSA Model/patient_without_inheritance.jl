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

#############
# Variables #
#############

# Initial concentrations of wild-type and mutant mtDNA.
u0 = [200, 0]

# 80 Years = 960 Months.
tend = 960.0

# Kinetic rates of reactions.
parameters = (r=0.01, m=0.001, d=0.01)

# Number of simulation repeats.
loops = 1000

#############
# Kickstart #
#############

# EXECUTE MODEL

# Run Gillespie SSA model with above parameters.
include("$(pwd())/gillespie_model.jl")


# DEFINE PLOT

# Define plot axis elements.
x = times / 12.0  # Convert months to years.
y = [mean_wild, mean_mutant]

print("\nWriting plot to '$(pwd())/patient_without_inheritance'... ")
fig = plot(
    x,  # Temporal States
    y,  # Molecule Concentration [Wild-type, Mutant]
    xlabel="Time (Years)",
    ylabel="Number of Molecules (%)",
    ylims=(0, 1),
    title="mtDNA Population Dynamics (SSA Model)",
    label=["Wild-type" "Mutant"],
    dpi=1200
)

# Save plot in current working directory.
savefig(fig, "$(pwd())/patient_without_inheritance.png")
println("Done")

# End of File.
