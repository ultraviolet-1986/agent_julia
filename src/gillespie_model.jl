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

# <https://nextjournal.com/bebi5009/gillespie-julia>

#################
# Prerequisites #
#################

import Pkg

Pkg.add(["CSV",
         "DataFrames",
         "Plots"])

using CSV,
      DataFrames,
      Random,
      Plots

Random.seed!(41269)

#############
# Variables #
#############

# Initial concentrations of species 'A' and 'B' respectively.
u0 = [200, 0]

# Time at which the simulation will stop.
tend = 10.0

# Kinetic rates of reactions.
parameters = (k1=1.0, k2=0.5)

#############
# Functions #
#############

"""
ssa(model, u0, tend, p, choose_stoich, tstart)

Adapted from: Chemical and Biomedical Enginnering Calculations Using
Python Ch.4-3
"""
function ssa(model, u0, tend, p, choose_stoich, tstart=zero(tend))
    t = tstart    # Current time.
    ts = [t]      # List of reaction times.
    u = copy(u0)  # Current state.
    us = copy(u)  # Record of states.

    # Cast 'ts' and 'us' as DataFrames.
    # us = DataFrame(us)
    # ts = DataFrame(ts)

    while t < tend
        dx = model(u, p, t)              # propensities of reactions.
        dt = Random.randexp() / sum(dx)  # Time step.
        stoich = choose_stoich(dx)       # Get stoichiometry.
        u .+= stoich
        t += dt

        # If time > next sample, do this. Update sample to be +1 week.
        # Add to record
        us = [us u]              ## TO COMMENT FOR DATAFRAME METHOD
        push!(ts, t)  # Record t ## TO COMMENT FOR DATAFRAME METHOD
    end

    us = collect(us') ## TO COMMENT FOR DATAFRAME METHOD

    return (t = ts, u = us)  ## TO CHANGE FOR DATAFRAME METHOD
end


"""
model(u, p, t)

Propensity model for this reaction.
Reaction of A <-> B with rate constants k1 & k2.
"""
function model(u, p, t)
    [p.k1 * u[1],  p.k2 * u[2]]
end


"""
choose_stoich(dx, dxsum)

Choose and return which Stoichiometry to update the state.
"""
function choose_stoich(dx, dxsum = sum(dx))
    sections = cumsum(dx ./ dxsum)
    roll = rand()

    if roll <= sections[1]
        stoich = [-1, 1]
    else
        stoich = [1, -1]
    end

    return stoich
end

#############
# Kickstart #
#############

# Perform the simulation and assign results to 'sol'.
print("\nNow running sumulation... ")
sol = ssa(model, u0, tend, parameters, choose_stoich)
println("Done!")

# Define the plot.
print("Now creating the plot... ")
fig = plot(
    sol.t,
    sol.u,
    xlabel="time",
    ylabel="# of molecules",
    title = "SSA",
    label=["A" "B"],
    dpi=300)

# Output plot to 'plot.png'.
savefig(fig, "plot.png")
println("Done!")

# BUILD DATAFRAME FROM COMPLETED DATA
# - This will help me to create a cohesive table of the completed data.
# - From here, it should be possible to ensure that changes to function
#   code will work correctly when the time comes.
# - The desired output of this model is a time-index and values such as
#   [0, 1, -1] for reaction rows for each species per-column.

# Print each Array on separate lines.
println("\n$(sol.t)")
println("\n$(sol.u)")

# CURRENTLY BROKEN, REQUIRES INDEX AND COLUMN NAMES
# data_frame = DataFrame(sol.t, sol.u)
# println(data_frame)

# End of File.
