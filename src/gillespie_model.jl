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

###########
# Imports #
###########

import Pkg

Pkg.add("Plots")

using Plots,
      Random

#################
# Prerequisites #
#################

# Ensure results are reproducible.
Random.seed!(41269)

#############
# Variables #
#############

# Initial concentrations of species 'A' and 'B'.
u0 = [200, 50]

# Time at which the simulation will stop.
tend = 10.0

# Kinetic rates of reactions.
parameters = (r=1.0, m=0.0, d=1.0)

#############
# Functions #
#############

"""
ssa(model, u0, tend, p, choose_stoich, tstart)

Adapted from: Chemical and Biomedical Enginnering Calculations Using
Python Ch.4-3
"""
function ssa(model, u0, tend, p, choose_stoich, tstart=zero(tend); delta=0.1)
    t = tstart    # Current time.
    ts = [t]      # List of reaction times.
    u = copy(u0)  # Current state.
    us = copy(u)  # Record of states.

    times = [tstart: delta: tend;] # Sequence from 0.0 - 10.0
    tindex = 2

    while t < tend
        dx = model(u, p, t)              # propensities of reactions.
        total_hazard = sum(dx)
        dt = Random.randexp() / total_hazard  # Time step.
        stoich = choose_stoich(dx, total_hazard)       # Get stoichiometry.
        u .+= stoich
        t += dt

        # If time > next sample, do this. Update sample to be +1 week.
        # Add to record
        if t >= times[tindex]
            us = [us u]
            push!(ts, t)  # Record t
            tindex = tindex + 1
        end
    end

    us = collect(us')

    return (t = ts, u = us)
end


"""
model(u, p, t)

Propensity model for this reaction.
Reaction of A <-> B with rate constants k1 & k2.
"""
function model(u, p, t)
    return [p.r * u[1], p.d * u[1], p.m * u[1], p.r * u[2], p.d * u[2]]
end


"""
choose_stoich(dx, dxsum)

Choose and return which Stoichiometry to update the state.
"""
function choose_stoich(dx, dxsum = sum(dx))
    sections = cumsum(dx / dxsum)  # Create probability.
    roll = rand()

    #=
    Reaction 1 [1, 0], rM_1
    Reaction 2 [-1, 0], dM_1
    Reaction 3 [0, 1], mM_1
    Reaction 4 [0, 1], rM_2
    Reaction 5 [0, -1], dM_2
    =#

    if roll <= sections[1]
        stoich = [1, 0]
    elseif roll <= sections[2]
        stoich = [-1, 0]
    elseif roll <= sections[3]
        stoich = [0, 1]
    elseif roll <= sections[4]
        stoich = [0, 1]
    elseif roll < 1.0
        stoich = [0, -1]
    else
        println("Whoops!")
    end

    return stoich
end

#############
# Kickstart #
#############

# Perform the simulation and assign results to 'sol'.
print("\nNow performing simulation... ")
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
    label=["Wild-type" "Mutant"],
    dpi=300)

# Output plot to 'plot.png'.
savefig(fig, "plot.png")
println("Done!")

# End of File.
