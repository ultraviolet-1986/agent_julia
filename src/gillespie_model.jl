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

# - Code written to detect the time and concentration arrays will not
#   function correctly for Julia versions older than 1.6.x.

##############
# References #
##############

# <https://nextjournal.com/bebi5009/gillespie-julia>

###########
# Imports #
###########

import Pkg

Pkg.add(["IterTools",
         "Plots",
         "StatsBase"])

using IterTools,
      Plots,
      Random,
      Statistics,
      StatsBase

#################
# Prerequisites #
#################

# Ensure results are reproducible.
Random.seed!(41269)

#############
# Variables #
#############

# Initial concentrations of species 'A' and 'B'.
u0 = [175, 25]

# Time at which the simulation will stop.
tend = 60.0  # tend > 27 will crash.

# Kinetic rates of reactions.
parameters = (r=1.0, m=0.1, d=1.0)

#############
# Functions #
#############


"""
`ssa(model, u0, tend, p, choose_stoich, tstart)`

Adapted from: Chemical and Biomedical Enginnering Calculations Using
Python Ch.4-3
"""
function ssa(model, u0, tend, p, choose_stoich, tstart=zero(tend); delta=0.1)
    t = tstart    # Current time.
    ts = [t]      # List of reaction times.
    u = copy(u0)  # Current state.
    us = copy(u)  # Record of states.

    times = [tstart: delta: tend;] # Sequence from 0.0 - 10.0
    tindex = 2  # Initial step is already defined.

    while t < tend
        # If time > next sample, do this. Update sample to be +1 week.
        if t >= times[tindex]
            us = [us u]
            push!(ts, t)  # Record t
            tindex = tindex + 1
        end

        dx = model(u, p, t)  # Propensity of reactions.
        total_hazard = sum(dx)
        dt = Random.randexp() / total_hazard  # Time step.
        stoich = choose_stoich(dx, total_hazard)  # Get stoichiometry.
        u .+= stoich
        t += dt
    end

    us = collect(us')

    return (t = ts, u = us)
end


"""
`model(u, p, t)`

Propensity model for this reaction.
Reaction of `A <-> B` with rate constants `k1` & `k2`.
"""
function model(u, p, t)
    return [p.r * u[1],  # Reaction 1: Replicate wild-type mtDNA.
            p.d * u[1],  # Reaction 2: Degrade wild-type mtDNA.
            p.m * u[1],  # Reaction 3: Mutate wild-type mtDNA.
            p.r * u[2],  # Reaction 4: Replicate mutant mtDNA.
            p.d * u[2]]  # Reaction 5: Degrade mutant mtDNA.
end


@doc raw"""
`choose_stoich(dx, dxsum)`

This function will choose and return the Stoichiometry which dictates
what reaction has occurred.

# Reactions

1. Replicate wild-type mtDNA.
```math
M_1 \begin{array}{c}rM_1\\ \Rightarrow\end{array} 2M_1
```

2. Degrade wild-type mtDNA.
```math
M_1 \begin{array}{c}dM_1\\ \Rightarrow\end{array} \emptyset
```

3. Mutate wild-type mtDNA.
```math
M_1 \begin{array}{c}mM_1\\ \Rightarrow\end{array} M_1 + M_2
```

4. Replicate mutant mtDNA.
```math
M_2 \begin{array}{c}rM_2\\ \Rightarrow\end{array} 2M_2
```

5. Degrade mutant mtDNA.
```math
M_2 \begin{array}{c}dM_2\\ \Rightarrow\end{array} \emptyset
```
"""
function choose_stoich(dx, dxsum = sum(dx))
    # Calculate probability.
    sections = cumsum(dx / dxsum)

    # Select pseudo-random number.
    # NOTE This can be influenced by seed.
    roll = rand()

    # Reaction 1: Replicate wild-type mtDNA.
    if roll <= sections[1]
        stoich = [1, 0]

    # Reaction 2: Degrade wild-type mtDNA.
    elseif roll <= sections[2]
        stoich = [-1, 0]

    # Reaction 3: Mutate wild-type mtDNA.
    elseif roll <= sections[3]
        stoich = [0, 1]

    # Reaction 4: Replicate mutant mtDNA.
    elseif roll <= sections[4]
        stoich = [0, 1]

    # Reaction 5: Degrade mutant mtDNA.
    elseif roll < 1.0
        stoich = [0, -1]

    # Error: catch all.
    else
        error("Stoichiometry is out of bounds.")
    end

    return stoich
end


"""
`loop_simulation(n::Int64)`

Perform the simulation `n` number of times and return the result. Note
that this function is hard-coded to run a single model defined eslwhere
in the `gillespie_model.jl` file.

# Examples

```julia-repl
julia> loop_simulation(10)
```

"""
function loop_simulation(n::Int64)
    # Create new arrays for holding all results.
    time_states = []
    concentration_states = []

    println("Simulation will be performed $(n) time(s).")

    # Start at loop number 1 in intervals of 1 until 'n'.
    for i in 1:1:n
        # Define and run the simulation.
        sol = ssa(model, u0, tend, parameters, choose_stoich)

        # - Take the results and store them in 'time_states' and
        #   'concentration_states' respectively.
        push!(time_states, sol.t)
        push!(concentration_states, sol.u)
    end

    return(t = time_states, c = concentration_states)
end


#############
# Kickstart #
#############

# Run simulation 'n' time(s) and assign returns to 'results'.
results = loop_simulation(1000)

# Calculate total number of elements.
num_simulations = size(results.c)[1]
num_times = size(results.c[1])[1]
num_species = size(results.c[1])[2]

# Create a list of times.
times = results.t[1]

# Preallocate array for states.
molecules = Array{Float64}(undef, num_simulations, num_times, num_species)

for i in 1:num_simulations  # Error if tend > 27.
    for j in 1:num_times
        for k in 1:num_species
            molecules[i, j, k] = results.c[i][j, k]  # Error if tend > 27.
        end
    end
end

# Calculate total number of molecules.
total = sum(molecules, dims=3)

# Calculate percentage of mutant mtDNA.
mutation_load = molecules[:, :, 2] ./ total

mean_mutant = mean(mutation_load, dims=1)
median_mutant = median(mutation_load, dims=1)

# EXPERIMENTAL Flatten objects to simple arrays.
mean_mutant = collect(Iterators.flatten(mean_mutant))
median_mutant = collect(Iterators.flatten(median_mutant))

# Calculate percentage of wild-type mtDNA.
wild_load = molecules[:, :, 1] ./total

mean_wild = mean(wild_load, dims=1)
median_wild = median(wild_load, dims=1)

# EXPERIMENTAL Flatten objects to simple arrays.
mean_wild = collect(Iterators.flatten(mean_wild))
median_wild = collect(Iterators.flatten(median_wild))

upper_quantile = Array{Float64}(undef, num_times)
lower_quantile = Array{Float64}(undef, num_times)

for j in 1:num_times
    mutation_loads = mutation_load[:, j]
    upper_quantile[j] = quantile(mutation_loads, 0.975)
    lower_quantile[j] = quantile(mutation_loads, 0.025)
end

# Create plot axis elements.
x = times
# y = [mean_wild, mean_mutant]
y = [median_wild, median_mutant]

# Create plot of median values.
fig = plot(
    x,  # Temporal States
    y,  # Molecule Concentration
    xlabel="Time",
    ylabel="Number of Molecules (%)",
    # xlims=(0, 10), # Between 0 and 10 for time states.
    ylims=(0, 1),  # Between 0 and 1 for percentage.
    title="mtDNA Population Dynamics (SSA Model)",
    label=["Wild-type" "Mutant"],
    dpi=300)

savefig(fig, "plot_2.png")

# End of File.
