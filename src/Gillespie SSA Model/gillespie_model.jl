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

##############
# References #
##############

# - <https://nextjournal.com/bebi5009/gillespie-julia>

###########
# Imports #
###########

# TODO Remove unused libraries.

import Pkg

Pkg.add([
    "IterTools",
    "Plots",
    "StatsBase",
    "StatsPlots"
])

using IterTools,
      Plots,
      Random,
      Statistics,
      StatsBase,
      StatsPlots

#################
# Prerequisites #
#################

# Ensure results are reproducible.
Random.seed!(41269)

#############
# Variables #
#############

# - Variables are to be defined within a simulation file.
#   - e.g. `patient_with_inheritance.jl`.

#############
# Functions #
#############


# FUNCTIONS > GILLESPIE MODEL

"""
`ssa(model, u0, tend, p, choose_stoich, tstart; delta=1.0)`

Adapted from: Chemical and Biomedical Enginnering Calculations Using
Python Ch.4-3
"""
function ssa(model, u0, tend, p, choose_stoich, tstart=zero(tend); delta=1.0)
    t = tstart    # Current time.
    ts = [t]      # List of reaction times.
    u = copy(u0)  # Current state.
    us = copy(u)  # Record of states.

    times = [tstart: delta: tend;]
    tindex = 2

    while t < tend && u[1] + u[2] > 0
        if t >= times[tindex]
            us = [us u]
            push!(ts, t)
            tindex = tindex + 1
        end

        dx = model(u, p, t)
        total_hazard = sum(dx)
        dt = Random.randexp() / total_hazard
        stoich = choose_stoich(dx, total_hazard)
        u .+= stoich
        t += dt
    end

    us = collect(us')

    return (t = ts, u = us)
end


"""
`model(u, p, t)`

Propensity model for this reaction.
Reaction of `Wild-type <-> Mutant` with rate constants `u[1]` & `u[2]`.
"""
function model(u, p, t)
    # Define copy number.
    cn = u[1] + u[2]

    # Prevent early extinction.
    if cn < 50
        r1 = p.r * u[1]
        r4 = p.r * u[2]
    else
        r1 = 0
        r4 = 0
    end

    return [r1,          # Reaction 1: Replicate wild-type mtDNA.
            p.d * u[1],  # Reaction 2: Degrade wild-type mtDNA.
            p.m * u[1],  # Reaction 3: Mutate wild-type mtDNA.
            r4,          # Reaction 4: Replicate mutant mtDNA.
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
    else
        stoich = [0, -1]
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

    # Start at loop number 1 in intervals of 1 until 'n'.
    for i in 1:1:n
        # Define and run the simulation.
        sol = ssa(model, u0, tend, parameters, choose_stoich)

        # - Take the results and store them in 'time_states' and
        #   'concentration_states' respectively.
        push!(time_states, sol.t)
        push!(concentration_states, sol.u)
    end

    return(t = time_states, u = concentration_states)
end


# FUNCTIONS > NaN HANDLING

nanmean(x) = mean(filter(!isnan, x))
nanmean(x, y) = mapslices(nanmean, x, dims=y)

nanmedian(x) = median(filter(!isnan, x))
nanmedian(x, y) = mapslices(nanmedian, x, dims=y)

nanquantile(x, y) = quantile(filter(!isnan,x), y)

###########
# Exports #
###########

# Make functions available to Julia REPL.

# FUNCTIONS > GILLESPIE MODEL

export ssa
export model
export choose_stoich
export loop_simulation

# FUNCTIONS > NaN HANDLING

export nanmean
export nanmedian
export nanquantile

#############
# Kickstart #
#############

# Run simulation 'loops' time(s) and assign output to 'results'.
results = loop_simulation(loops)

# Calculate total number of elements.
num_simulations = size(results.u)[1]
num_times = size(results.u[1])[1]
num_species = size(results.u[1])[2]

# Create a list of times.
times = results.t[1]

# Preallocate array for states.
molecules = fill(NaN, num_simulations, num_times, num_species)

# Populate 'molecules' array.
# NOTE Empty spaces will contain NaN.
for i in 1:num_simulations
    for j in 1:num_times
        for k in 1:num_species
            if size(results.u[i])[1] >= j
                molecules[i, j, k] = results.u[i][j, k]
            end
        end
    end
end

# Calculate total number of molecules.
total = sum(molecules, dims=3)

# Calculate percentage of wild-type and mutant mtDNA.
wild_load = molecules[:, :, 1] ./total
mutation_load = molecules[:, :, 2] ./ total

# Calculate metrics for plotting.
mean_wild = nanmean(wild_load, 1)
mean_mutant = nanmean(mutation_load, 1)
median_wild = nanmedian(wild_load, 1)
median_mutant = nanmedian(mutation_load, 1)

# Flatten objects to simple arrays.
mean_wild = collect(Iterators.flatten(mean_wild))
mean_mutant = collect(Iterators.flatten(mean_mutant))
median_wild = collect(Iterators.flatten(median_wild))
median_mutant = collect(Iterators.flatten(median_mutant))

# Trends for mean and median concentrations.  TODO Remove this.
# - Ensure that the smallest is subtracted from the largest.

if mean_wild > mean_mutant
    mean_trend = mean_wild - mean_mutant
else
    mean_trend = mean_mutant - mean_wild
end

if median_wild > median_mutant
    median_trend = median_wild - median_mutant
else
    median_trend = median_mutant - median_wild
end


# DEFINE QUANTILES

upper_quantile = fill(NaN, num_times)
middle_quantile = fill(NaN, num_times)
lower_quantile = fill(NaN, num_times)

for j in 1:num_times
    mutation_loads = mutation_load[:, j]
    upper_quantile[j] = nanquantile(mutation_loads, 0.975)
    middle_quantile[j] = nanquantile(mutation_loads, 0.5)
    lower_quantile[j] = nanquantile(mutation_loads, 0.025)
end

# Trend of quantiles. TODO Remove this.

if upper_quantile > lower_quantile
    quantile_trend = upper_quantile - lower_quantile
else
    quantile_trend = lower_quantile - upper_quantile
end


# DEFINE PLOT

# - Plot and elements are to be defined within a simulation file.
#   - e.g. `patient_with_inheritance.jl`.


# TEXT REPORT

# Generate specific time-frame elements.
runtime_in_years = length(mean_mutant) / 12.0

int_years = Int(floor(runtime_in_years))

int_months = runtime_in_years - floor(runtime_in_years)
int_months = Int(floor(int_months * 12))

# Print final report.
println("\nRESULTS\n=======\n")
println("Simulation was looped $(loops) time(s).")
println("Initial concentration of wild-type mtDNA: $(u0[1]) molecule(s).")
println("Initial concentration of mutant mtDNA: $(u0[2]) molecule(s).")
println("Population went extinct at $(num_times)/$(Int(floor(tend))) cycle(s) ",
    "or $(int_years) year(s) and $(int_months) month(s).")
println("97.5th percentile (Head): $(upper_quantile[1:5])")
println("50th percentile (Head):   $(middle_quantile[1:5])")
println("2.5th percentile (Head):  $(lower_quantile[1:5])")

# End of File.
