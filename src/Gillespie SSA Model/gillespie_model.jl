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

# - <https://nextjournal.com/bebi5009/gillespie-julia>

###########
# Imports #
###########

using Distributions
using IterTools
using Random
using Statistics

#################
# Prerequisites #
#################

# Ensure results are reproducible.
Random.seed!(41269)

#############
# Variables #
#############

# COLOURS > TEXT OUTPUT

green = "\e[32m"
red = "\e[31m"
yellow = "\e[33m"
reset = "\e[0m"

# VARIABLES > INITIAL CONDITIONS

# The `u0` variable is to be defined within a simulation file e.g.
# `patient_with_inheritance.jl`. Terminate execution if variable does
# not exist.

if (! @isdefined u0)
    error("Please run a simulation file e.g.: 'patient_with_inheritance.jl'.")
end

agent_max = u0[1] + u0[2]

# Number of simulation repeats.
loops = 1000

# VARIABLES > TEMPORAL UNITS

α = Float64(1.0e6)

# Apply α to reduce data-size.
day = Float64(24.0 * 3600.0 / α)

week = Float64(day * 7.0)
year = Float64(day * 365.0)
month = Float64(year / 12.0)
hour = Float64(day / 24.0)

δ = month

# Target end time.
tend = Float64(year * 80.0)

# VARIABLES > KINETIC RATES

λ = rand(Normal(260.0, 1), 1)[1]

reaction_rate = log(2) / λ
reaction_rate = reaction_rate ./ (day)

# mutation_rate = 1.157e-12  # Mutation enabled.
mutation_rate = 0.0        # Mutation disabled.

parameters = (r=reaction_rate, m=mutation_rate, d=reaction_rate)

#############
# Functions #
#############

# FUNCTIONS > GILLESPIE MODEL

"""
`ssa(model, u0, tend, p, choose_stoich, tstart=zero(tend); δ=month::Float64)`

Adapted from: Chemical and Biomedical Enginnering Calculations Using
Python Ch.4-3
"""
function ssa(model, u0, tend, p, choose_stoich, tstart=zero(tend); δ=month::Float64)
    t = tstart    # Current time.
    ts = [t]      # List of reaction times.
    u = copy(u0)  # Current state.
    us = copy(u)  # Record of states.

    times = [tstart: δ: tend;]
    tindex = 2

    while t < tend
        if t ≥ times[tindex]
            us = [us u]
            push!(ts, t)
            tindex = tindex + 1
        end

        dx = model(u, p)
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
`model(u, p; target_upper=250::Int64, target_lower=150::Int64)`

Propensity model for this reaction.
Reactions of `wild-type mtDNA` and `mutant mtDNA` using kinetic rates
`p.d`, `p.m` and `p.r`.
"""
function model(u, p; target_upper=250::Int64, target_lower=150::Int64)
    cn = u[1] + u[2]  # Define copy number.

    rxn1 = p.r * u[1]  # Reaction 1: Replicate wild-type mtDNA.
    rxn2 = p.d * u[1]  # Reaction 2: Degrade wild-type mtDNA.
    rxn3 = p.m * u[1]  # Reaction 3: Mutate wild-type mtDNA.
    rxn4 = p.r * u[2]  # Reaction 4: Replicate mutant mtDNA.
    rxn5 = p.d * u[2]  # Reaction 5: Degrade mutant mtDNA.

    # Disable replication to prevent population explosion.
    if cn ≥ target_upper
        rxn1 = 0.0
        rxn4 = 0.0

    # Disable degradation to prevnet population extinction.
    elseif cn ≤ target_lower
        rxn2 = 0.0
        rxn5 = 0.0
    end

    return [rxn1, rxn2, rxn3, rxn4, rxn5]
end


@doc raw"""
`choose_stoich(dx, dxsum=sum(dx))`

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
function choose_stoich(dx, dxsum=sum(dx))
    # Calculate probability.
    sections = cumsum(dx / dxsum)

    # Select pseudo-random number.
    # NOTE This can be influenced by seed.
    roll = rand()

    # Reaction 1: Replicate wild-type mtDNA.
    if roll ≤ sections[1]
        stoich = [1, 0]

    # Reaction 2: Degrade wild-type mtDNA.
    elseif roll ≤ sections[2]
        stoich = [-1, 0]

    # Reaction 3: Mutate wild-type mtDNA.
    elseif roll ≤ sections[3]
        stoich = [0, 1]

    # Reaction 4: Replicate mutant mtDNA.
    elseif roll ≤ sections[4]
        stoich = [0, 1]

    # Reaction 5: Degrade mutant mtDNA.
    else
        stoich = [0, -1]
    end

    return stoich
end


"""
`loop_simulation(n=1::Int64)`

Perform the simulation `n` number of times and return the result. Note
that this function is hard-coded to run a single model defined eslwhere
in the `gillespie_model.jl` file and will run the simulation only once
if an argument is not provided.

# Examples

```julia-repl
julia> loop_simulation(10)
```
"""
function loop_simulation(n=1::Int64)
    # Create new arrays for holding all results.
    time_states = []
    concentration_states = []

    # Start at loop number 1 in intervals of 1 until `n`.
    for i in 1:1:n
        # Define and run the simulation.
        sol = nothing
        sol = ssa(model, u0, tend, parameters, choose_stoich)

        # Take the results and store them in `time_states` and
        # `concentration_states` respectively.
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

#############
# Kickstart #
#############

# Run simulation `loops` time(s) and assign output to `results`.
results = nothing
results = loop_simulation(loops)

# Calculate total number of elements.
num_simulations = size(results.u)[1]
num_times = size(results.u[1])[1]
num_species = size(results.u[1])[2]

# Create a list of times.
# NOTE Each row will be exactly the same.
times = results.t[1]

# Preallocate array for states.
molecules = fill(NaN, num_simulations, num_times, num_species)

# Populate `molecules` array.
# NOTE Empty spaces will contain NaN.
for i in 1:num_simulations
    for j in 1:num_times
        for k in 1:num_species
            if size(results.u[i])[1] ≥ j
                molecules[i, j, k] = results.u[i][j, k]
            end
        end
    end
end

# Calculate total number of molecules.
total = sum(molecules, dims=3)

# Calculate percentage of wild-type and mutant mtDNA.
wild_load = molecules[:, :, 1] ./ total
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

# Create equivalent count numbers based on mean and median.
wild_copy_mean = mean_wild .* agent_max
mutant_copy_mean = mean_mutant .* agent_max
total_copy_mean = (wild_copy_mean .+ mutant_copy_mean)

wild_copy_median = (median_wild .* agent_max)
mutant_copy_median = (median_mutant .* agent_max)
total_copy_median = (wild_copy_median .+ mutant_copy_median)

# Define time list

time_steps = tend / δ
time_list = Int.([1:1:time_steps;])

# DEFINE QUANTILES

upper_quantile = fill(NaN, num_times)
middle_quantile = fill(NaN, num_times)
lower_quantile = fill(NaN, num_times)

for j in 1:num_times
    mutation_loads = mutation_load[:, j]
    upper_quantile[j] = nanquantile(mutation_loads, 0.95)  # 95%
    middle_quantile[j] = nanquantile(mutation_loads, 0.5)  # 50%
    lower_quantile[j] = nanquantile(mutation_loads, 0.05)  #  5%
end

# TEXT REPORT

# Print final report.

println("\nRUNNING GILLESPIE SSA SIMULATION")
println("================================\n")

println("Simulation was looped $(yellow)$(loops)$(reset) time(s).")
println("Initial concentration of wild-type mtDNA: $(yellow)$(u0[1])$(reset) molecule(s).")
println("Initial concentration of mutant mtDNA: $(yellow)$(u0[2])$(reset) molecule(s).")
println("95th percentile (Head): $(yellow)$(upper_quantile[1:5])$(reset)")
println("50th percentile (Head): $(yellow)$(middle_quantile[1:5])$(reset)")
println(" 5th percentile (Head): $(yellow)$(lower_quantile[1:5])$(reset)")

# End of File.
