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

    # tstart / Start of simulation.
    # delta  / Size of step to take.
    # tend   / Time at which simulation will be terminated.
    times = [tstart: delta: tend;] # Sequence from 0.0 - 10.0
    tindex = 2  # Initial step is already defined.

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


@doc raw"""
`model(u, p, t)`

Propensity model for this reaction.
Reaction of `A <-> B` with rate constants `k1` & `k2`.

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
function model(u, p, t)
    return [p.r * u[1],  # Reaction 1: Replicate wild-type mtDNA.
            p.d * u[1],  # Reaction 2: Degrade wild-type mtDNA.
            p.m * u[1],  # Reaction 3: Mutate wild-type mtDNA.
            p.r * u[2],  # Reaction 4: Replicate mutant mtDNA.
            p.d * u[2]]  # Reaction 5: Degrade mutant mtDNA.
end


"""
choose_stoich(dx, dxsum)

Choose and return which Stoichiometry to update the state.
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
    for i in [1: 1: n;]
        # Define and run the simulation.
        sol = ssa(model, u0, tend, parameters, choose_stoich)

        # - Take the results and store them in 'time_states' and
        #   'concentration_states' respectively.
        # TODO Move these operations to separate threads.
        push!(time_states, sol.t)
        push!(concentration_states, sol.u)
    end

    # TODO Cast concentrations as Int.
    # TODO Change instances of Int64 to Int.

    return(t = time_states, c = concentration_states)
end

"""
create_composite!(a::Array)

Mutate a multi-dimensional array and return a single array made up of
composite values.
"""
function create_composite!(a::Array)
    # TODO For each list in 't', calculate average and return single instance.
    # TODO For each list in 'c', calculate average and return single instance.
    # TODO Choose metric to return: median, mean, etc.

    a = mean(a) # FIXME Process converts Int to Float.
    # a = collect(Iterators.flatten(a))
    # a = sum(a, 1) ./ sum(a, 2)

    println("******************** INPUT SUMMARY ********************\n")
    println(summary(a))
    println(typeof(a))
    println(size(a))

    # User passed the concentration matrix.
    if isa(a, Matrix)
        println("Detected a matrix.")

        for i in eachindex(a)
            # TODO Force conversion of floating points into integers.
            # a[i] = round(a[i])  # Round to whole number.
            # a[i] = convert(Int64, a[i])  # Force conversion to Int.
            print("I=$(a[i, 1]) ")
            print("J=$(a[i, 2]) ")
        end

    # User passed the time-state vector.
    elseif isa(a, Vector)
        println("Detected a vector.")

        # Normalize the list of time states.
        a = sum(a[1:length(a)]) ./ length(a)
    end

    println()

    return(a)
end

#############
# Kickstart #
#############

## OLD CODE ################################################################################

# # Perform the simulation and assign results to 'sol'.
# print("\nNow performing simulation... ")
# sol = ssa(model, u0, tend, parameters, choose_stoich)
# println("Done!")

# # Define the plot.
# print("Now creating the plot... ")
# fig = plot(
#     sol.t,
#     sol.u,
#     xlabel="time",
#     ylabel="# of molecules",
#     title = "SSA",
#     label=["Wild-type" "Mutant"],
#     dpi=300)

# # Output plot to 'plot.png'.
# savefig(fig, "plot.png")
# println("Done!")

## NEW CODE ################################################################################

# Run simulation 'n' time(s).
results = loop_simulation(1000)

# # Print head of new arrays.
# # NOTE Assumes at least five entries exist.
# println(results.t[1:5])
# println(results.c[1:5])

# println("********** RESULTS (TIMES, UNEDITED) **********")
# println(results.t[1])

# println("********** RESULTS (CONCENTRATIONS, UNEDITED) **********")
# println(results.c[1])

x = create_composite!(results.t)  # Temporal states
y = create_composite!(results.c)  # Concentration states

println("******************** RESULTS (TIMES) ********************\n")
println(x)
println()

println("******************** RESULTS (CONCENTRATIONS) ********************\n")
println(y)
println()

print("Now creating the plot... ")
fig = plot(
    x,
    y,
    xlabel="Time",
    ylabel="Number of Molecules",
    xlims=(0, 10),
    ylims=(10, 300),
    title="mtDNA Population Dynamics (SSA Model)",
    label=["Wild-type" "Mutant"],
    dpi=300)
println("Done")

print("Now writing plot to 'plot_2.png'... ")
savefig(fig, "plot_2.png")
println("Done")

# End of File.
