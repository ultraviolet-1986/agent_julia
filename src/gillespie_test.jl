#!/usr/bin/env julia

# Reference: https://nextjournal.com/bebi5009/gillespie-julia

#################
# Prerequisites #
#################

import Pkg

Pkg.add("Plots")

using Random
using Plots

#############
# Functions #
#############

#=
Stochastic chemical reaction: Gillespie Algorithm
Adapted from: Chemical and Biomedical Enginnering Calculations Using Python Ch.4-3
=#
function ssa(model, u0, tend, p, choose_stoich, tstart=zero(tend))
    t = tstart    # Current time
    ts = [t]      # Time points
    u = copy(u0)  # Current state
    us = copy(u)  # Record of states

    while t < tend
        dx = model(u, p, t)              # propensities of reactions
        dt = Random.randexp() / sum(dx)  # Time step
        stoich = choose_stoich(dx)       # Get stoichiometry
        u .+= stoich
        t += dt

        # Add to record
        us = [us u]
        push!(ts, t)  # Record t
    end

    us = collect(us')

    return (t = ts, u = us)
end

#=
Reaction of A <-> B with rate constants k1 & k2
=#
"Propensity model for this reaction"
model(u, p, t) = [p.k1 * u[1],  p.k2 * u[2]]

"Choose and return which Stoichiometry to update the state"
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

u0 = [200, 0]
tend = 10.0
parameters = (k1=1.0, k2=0.5)

sol = ssa(model, u0, tend, parameters, choose_stoich)

fig = plot(sol.t, sol.u, xlabel="time", ylabel="# of molecules", title = "SSA", label=["A" "B"])

savefig(fig, "plot.png")

# End of File.
