#!/usr/bin/env julia

# This is an example file showing basic functionality from the
# Agents.jl tutorial page.

# This script will generate a "game_of_life.mp4" file containing the
# results of the model.

###########
# Imports #
###########

# Import and add packages.

import Pkg

Pkg.add("Agents")
Pkg.add("InteractiveDynamics")
Pkg.add("CairoMakie")

# Load imports and packages.

import CairoMakie

using Random
using Agents
using InteractiveDynamics

#############
# Variables #
#############

rules = (2, 3, 3, 3)

#############
# Functions #
#############

mutable struct Cell <: AbstractAgent
    id::Int
    pos::Dims{2}
    status::Bool
end

# function build_model(; rules::Tuple, dims = (100, 100), metric = :chebyshev, seed = 120)
function build_model(; rules::Tuple, dims = (100, 100), metric = :chebyshev, seed = 41269)
    space = GridSpace(dims; metric)
    properties = Dict(:rules => rules)
    model = ABM(Cell, space; properties, rng = MersenneTwister(seed))
    idx = 1
    for x in 1:dims[1]
        for y in 1:dims[2]
            add_agent_pos!(Cell(idx, (x, y), false), model)
            idx += 1
        end
    end
    return model
end

function ca_step!(model)
    new_status = fill(false, nagents(model))
    for agent in allagents(model)
        n = alive_neighbors(agent, model)
        if agent.status == true && (n ≤ model.rules[4] && n ≥ model.rules[1])
            new_status[agent.id] = true
        elseif agent.status == false && (n ≥ model.rules[3] && n ≤ model.rules[4])
            new_status[agent.id] = true
        end
    end

    for id in allids(model)
        model[id].status = new_status[id]
    end
end

function alive_neighbors(agent, model) # count alive neighboring cells
    c = 0
    for n in nearby_agents(agent, model)
        if n.status == true
            c += 1
        end
    end
    return c
end

#############
# Kickstart #
#############

model = build_model(rules = rules, dims = (50, 50))

for i in 1:nagents(model)
    if rand(model.rng) < 0.2
        model.agents[i].status = true
    end
end

ac(x) = x.status == true ? :black : :white
am(x) = x.status == true ? '■' : '□'
abm_video(
    "game of life.mp4",
    model,
    dummystep,
    ca_step!;
    title = "Game of Life",
    ac = :black,
    as = 12,
    am,
    framerate = 5,
    scatterkwargs = (strokewidth = 0,),
)

# End of File.
