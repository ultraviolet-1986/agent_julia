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

# - [Agents.jl Tutorial](https://git.io/Jc1w6)

###########
# Imports #
###########

import Pkg

Pkg.add("Agents")
Pkg.add("DrWatson")
Pkg.add("CairoMakie")
Pkg.add("InteractiveDynamics")
Pkg.add("Random")

using Agents
using DrWatson: @dict
using CairoMakie
using InteractiveDynamics
using Random

#################
# Prerequisites #
#################

# Share seed throughout program.
const seed = 41269

# Ensure results are reproducible.
Random.seed!(seed)

#############
# Variables #
#############

const steps_per_day = 24

###########
# Structs #
###########

mutable struct Agent <: AbstractAgent
    id::Int
    pos::NTuple{2,Float64}
    vel::NTuple{2,Float64}
    mass::Float64
end

mutable struct PoorSoul <: AbstractAgent
    id::Int
    pos::NTuple{2,Float64}
    vel::NTuple{2,Float64}
    mass::Float64
    days_infected::Int  # number of days since is infected
    status::Symbol  # :S, :I or :R
    β::Float64
end

#############
# Functions #
#############

function ball_model(; speed = 0.002)
    space2d = ContinuousSpace((1, 1), 0.02)
    model = ABM(Agent, space2d, properties = Dict(:dt => 1.0), rng = MersenneTwister(42))

    # And add some agents to the model
    for ind in 1:500
        pos = Tuple(rand(model.rng, 2))
        vel = sincos(2π * rand(model.rng)) .* speed
        add_agent!(pos, model, vel, 1.0)
    end
    return model
end

using DrWatson: @dict
function sir_initiation(;
    infection_period = 30 * steps_per_day,
    detection_time = 14 * steps_per_day,
    reinfection_probability = 0.05,
    isolated = 0.0, # in percentage
    interaction_radius = 0.012,
    dt = 1.0,
    speed = 0.002,
    death_rate = 0.044, # from website of WHO
    N = 1000,
    initial_infected = 5,
    seed = 42,
    βmin = 0.4,
    βmax = 0.8,
)

    properties = @dict(
        infection_period,
        reinfection_probability,
        detection_time,
        death_rate,
        interaction_radius,
        dt,
    )
    space = ContinuousSpace((1,1), 0.02)
    model = ABM(PoorSoul, space, properties = properties, rng = MersenneTwister(seed))

    # Add initial individuals
    for ind in 1:N
        pos = Tuple(rand(model.rng, 2))
        status = ind ≤ N - initial_infected ? :S : :I
        isisolated = ind ≤ isolated * N
        mass = isisolated ? Inf : 1.0
        vel = isisolated ? (0.0, 0.0) : sincos(2π * rand(model.rng)) .* speed

        # very high transmission probability
        # we are modelling close encounters after all
        β = (βmax - βmin) * rand(model.rng) + βmin
        add_agent!(pos, model, vel, mass, 0, status, β)
    end

    return model
end


function transmit!(a1, a2, rp)
    # for transmission, only 1 can have the disease (otherwise nothing happens)
    count(a.status == :I for a in (a1, a2)) ≠ 1 && return
    infected, healthy = a1.status == :I ? (a1, a2) : (a2, a1)

    rand(model.rng) > infected.β && return

    if healthy.status == :R
        rand(model.rng) > rp && return
    end
    healthy.status = :I
end

function sir_model_step!(model)
    r = model.interaction_radius
    for (a1, a2) in interacting_pairs(model, r, :nearest)
        transmit!(a1, a2, model.reinfection_probability)
        elastic_collision!(a1, a2, :mass)
    end
end


function sir_agent_step!(agent, model)
    move_agent!(agent, model, model.dt)
    update!(agent)
    recover_or_die!(agent, model)
end

update!(agent) = agent.status == :I && (agent.days_infected += 1)

function recover_or_die!(agent, model)
    if agent.days_infected ≥ model.infection_period
        if rand(model.rng) ≤ model.death_rate
            kill_agent!(agent, model)
        else
            agent.status = :R
            agent.days_infected = 0
        end
    end
end

#############
# Kickstart #
#############

model = ball_model()

sir_model = sir_initiation()

sir_colors(a) = a.status == :S ? "#2b2b33" : a.status == :I ? "#bf2642" : "#338c54"

fig, abmstepper = abm_plot(sir_model; ac = sir_colors)
fig # display figure


sir_model = sir_initiation()

abm_video(
    "agent_julia_simulation.mp4",
    sir_model,
    sir_agent_step!,
    sir_model_step!;
    title = "mtDNA population dynamics",
    frames = 1000,
    ac = sir_colors,
    as = 10,
    spf = 1,
    framerate = 60,
)

# End of File.
