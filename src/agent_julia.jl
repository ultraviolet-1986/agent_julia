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
# Using #
#########

using Pkg

###########
# Imports #
###########

Pkg.add("Agents")

#############
# Variables #
#############

#############
# Functions #
#############

"""
agent_julia_main
Main function wrapper for the `agent_julia.jl` script.
"""
function agent_julia_main()
    println("Hello, World!")
    return
end

###########
# Exports #
###########

# By exporting these functions, it is possible to load their docstrings
# from the REPL by using the `include("agent_julia.jl")` command (for
# this file). From here, use the '?' key and type the name of a
# function to display its docstring.

export agent_julia_main

#############
# Kickstart #
#############

agent_julia_main()

# End of File.
