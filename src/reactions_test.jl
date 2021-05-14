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

# - This script provides a proof-of-concept for performing basic
#   reactions to be used within Agent Julia's model.
# - Because Julia changes the value of an object in-place, two separate
#   sets had to be declared for reaction 1 and 2, then for 3, 4, and 5
#   respectively.

#############
# Variables #
#############

# Used for:
# - Reaction 1: replicate_wildtype_mtdna
# - Reaction 2: degrade_wildtype_mtdna
set1 = ["M₁"]

# Used for:
# - Reaction 3: mutate_wildtype_mtdna
# - Reaction 4: replicate_mutant_mtdna
# - Reaction 5: degrade_mutant_mtdna
set2 = ["M₁"]

# Text formatting > Colours
GREEN = "\e[32m"

# Text formatting > Weights
BOLD = "\e[1m"
RESET = "\e[0m"

#############
# Functions #
#############

"""
replicate_wildtype_mtdna(M)

Take M₁ as input, perform replication, and return two instances of M₁.
"""
function replicate_wildtype_mtdna(M::Array)
    println("$(BOLD)Replication of wild-type mtDNA$(RESET)")
    println("$(BOLD)Input:$(RESET)    $(M)")
    println("$(BOLD)Process:$(RESET)  rM₁")

    # push!(M, "M₁")  # Append "M₁" to M.
    # push!(M, M[])   # Duplicate array contents.
    push!(M, M[1])    # Duplicate initial array content.

    println("$(BOLD)Output:$(RESET)   $(GREEN)$(M)$(RESET)\n")
    return M
end


"""
degrade_wildtype_mtdna(M)

Take M₁ and M₁ as input, perform degradation, and return ∅.
"""
function degrade_wildtype_mtdna(M::Array)
    println("$(BOLD)Degradation of wild-type mtDNA$(RESET)")
    println("$(BOLD)Input:$(RESET)    $(M)")
    println("$(BOLD)Process:$(RESET)  dM₁")

    # pop!(M)  # Remove single element.
    empty!(M)  # Drop entire list.

    println("$(BOLD)Output:$(RESET)   $(GREEN)$(M)$(RESET)\n")
    return M
end


"""
mutate_wildtype_mtdna(M)

Take M₁ as input, perform mutation, and return M₁ and M₂.
"""
function mutate_wildtype_mtdna(M::Array)
    println("$(BOLD)Mutation of wild-type mtDNA$(RESET)")
    println("$(BOLD)Input:$(RESET)    $(M)")
    println("$(BOLD)Process:$(RESET)  mM₁")

    # replace.(M, r"^M₁"=>"M₂")  # Replace 'M₁' with 'M₂' elements.
    push!(M, "M₂")

    println("$(BOLD)Output:$(RESET)   $(GREEN)$(M)$(RESET)\n")
    return M
end


"""
replicate_mutant_mtdna(M)

Take M₁ and M₂ as input, perform replication, and return two instances
of M₂.
"""
function replicate_mutant_mtdna(M::Array)
    println("$(BOLD)Replication of mutant mtDNA$(RESET)")
    println("$(BOLD)Input:$(RESET)    $(M)")
    println("$(BOLD)Process:$(RESET)  rM₂")

    push!(M, M[2])  # Replicate M₂.
    splice!(M, 1)   # Remove initial M₁.

    println("$(BOLD)Output:$(RESET)   $(GREEN)$(M)$(RESET)\n")
    return M
end


"""
degrade_mutant_mtdna(M)

Take M₂ and M₂ as input, perform degradation, and return ∅.
"""
function degrade_mutant_mtdna(M::Array)
    println("$(BOLD)Degration of mutant mtDNA$(RESET)")
    println("$(BOLD)Input:$(RESET)    $(M)")
    println("$(BOLD)Process:$(RESET)  dM₂")

    # pop!(M)  # Remove single element.
    empty!(M)  # Drop entire list.

    println("$(BOLD)Output:$(RESET)   $(GREEN)$(M)$(RESET)\n")
    return M
end

#############
# Kickstart #
#############

# Perform reactions on 'set1'.
replicate_wildtype_mtdna(set1)
degrade_wildtype_mtdna(set1)

# Perform reactions on 'set2'.
mutate_wildtype_mtdna(set2)
replicate_mutant_mtdna(set2)
degrade_mutant_mtdna(set2)

# End of File.
