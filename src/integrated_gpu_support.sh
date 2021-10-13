#!/usr/bin/env bash

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

# - [Makie Backends & Output](https://tinyurl.com/efrfwz2t)

#########
# Notes #
#########

# - This script requires 'sudo' permissions and assumes you are using a
#   machine that has an integrated GPU e.g. Intel HD/UHD chipsets.
# - Users of an NVIDIA/AMD GPU should not need to use this script.

#############
# Variables #
#############

# VARIABLES > COLOURS

readonly green="\e[32m"   # Success
readonly red="\e[31m"     # Error
readonly yellow="\e[33m"  # Warning / Note
readonly reset="\e[0m"    # Reset text

# VARIABLES > PATHS

readonly model_file="agent_julia_model.jl"

#############
# Kickstart #
#############

# Execute the simulation if the model file exists at this location.
if [ -f "$model_file" ]; then
  echo -e "${yellow}NOTE${reset} Using 'sudo' permissions to execute '$model_file'."
  echo -e "${yellow}NOTE${reset} Enforcing the use of integrated GPU for plotting."
  echo
  sudo DRI_PRIME=1 julia "$model_file"
else
  echo -e "${red}ERROR${reset} Simulation file '$model_file' does not exist."
fi

# End of File.
