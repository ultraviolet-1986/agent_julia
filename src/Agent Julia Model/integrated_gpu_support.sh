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

# - <https://makie.juliaplots.org/stable/backends_and_output.html>

#########
# Notes #
#########

# - This script requires 'sudo' permissions and assumes you are using a
#   machine that has an integrated GPU e.g. Intel HD/UHD chipsets.
# - Users of an NVIDIA/AMD GPU should not need to use this script.

#############
# Variables #
#############

readonly model_file='agent_julia_model.jl'
readonly simulation_video='agent_julia_simulation.mp4'

#############
# Kickstart #
#############

# Delete previous video if exists.
if [ -f "$simulation_video" ]; then
  rm "$simulation_video"
  echo -e "NOTE Deleted pre-existing '$simulation_video' file."
fi

# Execute the simulation if the model file exists at this location.
if [ -f "$model_file" ]; then
  echo -e "NOTE Using 'sudo' permissions to execute '$model_file'."
  echo -e "NOTE Enforcing the use of integrated GPU for plotting.\n"
  sudo DRI_PRIME=1 julia "$model_file"
else
  echo -e "ERROR: Simulation file '$model_file' does not exist."
fi

# Change ownership of 'agent_julia_simulation.mp4' from 'root' to the
# current user if the file exists.
if [ -f 'agent_julia_simulation.mp4' ]; then
  echo -e "\nNOTE Using 'sudo' to change 'agent_julia_simulation.mp4' permissions.\n"
  sudo chown "$USER":"$USER" 'agent_julia_simulation.mp4'
fi

# End of File.
