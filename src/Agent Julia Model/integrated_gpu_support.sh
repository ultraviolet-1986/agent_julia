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
# Kickstart #
#############

echo "NOTE Using 'sudo' permissions to execute 'agent_julia_model.jl'."
echo "NOTE Enforcing the use of integrated GPU for plotting."
sudo DRI_PRIME=1 julia 'agent_julia_model.jl'

# Change ownership of 'bacteria.mp4' from 'root' to the current user if
# the file exists.
if [ -f 'bacteria.mp4' ]; then
  echo "NOTE Using 'sudo' to change 'bacteria.mp4' permissions."
  sudo chown "$USER":"$USER" 'bacteria.mp4'
fi

# End of File.
