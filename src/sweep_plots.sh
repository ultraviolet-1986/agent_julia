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

#########
# Notes #
#########

# - This script will remove any and all PNG files within the current
#   folder only.
# - This script does not outright delete files, it places them within
#   the Wastbasket.

#############
# Functions #
#############

sweep_plots_main(){
  for x in *; do
    if [ -d "$x" ]; then
      cd "$x" || return 1

      for y in *.png; do
        if [ -f "$y" ]; then
          printf "Now moving '%s/%s' to Wastebasket... " "$x" "$y"
          gio trash "$y"
          sync
          echo "Done"
        fi
      done

      cd ~- || return 1
    fi
  done

  return 0
}

#############
# Kickstart #
#############

sweep_plots_main

# End of File.
