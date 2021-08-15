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

#############
# Variables #
#############

# VARIABLES > INITIAL CONDITIONS

initial_mutants = 25

# VARIABLES > PATHS

data_path = "$(Base.source_dir())/data"
plot_path = "$(Base.source_dir())/plots/patient_with_inheritance"
video_path = "$(Base.source_dir())/videos"

csv_path = "$(data_path)/patient_with_inheritance.csv"
plot_1_path = "$(plot_path)/01_mutation_levels.png"
plot_2_path = "$(plot_path)/02_total_copy_average.png"
mp4_path = "$(video_path)/patient_with_inheritance.mp4"

#############
# Kickstart #
#############

# Create output directories.
mkpath(data_path)
mkpath(plot_path)
mkpath(video_path)

# EXECUTE MODEL

# Run Agent Julia model with above parameters.
# NOTE Static seed will be assigned from this file.
include("$(Base.source_dir())/agent_julia_model.jl")

# VARIABLES > PLOT AXES

times = eachcol(data)[1] / year        # Years (Plot 1 and 2)
mutation_level = eachcol(data)[3]      # N     (Plot 1)
total_copy_numbers = eachcol(data)[4]  # N     (Plot 2)

# DEFINE MUTATION LEVELS PLOT (PLOT 1)

print("Building mutation levels plot. Please wait... ")
fig = plt.plot(
    times,
    mutation_level,
    xlims=(0, 80),
    ylims=(0, agent_max),
    title="Patient with mutant mtDNA inheritance",
    xlabel="Time (years)",
    ylabel="Mutation levels (n, mean)",
    legend=false,
    dpi=1200,
)
println("$(green)Done$(reset)")

print("Rendering plot to $(yellow)$(plot_1_path)$(reset)... ")
plt.savefig(fig, "$(plot_1_path)")
println("$(green)Done$(reset)")

# DEFINE TOTAL AVERAGE COPY NUMBER PLOT (PLOT 2)

print("Building total average copy number plot. Please wait...")
fig2 = plt.plot(
    times,
    total_copy_numbers,
    xlims=(0, 80),
    ylims=(0, agent_max),
    title="Patient with mutant mtDNA inheritance",
    xlabel="Time (years)",
    ylabel="Total agent levels (n, mean)",
    legend=false,
    dpi=1200,
)
println("$(green)Done$(reset)")

print("Rendering plot to $(yellow)$(plot_2_path)$(reset)... ")
plt.savefig(fig2, "$(plot_2_path)")
println("$(green)Done$(reset)")

# DEFINE CSV EXPORT FILE

print("Writing data to $(yellow)$(csv_path)$(reset) Please wait... ")
CSV.write("$(csv_path)", data)
println("$(green)Done$(reset)")

# CREATE OPTIONAL VIDEO FILE

print("Defining simulation colour palette... ")
model_colours(a) = a.status == :W ? green_hex : red_hex
println("$(green)Done$(reset)")

print("Rendering simulation output as $(yellow)$(mp4_path)$(reset)... ")
try
    abm_video(
        "$(mp4_path)",
        agent_julia_model,
        agent_step!,
        model_step!;
        title = "Patient with mutant inheritance",
        frames = tend,
        ac = model_colours,
        as = 10,
        spf = Int(round(δ)),
        framerate = 60,
    )
    println("$(green)Done$(reset)")
catch
    println("$(red)Error$(reset)")
    println("Please run the $(yellow)integrated_gpu_support.sh$(reset) script.")
end

# End of File.
