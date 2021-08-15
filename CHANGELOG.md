# Changelog

## Version 0.0.1

- Updated `gillespie_model.jl` to contain and plot results of a basic stochastic
  simulation of mtDNA population dynamics.
  - Included use-case scripts which contain elements simulating different
    situations. For example, the `patient_with_inheritance.jl` file contains
    parameters intended to mimic a patient who has inherited mutant mtDNA since
    birth.
- Included example plots to view what results may look like per-situation.
- Plots are now stored within a well-defined directory structure.
- Included `integrated_gpu_support.sh` to correct an update for the
  `CairoMakie` package, which no longer supports integrated GPUs by default.
  - This script will execute the `agent_julia_model.jl` with `sudo` permissions
    and will also be used to change permissions for the output. Review code
    before use.
- `Agent Julia` model can now loop simulations in the same manner as the
  Gillespie SSA model. Because of the complexity of the model, it does take much
  longer to complete the simulation at `n = 1000`.
