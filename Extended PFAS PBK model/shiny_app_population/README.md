# PFAS PBK Model — Population Shiny App

A population (virtual-individual) version of `../shiny_app`: instead of running the PBK
model once for a single fixed body weight, each run generates `n` popgen virtual
individuals (organ mass/flow, sex-specific hematocrit, age-corrected renal clearance) and
solves the model once per individual. Forward Dosimetry produces a *distribution* of
tissue-concentration trajectories; Reverse Dosimetry produces a *distribution* of
exposure estimates that reproduce a target point of departure (POD).

Built on `Extended_Model_Population.r`, `Forward_Dosimetry_Stochastic.r`, and
`Reverse_Dosimetry_Stochastic.r` (population/stochastic counterparts of the deterministic
`shiny_app`'s `Extended_Model.r`, `forward_dosimetry.r`, and `reverse_dosimetry.r`).

## Required Packages

Install all dependencies from the R console:

```r
install.packages(c(
  "shiny",
  "deSolve",
  "ggplot2",
  "tidyverse",
  "bslib",
  "nloptr",
  "plotly"
))
```

## Running the App

Navigate to this directory and run:

```bash
Rscript -e "shiny::runApp('.')"
```

To open the app automatically in your browser:

```bash
Rscript -e "shiny::runApp('.', launch.browser = TRUE)"
```

To run on a specific port:

```bash
Rscript -e "shiny::runApp('.', port = 8080)"
```

## Project Structure

```
shiny_app_population/
├── app.r                              # Main app (UI + server)
├── R/
│   ├── Extended_Model_Population.r    # Population ODE model (popgen physiology, params, inits, events, ODE function)
│   ├── Forward_Dosimetry_Stochastic.r # Population forward dosimetry driver
│   └── Reverse_Dosimetry_Stochastic.r # Population reverse dosimetry optimizer
├── popgen_consts.R                    # PopGen reference constants (ICRP/P3M tables)
├── popgen_port.R                      # PopGen population-generation algorithm
├── estimated_parameters.csv           # Chemical-specific PK parameters (MW, Free, Vmax/Km/RAF, ...)
└── README.md
```

The three files in `R/` are trimmed copies of the originals (everything up to, but not
including, their own `# --- Example usage` demo block) — Shiny auto-sources every `.R`
file in `R/` on startup, and the full originals each end with a runnable ~100-person demo
simulation (and, for `Extended_Model_Population.r`, a `ggsave()` call and a read of a CSV
that doesn't exist here) that would otherwise fire on every app launch. If you update the
originals in the parent `Extended PFAS PBK model/` directory, regenerate these copies by
re-copying everything before their `# --- Example usage` marker.

`popgen_consts.R`, `popgen_port.R`, and `estimated_parameters.csv` must stay at this
directory's root (not inside `R/`): `Extended_Model_Population.r`'s own `source()`/
`read.csv()` calls for these use paths relative to the working directory, and Shiny keeps
the working directory at the app root (not `R/`) throughout.

## Differences from `../shiny_app`

- **No single body weight.** Each tab instead exposes the population-generation inputs
  `Extended_Model_Population.r`'s `.generate_population()` requires: population size,
  reference dataset (ICRP/P3M), population variability (Realistic/HighVariation), age/
  BMI/height ranges, probability male, ethnicity probabilities, and a random seed.
- **No intravenous route.** `Extended_Model_Population.r`'s event mechanism no longer has
  an `admin_type` (iv/bolus/oral) switch — every dose is expressed per kg body weight and
  routed through `exp_type`: `"continuous"` for a steady intake rate, `"pharmacokinetics"`
  for a one-off dose. Both are scaled to each individual's own popgen-derived body weight.
- **Population outputs.** Forward Dosimetry plots show the population median with a
  2.5th–97.5th percentile band (or, for the all-compartments plot, the population median
  per tissue) instead of a single trajectory. Reverse Dosimetry reports the *distribution*
  of per-individual exposure estimates (summary stats + histogram + a per-individual
  table) rather than one point estimate, and does not plot concentration-time curves —
  the underlying search records each individual's final estimated exposure and how well
  it matched the POD, not the full time-course at that exposure.
- **Runtime.** Reverse Dosimetry runs a full numerical search per individual (up to
  ~1000 ODE solves each), so it is much slower per-individual than Forward Dosimetry.
  Start with a small population (≤ 15) there before scaling up.
- **Interactive plots.** All plots use `plotly::ggplotly()` (hover tooltips, zoom/pan,
  legend click-to-toggle) instead of static `ggplot2` images.
- **POD measurement uncertainty (Reverse Dosimetry).** The "POD standard deviation" input
  (default 0) lets the POD itself carry uncertainty: above 0, each individual searches
  against its own random POD draw (lognormal, mean = POD, this SD) instead of the exact
  same target for everyone, so the exposure distribution reflects POD uncertainty on top
  of physiological variability. See `Reverse_Dosimetry_Stochastic.r`'s `POD_sd` parameter
  and its header's "POD uncertainty note".
