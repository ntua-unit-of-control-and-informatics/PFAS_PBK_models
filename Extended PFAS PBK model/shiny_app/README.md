# PFAS PBK Model — Shiny App

A physiologically based kinetic (PBK) model for PFAS compounds with forward and reverse dosimetry modes.

## Required Packages

Install all dependencies from the R console:

```r
install.packages(c(
  "shiny",
  "deSolve",
  "ggplot2",
  "tidyverse",
  "bslib",
  "nloptr"
))
```

## Running the App

Navigate to the project directory and run:

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
shiny_apps/
├── app.r                      # Main app (UI + server)
├── R/
│   ├── Extended_Model.r       # ODE model (params, inits, events, ODE function)
│   ├── forward_dosimetry.r    # Forward dosimetry functions
│   └── reverse_dosimetry.r    # Reverse dosimetry optimizer
└── README.md
```
