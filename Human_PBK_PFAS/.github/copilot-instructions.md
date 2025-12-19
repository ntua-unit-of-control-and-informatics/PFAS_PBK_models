# PFAS Pharmacokinetic (PBK) Model - AI Coding Instructions

## Project Overview
This is an R-based physiologically-based pharmacokinetic (PBK) modeling system simulating PFAS (per- and polyfluoroalkyl substances) distribution in human and animal organisms. The model tracks chemical fate through absorption, distribution, metabolism, and excretion across multiple tissue compartments.

## Architecture & Data Flow

### Core Model Structure
**Main simulation files:**
- [PFAS_PBK_simulation.R](PFAS_PBK_simulation.R) - Generic PFAS model template with manual parameter adjustment
- [PFOA_PBK_simulation.R](PFOA_PBK_simulation.R) - PFOA-specific variant (identical structure, different defaults)
- [Abraham_Worley_PFAS_trial.R](Abraham_Worley_PFAS_trial.R) - Human trial data fitted model

Each follows the same three-step pattern:
1. **create.params()** - Builds parameter list from user inputs and constants
2. **create.inits()** - Initializes ODE state variables (all compartments start at 0)
3. **ode.func()** - System of differential equations defining tissue kinetics

### Tissue Compartments & Flow
The model divides the body into interconnected compartments with defined blood flows (fraction of cardiac output):
- **Liver (L)** - Primary sink; 25% cardiac output (`QLC = 0.25`)
- **Kidney (K)** - Active transport site; 17.5% cardiac output (`QKC = 0.175`)
- **Proximal Tubule Cells (PTC)** - Kidney subunit; site of active secretion via OAT transporters
- **Rest of Body (R)** - Remaining tissue; receives remaining cardiac output
- **GI Tract** - Stomach (ST) → Small Intestine (SI) for oral absorption

All exchanges follow Michaelis-Menten kinetics for active transport; passive diffusion and blood-flow limited for tissue distribution.

### Parameter Scaling Pattern
Parameters use allometric scaling to body weight (BW):
- **Flow rates** scale by `BW^0.75` (e.g., `QC = QCC * BW^0.75`)
- **Rate constants** scale by `BW^(-0.25)` (e.g., `kabs = kabsc * BW^(-0.25)`)
- **Tissue volumes** scale linearly by `BW`

This pattern is applied consistently across all rate constants (absorption, elimination, transport).

## Key Parameters & Conventions

### Always Present Parameters
- **BW** - Body weight (kg); drives all scaling
- **MW** - Molecular weight (g/mol); PFOA = 414.07
- **Free** - Free fraction in plasma (0.001 for PFOA; small value reflects high protein binding)
- **PL, PK, PR** - Partition coefficients (tissue:blood); stored in [PCs.csv](PCs.csv)

### Transport Kinetics (Michaelis-Menten)
- **Vmax_baso** - Basolateral transporter max velocity (OAT1/OAT3)
- **Vmax_apical** - Apical transporter max velocity (OAT4)
- **Km_baso, Km_apical** - Saturation constants
- **RAFbaso, RAFapi** - Relative Activity Factors; males differ from females/females differ from males (see sensitivity analysis)

### Elimination Routes
- **Biliary** (`kbile`) - Via liver to feces
- **Renal/Urinary** (`kurine`) - Via kidney filtrate and PTC secretion
- **Unabsorbed** (`kunabs`) - Oral dose not absorbed in GI tract

## Code Patterns & Conventions

### Parameter Organization
Group related parameters with inline comments citing literature sources (e.g., `# Davies 1993`). Balance checks included:
```r
QBal <- QC - (QK + QL + QR)  # Should equal zero
VBal <- (0.93*BW) - (VR + VL + VPTC + Vfil + VPlas)  # Should equal zero
```

### ODE Function Output
Returns list with two elements:
1. **Rate equations** (derivatives starting with `d`)
2. **Auxiliary outputs** - Mass balances (`Atissue`, `Aloss`, `Atotal`) and derived concentrations (`CR`, `CA`, `CPTC`, etc.)

Mass balance check is critical: `Atotal = Atissue + Aloss` should equal input dose throughout simulation.

### Events-Based Input
Administration handled via `deSolve::events()`:
- **IV bolus** - Direct injection to `Aplas_free` (plasma compartment)
- **Oral bolus** - Injection to `AST` (stomach compartment)
- **Water ingestion** - Continuous via `Cwater` and `ingestion` state variables

Input structure: `admin_type` (character), doses/times as vectors with equal length.

## Development Workflows

### Running Simulations
```r
# 1. Define user inputs
user_input <- list("BW" = 82, "admin_type" = "bolus", "admin_dose_bolus" = c(3.96), ...)
# 2. Create parameters, initialize states
params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)
# 3. Solve ODE system
solution <- ode(times = seq(0, 60, 0.01), func = ode.func, y = inits, parms = params, 
                 events = events, method = "lsoda", rtol = 1e-05, atol = 1e-05)
```

### Sensitivity Analysis
[Abraham_Worley_PFAS_sensitivity.R](Abraham_Worley_PFAS_sensitivity.R) uses `sensitivity` package for parameter screening. When modifying kinetic parameters, update both trial and sensitivity scripts.

### Experimental Data Fitting
- Plasma: [exp_data_plasma.csv](exp_data_plasma.csv) - Multiple PFAS compounds (PFBA, PFOA, PFOS, etc.)
- Urine: [exp_data_urine.csv](exp_data_urine.csv)
- Feces: [exp_data_feces.csv](exp_data_feces.csv)

Use `cumulative_exp_data()` function to compute cumulative excretion from time-concentration pairs.

## Testing & Validation

- Check mass balance (`Atotal` should remain constant = input dose)
- Compare predicted vs. experimental concentrations across all timepoints
- Verify parameter scaling sensibility (doubling BW should scale flows by 2^0.75 ≈ 1.68×)
- Sensitivity analysis identifies high-impact parameters for optimization

## File Dependencies

- [Chemical names and doses.csv](Chemical%20names%20and%20doses.csv) - PFAS chemical metadata
- [Clearance_data.csv](Clearance_data.csv) - Reference clearance values
- [Data_on_kinetics_of_oral_absorption_during_the_first_days.csv](Data_on_kinetics_of_oral_absorption_during_the_first_days.csv) - Early absorption kinetics
- CSV files loaded via `read.csv()` in trial/sensitivity scripts

## Important Notes for Modifications

1. **Parameter coupling**: Changes to molecular weight (MW) must propagate through Vmax calculations (involves MW scaling)
2. **Allometric scaling**: Any new rate constants should follow `BW^(-0.25)` pattern unless literature specifies otherwise
3. **Units consistency**: All concentrations in µg/L; volumes in L; masses in g; times in days; flows in L/day
4. **Literature tracking**: Every non-fitted parameter includes citation (e.g., Brown 1997) for reproducibility
