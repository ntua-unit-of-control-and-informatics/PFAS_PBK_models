# Rainbow Trout PFAS Physiologically Based Kinetic (PBK) Model

## Overview

This model simulates the uptake, distribution, and elimination of per- and polyfluoroalkyl substances (PFAS) in rainbow trout (*Oncorhynchus mykiss*) following dietary exposure. It is a **Physiologically Based Kinetic (PBK) model**, meaning it describes how a chemical moves through the body based on realistic anatomical and physiological information — blood flows, organ sizes, and tissue-specific chemical affinity — rather than purely statistical relationships.

The model can simulate five PFAS compounds: **PFOS**, **PFOA**, **PFBS**, **PFHxS**, and **PFNA**.

---

## Model Structure and Compartments

The model represents the fish body as a network of interconnected compartments, each corresponding to a tissue or physiological space. Chemical mass is tracked in each compartment over time and moves between them via blood circulation, bile flow, urinary excretion, and dietary absorption.

### Tissue Compartments
| Compartment | Description |
|---|---|
| **Arterial blood** | Oxygenated blood delivering PFAS to tissues |
| **Venous blood** | Blood returning from tissues to the gills |
| **Gills** | Primary site of blood oxygenation; PFAS exchange with circulation |
| **Liver** | Central metabolic organ; receives blood from both arterial supply and viscera; source of bile secretion |
| **Viscera** | Abdominal organs (gut wall, mesentery); primary site of dietary absorption |
| **Gut lumen (absorbable fraction)** | Ingested PFAS available for absorption into the viscera |
| **Gut lumen (bile fraction)** | Bile-derived PFAS returned to the gut; only eliminated via feces |
| **Kidney** | Filters PFAS from blood; handles urinary excretion and reabsorption |
| **Urine storage** | Urinary bladder; temporary storage before voiding |
| **Muscle** | The largest tissue by mass; slow distribution kinetics |
| **Skin** | Peripheral tissue with direct routing of venous blood to kidney |
| **Carcass** | Remaining body mass (bone, connective tissue, etc.) |

### Elimination Routes
- **Urine**: PFAS filtered at the kidney, partially reabsorbed, and excreted via the urinary bladder.
- **Feces**: Unabsorbed dietary PFAS and bile-derived PFAS eliminated through the gut.
- **Enterohepatic circulation**: A fraction of PFAS secreted into bile is reabsorbed from the gut back into the viscera, prolonging tissue retention.

---

## Physiological Basis

Organ volumes and blood flow fractions are taken from **Vidal et al. (2019)**, measured in rainbow trout. Cardiac output scales with body weight and is corrected for water temperature using an Arrhenius function, with reference values from **Barron et al. (1987)** at 6, 12, and 18 °C. Chemical-specific parameters (protein binding, urinary clearance, intestinal absorption fractions) are sourced from **Sun et al. (2022)**, **Cao et al. (2022)**, and **Goeritz et al. (2013)**.

---

## Training Data and Parameter Estimation

The model was calibrated against experimental data from **Falk et al. (2015)**, who exposed rainbow trout to PFAS-contaminated food under controlled laboratory conditions. The experimental design consisted of:

- **28-day dietary uptake phase**: fish fed food spiked with 500 µg PFAS/kg at 2.6% of body weight per day.
- **28-day depuration phase**: fish returned to clean food.
- **Tissue samples** collected at days 7, 14, 28, 31, 35, 42, and 56.
- **Tissues measured**: liver, blood, skin, muscle, gills, kidney, and carcass.
- **Five substances** tested in parallel: PFOS, PFOA, PFBS, PFHxS, PFNA.

### Fitted Parameters

Parameter estimation was performed using the **Subplex (NLOPT_LN_SBPLX)** derivative-free optimization algorithm. The objective function was the **PBK Objective Function (PBKOF)**, which measures the geometric mean fold-error between model predictions and observations across all tissues and time points. A score of 1.0 is a perfect fit; values between 0.5 and 2.0 are generally considered acceptable.

The final PBKOF score across all five substances was **0.625**, indicating good agreement between model predictions and experimental observations.

Three parameters were estimated as **common across all five PFAS** (reflecting shared physiological processes):

| Parameter | Description | Fitted Value |
|---|---|---|
| **Ku** | Intestinal absorption rate constant (1/h) | 1.467 |
| **CLU_coef** | Urinary clearance scaling coefficient | 5.72 × 10⁻⁴ |
| **Cl_feces** | Fecal elimination rate constant (1/h) | 1.306 |

Seven **tissue-to-blood partition coefficients** were estimated separately for each substance, reflecting the different affinity of each PFAS for each tissue:

| Tissue | PFOS | PFOA | PFBS | PFHxS | PFNA |
|---|---|---|---|---|---|
| Liver | 1.569 | 2.004 | 1.741 | 1.698 | 0.803 |
| Muscle | 0.113 | 0.037 | 0.139 | 0.043 | 0.065 |
| Kidney | 0.440 | 0.851 | 0.763 | 0.361 | 0.246 |
| Skin | 0.272 | 0.319 | 0.224 | 0.293 | 0.234 |
| Gills | 0.229 | 0.343 | 0.277 | 0.154 | 0.220 |
| Carcass | 0.107 | 0.180 | 0.116 | 0.041 | 0.114 |
| Viscera | 3.699 | 0.564 | ~0 | ~0 | 1.275 |

---

## How to Use the Model

The model requires a small set of inputs that describe the fish and the exposure scenario. No programming knowledge is needed.

### Selecting the Substance
Choose which PFAS compound to simulate: **PFOS**, **PFOA**, **PFBS**, **PFHxS**, or **PFNA**. Each compound has its own chemical properties and fitted partition coefficients, so only one substance is simulated per run.

### Fish Body Weight
Provide the body weight of the fish in grams. The model was calibrated on fish weighing between **314 g and 808 g**, which represents the validated range. If a weight outside this range is entered, the model will automatically use the nearest boundary value to ensure predictions remain within the calibrated domain.

### Defining the Exposure Scenario
The exposure is described by two paired lists:
- **Dose** (`admin.dose`): the amount of PFAS ingested at each feeding event, in micrograms (µg).
- **Time** (`admin.time`): the time of each feeding event, in hours from the start of the simulation.

This flexible format allows you to simulate a wide range of exposure scenarios. For example:
- **Single acute dose**: one entry in each list representing a single exposure event.
- **Repeated daily feeding**: 28 entries with doses spaced 24 hours apart, mimicking a dietary uptake experiment.
- **Variable intake**: different dose values at each time point to reflect changing food contamination levels or feeding rates over time.

The two lists must always have the same number of entries.

### Simulation Time Window
Set the start time, end time, and output time step (all in hours) to define how long the simulation runs and how finely the results are reported.

### Model Outputs
The model returns time-series predictions of PFAS concentration (µg/g tissue) in all compartments, as well as the total mass of PFAS in each compartment (µg). These can be used to estimate tissue residues at any time point during or after exposure.

---

## References

- Falk, S. et al. (2015). *Bioconcentration of perfluoroalkyl acids in rainbow trout.* Environmental Toxicology and Chemistry.
- Vidal, A. et al. (2019). *Physiologically based toxicokinetic modelling of PFAS in rainbow trout.* Science of the Total Environment.
- Barron, M.G. et al. (1987). *Cardiovascular parameters in rainbow trout.* Journal of Fish Biology.
- Sun, M. et al. (2022). *Urinary and fecal elimination of PFAS in fish.* Environmental Science & Technology.
- Cao, X. et al. (2022). *Enterohepatic circulation of PFAS in fish.* Environmental Health Perspectives.
- Goeritz, I. et al. (2013). *Intestinal absorption fractions of PFAS.* Archives of Toxicology.
- Grosell, M. et al. (2000). *Biliary secretion in teleost fish.* Journal of Experimental Biology.
- Nichols, J.W. et al. (1996). *Physiological model for gill uptake of organic chemicals in fish.* Environmental Toxicology and Chemistry.
