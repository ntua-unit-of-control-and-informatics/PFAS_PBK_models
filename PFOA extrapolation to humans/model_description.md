# PFOA Physiologically Based Pharmacokinetic (PBPK) Model — Human Extrapolation

## Overview

This service provides a **whole-body physiologically based pharmacokinetic (PBPK) model** for
**perfluorooctanoic acid (PFOA)** in humans. PFOA (C₇F₁₅COOH, MW = 414.07 g/mol) is a
per- and polyfluoroalkyl substance (PFAS) that is persistent in the environment and
accumulates in the human body, primarily in serum, liver, and kidneys. The model predicts
time-dependent PFOA concentrations in blood, plasma, and 14 major organs/tissues following
different routes of exposure.

The model was developed by extrapolating from rodent data to humans using in-vitro to in-vivo
extrapolation (IVIVE) of transporter kinetics and was calibrated against human biomonitoring data
(Abraham et al., 2024).

---

## Model Structure

### Compartments

The body is divided into **14 tissue compartments** plus blood:

| Compartment | Abbreviation |
|---|---|
| Kidney | Ki |
| Liver | Li |
| Intestine | In |
| Stomach | St |
| Muscle | Mu |
| Adipose | Ad |
| Lungs | Lu |
| Spleen | Sp |
| Heart | Ht |
| Brain | Br |
| Gonads | Go |
| Skin | Sk |
| Bones | Bo |
| Rest (lumped remaining tissues) | Re |

Each compartment is further subdivided into three sub-compartments:

- **Blood capillaries** — receives arterial blood and returns venous blood
- **Interstitial fluid** — exchanges PFOA with capillaries and cells via diffusion
- **Intracellular / tissue space** — intracellular accumulation

The kidney additionally includes the **tubular filtrate compartments** (proximal tubule, loop of
Henle, distal tubule, collecting duct) and a **urinary bladder** compartment.
The liver includes a **bile** compartment linked to enterohepatic recirculation.
The lung includes an **alveolar lining fluid** compartment and an **upper airway** compartment.

### Blood

Blood is separated into arterial and venous pools. **Albumin–PFOA binding** and, in the
kidney tubular cells and liver, **L-FABP–PFOA binding** are modelled explicitly using
reversible saturable kinetics, tracking both the free and bound fractions at each site.

### Transport Processes

PFOA is transported by saturable carrier-mediated transporters scaled from in-vitro data
via IVIVE:

| Transporter | Location | Direction |
|---|---|---|
| OATP (SLCO) | Kidney proximal tubule | Uptake into cells |
| OAT1 (SLC22A6) | Kidney proximal tubule | Uptake into cells |
| OAT3 (SLC22A8) | Kidney proximal tubule | Uptake into cells |
| URAT1 (SLC22A12) | Kidney proximal tubule | Uptake into cells |
| OATP1B1/1B3 | Liver (basolateral) | Uptake into hepatocytes |
| OATP2B1 | Liver (basolateral) | Uptake into hepatocytes |
| NTCP (SLC10A1) | Liver (basolateral) | Uptake into hepatocytes |
| OATP | Lung epithelium (apical/basolateral) | Uptake |
| OATP2B1 | Intestinal epithelium | Uptake |

Passive transcellular diffusion and paracellular (gap-junction) transport are modelled for
every capillary–interstitium interface using the Renkin pore model.

### Elimination Routes

- **Renal excretion**: glomerular filtration followed by partial tubular reabsorption via
  the transporters listed above; free filtrate reaches the bladder.
- **Biliary/fecal excretion**: hepatobiliary clearance transports PFOA from hepatocytes into
  bile; bile empties into the intestinal lumen; a fecal clearance term accounts for
  direct intestinal elimination.
- **Enterohepatic recirculation**: PFOA secreted into bile is reabsorbed from the intestinal
  lumen back into the portal circulation.

### Exposure Routes

The model supports five routes of exposure, selectable via `admin.type`:

| `admin.type` | Description |
|---|---|
| `"oral"` | Oral ingestion — dose enters the stomach lumen |
| `"iv"` | Intravenous bolus — dose enters venous blood |
| `"inh"` | Inhalation of dust/aerosol — dose deposited in alveolar lining fluid |
| `"nasal"` | Nasal exposure — dose split between upper airway mucosa and alveolar fluid |
| `"dermal"` | Dermal absorption — dose enters the skin tissue compartment |

Multiple doses at different time points are supported by providing vectors for `admin.dose`
and `admin.time`.

---

## Inputs

| Parameter | Description | Unit / Values |
|---|---|---|
| `BW` | Body weight | kg |
| `sex` | Biological sex (affects blood volume, haematocrit, and urine flow) | `"M"` or `"F"` |
| `admin.dose` | Administered dose(s); provide a vector for multiple doses | µg |
| `admin.time` | Time(s) of dose administration; must match length of `admin.dose` | h |
| `admin.type` | Route of exposure | `"oral"`, `"iv"`, `"inh"`, `"nasal"`, `"dermal"` |
| `depfr_AF` | Deposited fraction of inhaled dose reaching alveolar fluid (inhalation/nasal only) | dimensionless (0–1) |
| `depfr_head` | Deposited fraction of inhaled dose retained in the head/upper airways (nasal only) | dimensionless (0–1) |

---

## Outputs

All concentration variables (`C*`) are in **ng/mL** (equivalent to µg/L).
All mass variables (`M*`) are in **µg**.
Volume variables (`V*`) are in **L**.

The background PFOA plasma level hardcoded in the model is **1.48 ng/mL**, corresponding to
typical population-level serum concentrations; `Cblood_background` converts this to a
whole-blood value.

| Variable | Description | Unit |
|---|---|---|
| `Cblood` | Total PFOA concentration in whole blood (arterial + venous) | ng/mL |
| `Cplasma` | PFOA concentration in plasma, derived from `Cblood` corrected for haematocrit | ng/mL |
| `Cblood_background` | Background whole-blood PFOA concentration (based on 1.48 ng/mL plasma reference) | ng/mL |
| `Ckidney` | Average PFOA concentration across the entire kidney (capillaries + interstitium + cells + tubular filtrate) | ng/mL |
| `Cliver` | Average PFOA concentration across the entire liver (capillaries + interstitium + cells + bile) | ng/mL |
| `Cintestine` | Average PFOA concentration across the entire intestine (capillaries + interstitium + cells + lumen) | ng/mL |
| `Cstomach` | Average PFOA concentration across the entire stomach (capillaries + interstitium + cells + lumen) | ng/mL |
| `Cmuscle` | Average PFOA concentration across the muscle compartment | ng/mL |
| `Cadipose` | Average PFOA concentration across the adipose compartment | ng/mL |
| `Clungs` | Average PFOA concentration across the entire lung (including alveolar lining fluid and retained dust) | ng/mL |
| `Clungtissue` | PFOA concentration in lung tissue only (capillaries + interstitium + cells, excluding free alveolar fluid) | ng/mL |
| `CUpperair` | PFOA concentration in the upper airway (nasal) mucosa fluid | ng/mL |
| `CalveolarLF` | PFOA concentration in alveolar lining fluid (free + dust-associated) | ng/mL |
| `CBALF` | PFOA concentration in bronchoalveolar lavage fluid (BALF); reference volume 5 mL | ng/mL |
| `Cspleen` | Average PFOA concentration in the spleen | ng/mL |
| `Cheart` | Average PFOA concentration in the heart | ng/mL |
| `Cbrain` | Average PFOA concentration in the brain | ng/mL |
| `Cgonads` | Average PFOA concentration in the gonads | ng/mL |
| `Cskin` | Average PFOA concentration in skin | ng/mL |
| `Cbones` | Average PFOA concentration in bone | ng/mL |
| `Crest` | Average PFOA concentration in the lumped "rest" compartment | ng/mL |
| `Ccarcass` | Average PFOA concentration across carcass tissues (muscle + adipose + rest + bone + skin) | ng/mL |
| `CBile` | PFOA concentration in the bile compartment | ng/mL |
| `CBladder` | PFOA concentration in the urinary bladder | ng/mL |
| `Cfeces` | PFOA concentration in feces | ng/mL |
| `Curine` | PFOA concentration in urine | ng/mL |
| `Murine` | Cumulative mass of PFOA excreted via urine | µg |
| `Mfeces` | Cumulative mass of PFOA excreted via feces | µg |
| `Vurine` | Cumulative volume of urine produced | L |
| `Vfeces` | Cumulative volume of feces produced | L |
| `Mblood` | Total PFOA mass in blood (arterial + venous) | µg |
| `Mkidney` | Total PFOA mass in kidney (all sub-compartments + tubular filtrate) | µg |
| `Mliver` | Total PFOA mass in liver (all sub-compartments + bile) | µg |
| `Mbrain` | Total PFOA mass in brain (all sub-compartments) | µg |
