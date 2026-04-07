# Meropenem Bayesian Validation Analysis

This repository contains code to reproduce the analysis, tables, and figures for external validation of a meropenem population PK model using Pmetrics and a Bayesian Dosing Simulator (BDS).

---

## Repository Structure
```
App/
  app_mod_up_res.R
  models/
    mem_plasma_model_state2.stan
    mem_plasma_model_state2.exe

Meropenem/
  Rscript/
    Analysis.R
  Runs/        # Pmetrics outputs (NOT tracked in git)
    22/
    23/
    MEM PLA EXT VALID Plot.svg
  Sim/
  src/         # input + BDS export CSVs (NOT tracked in git)
    MEM1.csv
    valid.csv
    valid_full.csv
    meropenem_all_subjects_predictions_at_obs_1.csv
    meropenem_all_subjects_predictions_at_obs_2.csv
    meropenem_all_subjects_predictions_at_obs_3.csv

bds-dat.Rproj
prior.csv
```
---

## Overview

Workflow:

1. Prior model fit (run22) using Pmetrics  
2. Bayesian update (run23, cycles = 0)  
3. Simulation benchmarking  
4. BDS predictions via Shiny app (Stan state2 model)  
5. Generation of tables and figures  

Outputs:

- Table 1: Baseline demographics  
- Table 2: Prediction error (Pmetrics + simulation)  
- Table 3: BDS performance  
- Figure S7: BDS WRES  
- Figure S8: Pmetrics WRES  

---

## Data Requirements

Data are not included due to DUA restrictions.

Place the following files in:

Meropenem/src/

Required:

- MEM1.csv  
- valid.csv  
- valid_full.csv  

BDS exports (generated via app):

- meropenem_all_subjects_predictions_at_obs_1.csv  
- meropenem_all_subjects_predictions_at_obs_2.csv  
- meropenem_all_subjects_predictions_at_obs_3.csv  

---

## Running the Analysis

Open:

bds-dat.Rproj

Run:

Meropenem/Rscript/Analysis.R

This will:

- fit models (run22, run23)
- run simulations
- compute tables and figures

---

## BDS Workflow

Run:

App/app_mod_up_res.R

This uses:

App/models/mem_plasma_model_state2.stan

Required output columns:

ID, time, obs_conc, IPRED, PRED, SD, IWRES, PRED_SD, PWRES

---

## Prediction Error Definition

All analyses use:

(P - O) / O

rMPE = median relative prediction error  
rMAPE = median absolute relative prediction error  

---

## Residuals

Pmetrics:
- wd

BDS:
- Posterior: IWRES  
- Population: PWRES  

---

## Dependencies

R packages:

Pmetrics (v3.0.9 required)
tidyverse
dplyr
tibble
kableExtra
plotly
readr
reticulate (optional)

---

## Pmetrics Installation (Version Locked)

This analysis requires Pmetrics version 3.0.9.

The analysis script installs it automatically from LAPKB r-universe if needed.

---

## CmdStan / BDS App Notes

The Shiny app:

- will install CmdStan automatically if not present  
- will compile the Stan model on startup  

This may take several minutes on first run.

---

## SVG Export (Optional)

Static SVG export requires Python + kaleido.

If not configured:
- analysis runs normally  
- plots render  
- SVG export is skipped  

---

## Reproducibility Notes

- Pmetrics fitting and simulation may vary slightly across machines due to search behavior  
- Results correspond to a specific execution environment  
- Code structure is deterministic but not bitwise reproducible across systems  

---

## Notes

- Data are excluded due to DUA  
- Repository is code-only  
- Requires local data and BDS exports to run  

