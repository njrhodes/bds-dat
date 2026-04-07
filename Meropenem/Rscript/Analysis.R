# Run this on the first start up

if (!requireNamespace("Pmetrics", quietly = TRUE) ||
    packageVersion("Pmetrics") != "3.0.9") {
  
  repos <- c(
    LAPKB = "https://lapkb.r-universe.dev",
    CRAN  = "https://cloud.r-project.org"
  )
  
  install.packages("Pmetrics", repos = repos)
}

cran_pkgs <- c(
  "tidyverse",
  "dplyr",
  "tibble",
  "kableExtra",
  "plotly",
  "readr",
  "reticulate"
)

missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

# Optional: required for Pmetrics static SVG export via reticulate if desired
if (requireNamespace("reticulate", quietly = TRUE)) {
  try({
    reticulate::use_condaenv("r-reticulate", required = FALSE)
    reticulate::py_run_string("import sys")
  }, silent = TRUE)
}

# -----------------------------------------------------------------------------
# DATA REQUIREMENTS
# -----------------------------------------------------------------------------
# This project expects input data to be placed in:
#   Meropenem/src/
#
# Required files (not included due to DUA restrictions):
#   - MEM1.csv
#   - valid.csv
#   - valid_full.csv
#   - BDS export files (e.g., meropenem_all_subjects_predictions_at_obs_*.csv)
#
# If these files are present and the RProject root is used, the code will run as-is.
# -----------------------------------------------------------------------------

# Package loading ####
library(Pmetrics)
library(tidyverse)
library(dplyr)
library(tibble)
library(kableExtra)
library(plotly)
library(readr)

# Run only once per project
#PM_tree("Meropenem")
stopifnot(dir.exists("Meropenem"))

# Rust re-fit of Rohani et al. Run 15 model ####

dat <- PM_data$new("Meropenem/src/MEM1.csv")

dat$standard_data %>% filter(out!= -99) %>% group_by(outeq) %>% summarise(loq = min(out))

dat <- dat$standard_data %>%
  filter(!id %in% c(1049, 1066, 1086, 1094, 1165, 1294, 1340, 1345, 1354, 1358)) %>%
  arrange(id, time) # remove subjects missing plasma and BAL observations

dat <- PM_data$new(dat, loq=c(0.5,1))

mod1 <- PM_model$new(
  pri = list(
    CL1  = ab(0.1, 20),
    CL2  = ab(0.1, 20),
    CL3  = ab(0.1, 10),
    V    = ab(0.01, 100),
    V2   = ab(0.01, 150),
    K12  = ab(0, 10),
    K21  = ab(0, 20),
    pwr  = ab(0, 2)
  ),
  cov = list(
    age      = interp(),
    male     = interp("none"),
    ht       = interp(),
    wt       = interp(),
    scr      = interp(),
    crcl     = interp(),
    bsa      = interp(),
    crcl_bsa = interp(),
    ecmo     = interp("none"),
    hd       = interp("none"),
    crrt     = interp("none"),
    cvvh     = interp("none"),
    cvvhd    = interp("none"),
    cvvhdf   = interp("none"),
    flow     = interp("none"),
    qb       = interp("none")
  ),
  eqn = function(){
    CLT = CL1 * (crcl/120)**pwr  * (1 - crrt) * (1 - hd) + CL2 * hd + CL3 * crrt * flow/1000
    Ke = CLT / V
    
    dX[1] = R[1] + b[1] - Ke * X[1] - K12 * X[1] + K21 * X[2]
    dX[2] = K12 * X[1] - K21 * X[2]
  },
  out = function(){
    Y[1] = X[1] / V
    Y[2] = X[2] / V2
  },
  err = list(
    proportional(2, c(1, 0.15, 0, 0)),
    proportional(2, c(1, 0.15, 0, 0))
  )
)


#re-fit model in rust for use as a Bayesian prior
run22 <- mod1$fit(data = dat,
                  path = "Meropenem/Runs",
                  cycles = 750,
                  run = 22,
                  points=350 #,
                  #overwrite=T
                  )

run22 <- PM_load(22, path="Meropenem/Runs")

# Apply model to 1 sample Plasma validation cohort ####

dat2 <- PM_data$new("Meropenem/src/valid.csv")

dat2 <- dat2$standard_data

dat2 <- PM_data$new(dat2,loq=c(1,0.5))

print(dat2$summary())

table1_out <- dat2$standard_data %>%
  group_by(id) %>%
  summarise(
    # collapse to first for continuous/demographic
    age  = first(age),
    male = first(male),
    ht   = first(ht),
    wt   = first(wt),
    scr  = first(scr),
    crcl = first(crcl),
    bsa  = first(bsa),
    # take max across subject rows for indicators
    ecmo = max(ecmo, na.rm = TRUE),
    hd   = max(hd,   na.rm = TRUE),
    crrt = max(crrt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  summarise(
    Age = sprintf("%.1f (%.1f), %d–%d",
                  mean(age, na.rm=TRUE),
                  sd(age, na.rm=TRUE),
                  min(age, na.rm=TRUE),
                  max(age, na.rm=TRUE)),
    Male = sprintf("%d (%.1f%%)",
                   sum(male == 1, na.rm=TRUE),
                   mean(male == 1, na.rm=TRUE) * 100),
    Height = sprintf("%.1f (%.1f), %.1f–%.1f",
                     mean(ht, na.rm=TRUE),
                     sd(ht, na.rm=TRUE),
                     min(ht, na.rm=TRUE),
                     max(ht, na.rm=TRUE)),
    Weight = sprintf("%.1f (%.1f), %.1f–%.1f",
                     mean(wt, na.rm=TRUE),
                     sd(wt, na.rm=TRUE),
                     min(wt, na.rm=TRUE),
                     max(wt, na.rm=TRUE)),
    SCr = {
      scr_sub <- filter(., crrt == 0)
      sprintf("%.2f (%.2f), %.2f–%.2f",
              mean(scr_sub$scr, na.rm=TRUE),
              sd(scr_sub$scr, na.rm=TRUE),
              min(scr_sub$scr, na.rm=TRUE),
              max(scr_sub$scr, na.rm=TRUE))
    },
    CrCl = {
      crcl_sub <- filter(., crrt == 0)
      sprintf("%.1f (%.1f), %.1f–%.1f",
              mean(crcl_sub$crcl, na.rm=TRUE),
              sd(crcl_sub$crcl, na.rm=TRUE),
              min(crcl_sub$crcl, na.rm=TRUE),
              max(crcl_sub$crcl, na.rm=TRUE))
    },
    BSA = sprintf("%.2f (%.2f), %.2f–%.2f",
                  mean(bsa, na.rm=TRUE),
                  sd(bsa, na.rm=TRUE),
                  min(bsa, na.rm=TRUE),
                  max(bsa, na.rm=TRUE)),
    ECMO = sprintf("%d (%.1f%%)",
                   sum(ecmo == 1, na.rm=TRUE),
                   mean(ecmo == 1, na.rm=TRUE) * 100),
    HD = sprintf("%d (%.1f%%)",
                 sum(hd == 1, na.rm=TRUE),
                 mean(hd == 1, na.rm=TRUE) * 100),
    CRRT = sprintf("%d (%.1f%%)",
                   sum(crrt == 1, na.rm=TRUE),
                   mean(crrt == 1, na.rm=TRUE) * 100)
  ) %>%
  pivot_longer(
    everything(),
    names_to = "Characteristic",
    values_to = "Summary"
  )

# Render Table 1 for paper ####
library(kableExtra)
table1_out %>%
  kable("html", caption = "Table 1. Baseline Demographics") %>%
  kable_styling(full_width = FALSE, position = "center", bootstrap_options = c("striped", "hover"))

# made sure data compatible
setdiff(names(dat$standard_data),names(dat2$standard_data))

# now run cycles = 0 for the new data
run23 <- mod1$fit(
  data  = dat2,
  path = "Meropenem/Runs",
  cycles = 0,
  prior = 22,
  run   = 23 #,
  #overwrite = TRUE
)

# load fit model 
run23 <- PM_load(23,path="Meropenem/Runs")

# plot the op for validation
run23$op$plot(pred.type="pop",outeq=1)
run23$op$plot(pred.type="post",outeq=1)

run23$op$plot(pred.type="pop",outeq=1,resid=T)
run23$op$plot(pred.type="post",outeq=1,resid=T)

# plot the validation results
library(plotly)
f <- list(size = 20)

a <- list(
  text = "A. Population",
  font=f,
  xref = "paper",
  yref = "paper",
  yanchor = "bottom",
  xanchor = "center",
  align = "right",
  x = 0.075,
  y = 1,
  showarrow = FALSE
)

b <- list(
  text = "B. Posterior",
  font=f,
  xref = "paper",
  yref = "paper",
  yanchor = "bottom",
  xanchor = "center",
  align = "right",
  x = 0.05,
  y = 1,
  showarrow = FALSE
)

vp1<-run23$op$plot(pred.type="pop",outeq=1,
                   stats = list(x = 0.75, y = 0.25, font = list(size = 10)))
vp2<-run23$op$plot(pred.type="post",outeq=1,
                   stats = list(x = 0.95, y = 0.25, font = list(size = 10)))

vp1 <- vp1 %>% layout(annotations = a)
vp2 <- vp2 %>% layout(annotations = b)

# Export Figure 1 plot
can_export_svg <- FALSE

if (requireNamespace("reticulate", quietly = TRUE)) {
  try({
    reticulate::use_condaenv("r-reticulate", required = FALSE)
    reticulate::py_run_string("import kaleido")
    can_export_svg <- TRUE
  }, silent = TRUE)
}

if (can_export_svg) {
  sub_plot(vp1, vp2, nrows = 1, margin=0.05, titleX=T, titleY=T,shareX=T,shareY=T) %>%
  export_plotly("MEM PLA EXT VALID Plot.svg", width = 3.5 * 300, height = 1.5 * 300)
} else {
  message("SVG export skipped (kaleido not available)")
}

# Generate Pmetrics predictions for external sample ####

crrt_lookup <- run23$data$data %>%
  group_by(id) %>%
  summarise(
    crrt_status = as.integer(any(crrt == 1, na.rm = TRUE)),
    .groups = "drop"
  )

subj_level <- run23$op$data %>%
  filter(outeq == 1, icen == "median", obs != 0) %>%
  left_join(crrt_lookup, by = "id") %>%
  group_by(id, pred.type, crrt_status) %>%
  summarise(
    subj_pe   = median((pred-obs) / obs, na.rm = TRUE),
    subj_ape  = median(abs( pred-obs) / obs, na.rm = TRUE),
    .groups = "drop"
  )

np_posterior <- subj_level %>%
  mutate(CRRT = factor(crrt_status, levels = c(1, 0), labels = c("CRRT", "No CRRT"))) %>%
  group_by(CRRT, pred.type) %>%
  summarise(
    mpe  = round(100 * median(subj_pe, na.rm = TRUE), 1),
    mape = round(100 * median(subj_ape, na.rm = TRUE), 1),
    n    = n_distinct(id),
    .groups = "drop"
  )

np_posterior


# Generate MCS from Pmetrics using first sample TDM and run 22 prior

# simulated median obs vs pred
sim1 <- run22$sim(
  data=dat2,
  limits=NA,
  nsim=1000
)

# 1) Stack simulation outputs and keep outeq==1
sim_df <- sim1$obs %>%
  filter(outeq == 1)

# 2) Median simulated prediction per simulation grid (one time per list element)
sim_medians <- sim_df %>%
  group_by(time) %>%
  summarise(
    pred_median = median(out, na.rm = TRUE),
    .groups = "drop"
  )

grid_unique <- sim_medians

# 3) Observation times (outeq==1) from the standard data
obs_tbl <- dat2$standard_data %>%
  filter(outeq == 1, !is.na(time)) %>%
  transmute(
    id_subject = id,
    time_obs   = time,
    out_obs    = out
  )

# 4) Nearest-time match using base R (no data.table)
closest_index <- function(x, vec_sorted) {
  # vec_sorted must be sorted ascending
  idx <- findInterval(x, vec_sorted)              # in [0, length-1]
  idx[idx == 0] <- 1                              # clamp to left boundary
  right <- pmin(idx + 1, length(vec_sorted))      # right neighbor
  left_time  <- vec_sorted[idx]
  right_time <- vec_sorted[right]
  pick_right <- abs(right_time - x) < abs(x - left_time)
  idx[pick_right] <- right[pick_right]
  idx
}

tgrid <- grid_unique$time
# ensure sorted (it already is from arrange(), but be explicit)
o <- order(tgrid)
tgrid <- tgrid[o]
grid_unique <- grid_unique[o, ]

match_idx <- closest_index(obs_tbl$time_obs, tgrid)

obs_matched <- obs_tbl %>%
  mutate(
    time_sim = tgrid[match_idx],
    dt_abs   = abs(time_sim - time_obs)
  )

# 5) Bring in median prediction at the matched sim time
sim_result <- obs_matched %>%
  left_join(
    grid_unique %>% select(time, pred_median),
    by = c("time_sim" = "time")
  ) %>%
  arrange(id_subject, time_obs) %>%
  left_join(crrt_lookup, by = c("id_subject" = "id")) %>%
      mutate(
        CRRT = factor(crrt_status, levels = c(1, 0), labels = c("CRRT", "No CRRT"))
      ) %>%
  group_by(CRRT, .drop = FALSE) %>%
  summarise(
    mpe  = 100 * median((pred_median - out_obs) / out_obs, na.rm = TRUE),
    mape = 100 * median(abs(pred_median - out_obs) / out_obs, na.rm = TRUE),
    n    = dplyr::n(),
    .groups = "drop"
  )

sim_result

# Table 2 out
table2_out <- bind_rows(
  np_posterior %>%
    mutate(
      Source = "MAP Bayesian (Pmetrics)",
      `Prediction type` = dplyr::case_when(
        pred.type == "post" ~ "Individual",
        pred.type == "pop"  ~ "Population",
        TRUE ~ pred.type
      )
    ) %>%
    transmute(
      Source,
      `CRRT status` = as.character(CRRT),
      `Prediction type`,
      `rMPE (%)` = as.numeric(mpe),
      `rMAPE (%)` = as.numeric(mape),
      N = n
    ),
  sim_result %>%
    mutate(
      Source = "Simulation (Pmetrics)",
      `Prediction type` = "Population"
    ) %>%
    transmute(
      Source,
      `CRRT status` = dplyr::recode(as.character(CRRT), `0` = "No CRRT", `1` = "CRRT"),
      `Prediction type`,
      `rMPE (%)` = as.numeric(mpe),
      `rMAPE (%)` = as.numeric(mape),
      N = n
    )
) %>%
  mutate(
    `CRRT status` = dplyr::recode(`CRRT status`, `0` = "No CRRT", `1` = "CRRT")
  ) %>%
  arrange(Source, `CRRT status`, `Prediction type`)

table2_out %>%
  mutate(
    `rMPE (%)` = sprintf("%.1f", `rMPE (%)`),
    `rMAPE (%)` = sprintf("%.1f", `rMAPE (%)`)
  ) %>%
  kable(
    caption = "Table 2. Prediction error in the one-sample-per-subject validation dataset, by platform, CRRT status, and prediction type",
    align = "lllrcc",
    col.names = c("Source", "CRRT status", "Prediction type", "rMPE (%)", "rMAPE (%)", "N")
  ) %>%
  kable_styling(full_width = FALSE, position = "center", bootstrap_options = c("striped", "hover"))

# Table 3 out from BDS Result: sourcing from app_mod_up_res.R file and stan: state2
summarize_bds <- function(path, scenario_label) {
  df <- read_csv(path, show_col_types = FALSE)
  
  if ("pred_conc" %in% names(df) && !"IPRED" %in% names(df)) df$IPRED <- df$pred_conc
  if ("obs_conc"  %in% names(df) && !"DV"    %in% names(df)) df$DV    <- df$obs_conc
  
  df %>%
    transmute(
      ID = as.character(ID),
      DV = as.numeric(DV),
      IPRED = as.numeric(IPRED)
    ) %>%
    filter(!is.na(ID), !is.na(DV), !is.na(IPRED), DV != 0) %>%
    summarise(
      `Sample scenario` = scenario_label,
      `N (patients)`    = n_distinct(ID),
      `N (samples)`     = n(),
      `rMPE (%)`        = round(100 * median((IPRED - DV) / DV, na.rm = TRUE), 1),
      `rMAPE (%)`       = round(100 * median(abs(IPRED - DV) / DV, na.rm = TRUE), 1),
      `F20 (%)`         = round(100 * mean(abs(IPRED - DV) / DV <= 0.20, na.rm = TRUE), 0),
      `F30 (%)`         = round(100 * mean(abs(IPRED - DV) / DV <= 0.30, na.rm = TRUE), 0)
    )
}

table3 <- bind_rows(
  summarize_bds("Meropenem/src/meropenem_all_subjects_predictions_at_obs_1.csv", "One observation per patient"),
  summarize_bds("Meropenem/src/meropenem_all_subjects_predictions_at_obs_2.csv", "Two observations per patient"),
  summarize_bds("Meropenem/src/meropenem_all_subjects_predictions_at_obs_3.csv", "Three observations per patient")
)

table3 %>%
  kable(
    caption = "Table 3. Predictive Performance of Bayesian Dosing Simulator (BDS) Across Varied Sampling Scenarios",
    align = c("l", "r", "r", "r", "r", "r", "r")
  ) %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))


# ---------------------------------------------------------
# Figure S7: BDS WRES
# uses one-sample-per-subject export from app
# expects columns: ID, obs_conc, IPRED, PRED, SD, IWRES
# ---------------------------------------------------------
bds <- read_csv("Meropenem/src/meropenem_all_subjects_predictions_at_obs_1.csv",
                show_col_types = FALSE)

if (!all(c("PRED", "PWRES") %in% names(bds))) {
  stop("BDS export must include PRED and PWRES (rerun app_mod_up_res.R)")
}

bds_long <- bind_rows(
  bds %>%
    transmute(
      ID = ID,
      model = "Post",
      Prediction = IPRED,
      WRES = IWRES
    ),
  bds %>%
    transmute(
      ID = ID,
      model = "Pop",
      Prediction = PRED,
      WRES = PWRES
    )
)

fig_s7 <- ggplot(bds_long, aes(x = Prediction, y = WRES)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", linewidth = 0.4, color = "firebrick") +
  geom_point(size = 3) +
  facet_wrap(~ model, scales = "free_x") +
  theme_bw() +
  labs(
    title = "Figure S7. Weighted residuals (WRES) for BDS predictions",
    x = "Predicted Concentration",
    y = "WRES"
  )

fig_s7

# optional ID audits
bds_ids <- bds_long %>%
  distinct(ID, model) %>%
  arrange(model, ID)

bds_outlier_ids <- bds_long %>%
  filter(abs(WRES) > 2) %>%
  distinct(ID, model, WRES) %>%
  arrange(model, desc(abs(WRES)))

bds_outlier_ids

# ---------------------------------------------------------
# Figure S8: Pmetrics WRES
# uses run23 validation object
# ---------------------------------------------------------

pm_long <- run23$op$data %>%
  filter(
    icen == "median",
    outeq == 1,
    pred.type %in% c("pop", "post"),
    !is.na(pred),
    !is.na(wd)
  ) %>%
  transmute(
    ID = as.character(id),
    model = recode(pred.type, pop = "Pop", post = "Post"),
    Prediction = as.numeric(pred),
    WRES = as.numeric(wd)
  )

fig_s8 <- ggplot(pm_long, aes(x = Prediction, y = WRES)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", linewidth = 0.4, color = "firebrick") +
  geom_point(size = 3) +
  facet_wrap(~ model, scales = "free_x") +
  theme_bw() +
  labs(
    title = "Figure S8. Weighted residuals (WRES) for Pmetrics predictions",
    x = "Predicted Concentration",
    y = "WRES"
  )

fig_s8

# optional ID audits
pm_ids <- pm_long %>%
  distinct(ID, model) %>%
  arrange(model, ID)

pm_outlier_ids <- pm_long %>%
  filter(abs(WRES) > 2) %>%
  distinct(ID, model, WRES) %>%
  arrange(model, desc(abs(WRES)))

pm_outlier_ids
