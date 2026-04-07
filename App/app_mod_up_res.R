# app.R  (FULLY WIRED: IWRES extraction + IWRES plots for manual + upload subject viewer
#         + exports include IWRES/IPRED/SD)

library(shiny)
library(cmdstanr)
library(ggplot2)
library(plotly)
library(shinycssloaders)
library(DT)
library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(zoo)    # for na.locf

if (is.null(cmdstanr::cmdstan_path())) {
  cmdstanr::install_cmdstan(cores = parallel::detectCores())
}

stan_path <- "models/mem_plasma_model_state2.stan"
mod <- cmdstan_model(
  stan_path,
  compile = FALSE,
  cpp_options = list(CXXFLAGS = "-Wno-overloaded-virtual")
)
mod$compile()

# ---------- Helpers ---------------------------------------------------------

cv2sd <- function(cv) sqrt(log(1 + (cv/100)^2))
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Build covariate state-change table for ONE subject
# Uses baseline + rows with AMT=0 & TINF=0; if HD column missing, toggles HD on each state-change
build_cov_changes <- function(df_subject) {
  # Detect whether an HD column is present/non-NA
  hd_present <- "HD" %in% names(df_subject) && any(!is.na(df_subject$HD))
  
  base_row <- df_subject %>%
    arrange(TIME) %>% slice_head(n = 1) %>%
    mutate(HD = if (hd_present) dplyr::coalesce(as.integer(HD), 0L) else 0L) %>%
    select(TIME, CRCL, CRRT, FLOW, HD)
  
  sc_rows <- df_subject %>%
    filter(AMT == 0, TINF == 0) %>%   # state changes
    arrange(TIME) %>%
    mutate(HD = if (hd_present) dplyr::coalesce(as.integer(HD), 0L) else NA_integer_) %>%
    select(TIME, CRCL, CRRT, FLOW, HD)
  
  cov_tbl <- bind_rows(base_row, sc_rows) %>%
    arrange(TIME)
  
  # If HD missing, toggle HD at each state-change time
  if (!hd_present) {
    hd <- integer(nrow(cov_tbl))
    if (nrow(cov_tbl) > 1) {
      for (k in 2:nrow(cov_tbl)) hd[k] <- 1L - hd[k-1]
    }
    cov_tbl$HD <- hd
  }
  
  cov_tbl <- cov_tbl %>%
    mutate(
      CRCL = na.locf(CRCL, na.rm = FALSE),
      CRRT = na.locf(CRRT, na.rm = FALSE),
      FLOW = na.locf(FLOW, na.rm = FALSE),
      HD   = na.locf(HD,   na.rm = FALSE)
    ) %>%
    tidyr::fill(CRCL, CRRT, FLOW, HD, .direction = "down")
  
  cov_tbl$CRCL[is.na(cov_tbl$CRCL)] <- 60
  cov_tbl$CRRT[is.na(cov_tbl$CRRT)] <- 0
  cov_tbl$FLOW[is.na(cov_tbl$FLOW)] <- 0
  cov_tbl$HD[is.na(cov_tbl$HD)]     <- 0
  cov_tbl$CRRT <- as.integer(cov_tbl$CRRT)
  cov_tbl$HD   <- as.integer(cov_tbl$HD)
  cov_tbl
}

# Evaluate covariates at requested times using piecewise-constant step function
eval_cov_at_times <- function(times, cov_tbl) {
  idx <- findInterval(times, cov_tbl$TIME, left.open = FALSE, rightmost.closed = TRUE)
  idx[idx < 1] <- 1L
  data.frame(
    CRCL_i = as.numeric(cov_tbl$CRCL[idx]),
    CRRT_i = as.integer(cov_tbl$CRRT[idx]),
    FLOW_i = as.numeric(cov_tbl$FLOW[idx]),
    HD_i   = as.integer(cov_tbl$HD[idx])
  )
}

RES_ADD_SIGMA  <- 0.5   # mg/L  (set to your validated population residual)
RES_PROP_SIGMA <- 0.07  # proportional

# Core runner used by both manual & upload paths
# Supports time-varying covariates (incl. HD) via cov_tbl/cov_at_obs; bolus (TINF==0) allowed
run_map_meropenem <- function(
    obs_df,           # cols: time, conc
    dosing_df,        # cols: dose, tinf, start
    cov_tbl,          # data.frame(TIME, CRCL, CRRT, FLOW, HD)
    cov_at_obs,       # data.frame(CRCL_i, CRRT_i, FLOW_i, HD_i) aligned to obs_df$time
    fu,
    priors, mic, tau,
    time_horizon = NULL
) {
  stopifnot(all(c("time","conc") %in% names(obs_df)))
  stopifnot(all(c("dose","tinf","start") %in% names(dosing_df)))
  stopifnot(all(c("CRCL_i","CRRT_i","FLOW_i","HD_i") %in% names(cov_at_obs)))
  
  od <- obs_df %>% arrange(time) %>% mutate(time=as.numeric(time), conc=as.numeric(conc))
  dz <- dosing_df %>% arrange(start)
  
  if (!nrow(od)) stop("No observations.")
  if (!nrow(dz)) stop("No dosing rows.")
  if (any(!is.finite(od$time)) || any(od$time < 0)) stop("Obs times must be non-negative/finite.")
  if (any(!is.finite(od$conc)) || any(od$conc < 0)) stop("Obs concentrations must be >=0/finite.")
  if (any(!is.finite(dz$start)) || any(dz$start < 0)) stop("Dose start times must be non-negative/finite.")
  # Allow bolus (0), forbid negatives
  if (any(!is.finite(dz$tinf))) stop("Infusion durations must be finite (use 0 for bolus).")
  if (any(dz$tinf < 0)) stop("Infusion durations must be >= 0 (0 = bolus).")
  if (fu <= 0 || fu > 1) stop("fu must be between 0 and 1.")
  
  # Stan data (per-observation covariates)
  dat <- list(
    N             = nrow(od),
    time          = od$time,
    conc          = od$conc,
    ndoses        = nrow(dz),
    dose_vec      = as.array(dz$dose),
    infusion_vec  = as.array(dz$tinf),
    start_vec     = as.array(dz$start),
    
    CRCL_i        = as.array(cov_at_obs$CRCL_i),
    CRRT_i        = as.integer(cov_at_obs$CRRT_i),
    FLOW_i        = as.array(cov_at_obs$FLOW_i),
    HD_i          = as.integer(cov_at_obs$HD_i),
    
    V_pop         = priors$V_pop$mean,   V_sd   = priors$V_pop$sd,
    CL1_pop       = priors$CL1_pop$mean, CL1_sd = priors$CL1_pop$sd,
    CL2_pop       = priors$CL2_pop$mean, CL2_sd = priors$CL2_pop$sd,
    pwr_pop       = priors$pwr_pop$mean, pwr_sd = priors$pwr_pop$sd,
    CLHD_pop      = priors$CLHD_pop$mean, CLHD_sd = priors$CLHD_pop$sd,
    
    add_sigma     = RES_ADD_SIGMA,
    prop_sigma    = RES_PROP_SIGMA
  )
  
  fit <- mod$optimize(data = dat, seed = 123)
  
  summ <- fit$summary()
  get_map <- function(v){
    r <- summ[summ$variable == v, , drop=FALSE]
    if (nrow(r) != 1) stop(v," missing")
    if ("mean" %in% names(r)) r$mean else r$estimate
  }
  
  eta_V    <- get_map("eta_V")
  eta_CL1  <- get_map("eta_CL1")
  eta_CL2  <- get_map("eta_CL2")
  eta_pwr  <- get_map("eta_pwr")
  eta_CLHD <- get_map("eta_CLHD")
  
  # ---- residual SD parameters (from your Stan) ----
  add_sigma  <- RES_ADD_SIGMA
  prop_sigma <- RES_PROP_SIGMA
  
  V_map    <- priors$V_pop$mean    * exp(eta_V)
  CL1_map  <- priors$CL1_pop$mean  * exp(eta_CL1)
  CL2_map  <- priors$CL2_pop$mean  * exp(eta_CL2)
  pwr_map  <- priors$pwr_pop$mean  * exp(eta_pwr)
  CLHD_map <- priors$CLHD_pop$mean * exp(eta_CLHD)
  
  # Simulation grid ---------------------------------------------------------
  end_times <- if (nrow(dz)) dz$start + dz$tinf else numeric(0)
  last_obs  <- if (nrow(od)) max(od$time, na.rm = TRUE) else 0
  base_h    <- max(c(end_times, last_obs, 48), na.rm = TRUE)
  
  if (!is.null(time_horizon) && is.finite(time_horizon)) {
    base_h <- max(base_h, time_horizon)
  }
  horiz_max <- base_h
  
  by_step <- 0.1
  times <- seq(0, horiz_max, by = by_step)
  
  kinks <- c(dz$start, dz$start + dz$tinf, od$time)
  kinks <- kinks[is.finite(kinks) & kinks >= 0]
  if (length(kinks)) times <- sort(unique(c(times, kinks)))
  
  if (nrow(od) && max(od$time, na.rm = TRUE) > max(times) + 1e-9) {
    stop("Internal: horizon does not include last observation time.")
  }
  
  cov_grid <- eval_cov_at_times(times, cov_tbl)
  
  sim_fn_tv <- function(V, CL1, CL2, pwr, CLHD, fu, dose_df) {
    sapply(seq_along(times), function(j) {
      t <- times[j]
      CRCLt <- cov_grid$CRCL_i[j]
      CRRTt <- cov_grid$CRRT_i[j]
      FLOWt <- cov_grid$FLOW_i[j]
      HDt   <- cov_grid$HD_i[j]
      
      sum(mapply(function(dose_i, t0, inf_t) {
        dt <- t - t0
        if (dt <= 0) return(0)
        
        CLT_t <- max(CL1*(CRCLt/120)^pwr*(1-CRRTt) +
                       CL2*CRRTt*(FLOWt/1000) +
                       CLHD*HDt, 1e-6)
        kel   <- CLT_t / V
        
        if (inf_t > 0) {
          # infusion
          rate  <- fu * (dose_i/inf_t) / CLT_t
          if (dt <= inf_t) rate * (1 - exp(-kel*dt))
          else             rate * (1 - exp(-kel*inf_t)) * exp(-kel*(dt - inf_t))
        } else {
          # bolus
          fu * (dose_i / V) * exp(-kel * dt)
        }
      }, dose_df$dose, dose_df$start, dose_df$tinf, SIMPLIFY = TRUE))
    })
  }
  
  indiv <- sim_fn_tv(V_map, CL1_map, CL2_map, pwr_map, CLHD_map, fu, dz)
  pop   <- sim_fn_tv(priors$V_pop$mean, priors$CL1_pop$mean, priors$CL2_pop$mean,
                     priors$pwr_pop$mean, priors$CLHD_pop$mean, fu, dz)
  
  # Diagnostics (align to obs times)
  pick_idx <- sapply(od$time, function(t) which.min(abs(times - t)))
  pred_at_obs <- indiv[pick_idx]
  pop_at_obs  <- pop[pick_idx]
  
  pe_vec <- (pred_at_obs - od$conc) / ((pred_at_obs + od$conc) / 2)
  abs_pe_vec <- abs((pred_at_obs - od$conc)) / ((pred_at_obs + od$conc) / 2)
  mpe    <- median(pe_vec) * 100
  mape   <- median(abs_pe_vec) * 100
  f0_24  <- mean(indiv[times<=24] >= mic)
  f24_48 <- if (max(times)>24) mean(indiv[times>24 & times<=48] >= mic) else NA_real_
  f_tau  <- mean(indiv[times<=tau] >= mic)
  
  # Compute CLT and Ke at t0 for reporting (informative only)
  CLT_t0 <- max(CL1_map*(cov_grid$CRCL_i[1]/120)^pwr_map*(1-cov_grid$CRRT_i[1]) +
                  CL2_map*cov_grid$CRRT_i[1]*(cov_grid$FLOW_i[1]/1000) +
                  CLHD_map*cov_grid$HD_i[1], 1e-6)
  Kel_t0 <- CLT_t0 / V_map
  
  # ---- IWRES (matches Stan likelihood) ----
  sd_post  <- sqrt(RES_ADD_SIGMA^2 + (pmax(pred_at_obs, 0) * RES_PROP_SIGMA)^2)
  iwres_post <- (pred_at_obs - od$conc ) / pmax(sd_post, 1e-9)
  
  sd_pop   <- sqrt(RES_ADD_SIGMA^2 + (pmax(pop_at_obs, 0) * RES_PROP_SIGMA)^2)
  iwres_pop  <- (pop_at_obs - od$conc ) / pmax(sd_pop, 1e-9)
  
  resid_df <- od %>%
    mutate(
      IPRED = pred_at_obs,
      PRED  = pop_at_obs,
      SD    = sd_post,
      IWRES = iwres_post,
      PRED_SD = sd_pop,
      PWRES   = iwres_pop
    )
  
  list(
    map_params = data.frame(
      Parameter = c("V (L)","CLT (L/hr)","Ke (1/hr)","MPE (%)","MAPE (%)"),
      Value     = c(V_map, CLT_t0, Kel_t0, mpe, mape)
    ),
    ftmic = data.frame(
      Metric = c("fT>MIC 0–24h","fT>MIC 24–48h","fT>MIC tau"),
      Value  = c(f0_24*100, ifelse(is.na(f24_48), NA, f24_48*100), f_tau*100)
    ),
    times = times,
    indiv = indiv,
    pop   = pop,
    obs   = od,
    pred_at_obs = pred_at_obs,
    pop_at_obs  = pop_at_obs,
    
    # residual outputs
    resid_df = resid_df,
    add_sigma = add_sigma,
    prop_sigma = prop_sigma
  )
}

# ---------- UI --------------------------------------------------------------

ui <- fluidPage(
  titlePanel("MAP Bayesian PK Estimation for Meropenem"),
  tabsetPanel(
    id = "mainTabs",
    
    tabPanel("Inputs",
             fluidRow(
               column(4,
                      h4("Dosing Inputs"),
                      wellPanel(
                        numericInput("dose",          "Dose (mg):",               1000),
                        numericInput("infusion_time", "Infusion duration (hr):",    3),
                        numericInput("tau",           "Dose Interval (hr):",        8),
                        numericInput("n_doses",       "Number of Doses:",           6),
                        numericInput("mic",           "Target MIC (mg/L):",         1),
                        numericInput("fu",            "Fraction unbound (0–1):",  0.98)
                      )
               ),
               column(4,
                      h4("Renal Inputs"),
                      wellPanel(
                        numericInput("crcl", "CrCL (mL/min):", 60),
                        radioButtons("crrt","CRRT:", choices=c("No"=0,"Yes"=1), selected=0),
                        numericInput("flow","CRRT Flow (mL/min):", 0),
                        checkboxInput("use_load", "Include loading dose?", FALSE),
                        conditionalPanel(
                          "input.use_load",
                          selectInput("ld_dose", "Loading dose (mg):", choices=c(1000,2000)),
                          helpText("Over 0.5 hr")
                        )
                      ),
                      actionButton("simulate","Simulate",class="btn-primary")
               ),
               column(4,
                      h4("Observed Concentrations"),
                      uiOutput("obs_inputs"),
                      actionButton("add_obs",    "+ Obs",    class="btn-sm btn-success"),
                      actionButton("remove_obs", "– Obs",    class="btn-sm btn-danger"),
                      br(), br(),
                      actionButton("toPriors", "Next: Priors", class="btn-secondary")
               )
             )
    ),
    
    tabPanel("Population Priors",
             fluidRow(
               column(12,
                      h4("Population Priors (Log-Normal)"),
                      helpText("SD from CV% assuming log-normal IIV"),
                      DTOutput("prior_table_ui"),
                      actionButton("toPlot","Next: Plot",class="btn-secondary")
               )
             )
    ),
    
    # UPDATED: Plot tab now includes IWRES subtab
    tabPanel("Plot",
             tabsetPanel(
               tabPanel("Concentration", withSpinner(plotlyOutput("concPlot"))),
               tabPanel("IWRES",         withSpinner(plotlyOutput("iwresPlot")))
             ),
             br(),
             actionButton("toParams","Next: Parameters",class="btn-secondary")
    ),
    
    tabPanel("Parameters",
             tableOutput("paramTable"), br(),
             actionButton("toSummary","Next: Summary",class="btn-secondary")
    ),
    
    tabPanel("Summary",
             tableOutput("ftmicTable")
    ),
    
    tabPanel("Upload & Run",
             fluidRow(
               column(4,
                      h4("Standardized CSV"),
                      fileInput("datafile", "Upload CSV", accept = c(".csv")),
                      helpText(HTML(
                        "Required columns (case-insensitive): ",
                        "<code>ID</code>, <code>EVID</code> (1=dose, 0=obs), ",
                        "<code>AMT</code> (mg), <code>TINF</code> (hr), <code>TIME</code> (hr), ",
                        "<code>DV</code> (mg/L), covariates: <code>CRCL</code>, <code>CRRT</code> (0/1), <code>FLOW</code> (mL/min). ",
                        "Rows with <code>AMT==0</code> & <code>TINF==0</code> act as covariate (and HD) state-changes."
                      )),
                      checkboxInput("upload_fu_override", "Override fu with Inputs tab value", TRUE),
                      actionButton("parse_btn","Parse file"),
                      hr(),
                      h5("Level selection"),
                      radioButtons("level_cap", "Use per subject:",
                                   choices = c("First level only" = "cap1",
                                               "Up to 2 levels"   = "cap2",
                                               "Up to 3 levels"   = "cap3",
                                               "All available samples" = "all"),
                                   selected = "cap1"
                      ),
                      radioButtons("min_levels", "Include subjects with at least:",
                                   choices = c("1 level" = 1, "2 levels" = 2, "3 levels" = 3), selected = 1
                      ),
                      actionButton("preview_cohort","Preview eligible subjects", class="btn-light"),
                      br(), br(),
                      actionButton("run_all_btn","Run All Subjects", class="btn-primary")
               ),
               column(8,
                      h4("Preview (raw & cohort)"),
                      tabsetPanel(id = "uploadTabs",
                                  tabPanel("Subject Viewer",
                                           uiOutput("subject_picker"),
                                           actionButton("run_subject_btn","Run Selected Subject", class="btn-secondary"),
                                           br(), br(),
                                           withSpinner(plotlyOutput("upload_concPlot")),
                                           br(),
                                           withSpinner(plotlyOutput("upload_iwresPlot")),  # NEW
                                           br(),
                                           tableOutput("upload_paramTable"),
                                           tableOutput("upload_ftmicTable")
                                  ),
                                  tabPanel("Cohort Preview",
                                           DTOutput("upload_preview"),
                                           br(),
                                           DTOutput("cohort_summary"),
                                           br(),
                                           actionButton("back_to_subject","Back to Subject Viewer", class="btn-light")
                                  )
                      ),
                      br(),
                      actionButton("toReports","Next: Reports", class="btn-secondary")
               )
             )
    ),
    
    tabPanel("Reports",
             fluidRow(
               column(6,
                      h4("Subject report (single)"),
                      uiOutput("report_subject_picker"),
                      downloadButton("download_subject_metrics", "Download Subject Metrics (CSV)"),
                      downloadButton("download_subject_predictions", "Download Subject Predictions (CSV)"),
                      br(), br(),
                      h5("Most recent subject run"),
                      tableOutput("recent_metrics")
               ),
               column(6,
                      h4("Global cohort results"),
                      withSpinner(DTOutput("global_metrics_dt")),
                      br(),
                      withSpinner(DTOutput("run_log_dt")),
                      br(),
                      downloadButton("download_all_metrics_csv", "Download All Metrics (CSV)"),
                      downloadButton("download_all_predictions_csv", "Download All Predictions at Obs Times (CSV)")
               )
             )
    )
  )
)

# ---------- Server ----------------------------------------------------------

server <- function(input, output, session) {
  
  # Priors table -------------------------------------------------------------
  # Include CLHD prior, but hide it from the Priors UI
  prior_vals <- reactiveVal(
    data.frame(
      StanName  = c("V_pop","CL1_pop","CL2_pop","pwr_pop","CLHD_pop"),
      Parameter = c("V (L)","CL1 (L/hr)","CL2 (L/hr)","Power","CL_HD during HD (L/hr)"),
      Mean      = c(29.7,   8.5,        2.2,         0.75,    9.5),
      CV        = c(75,     52,         100,         80,      75),
      stringsAsFactors = FALSE
    )
  )
  
  output$prior_table_ui <- renderDT({
    df <- prior_vals()
    df$SD <- round(cv2sd(df$CV), 3)
    df_disp <- df[df$StanName != "CLHD_pop", , drop = FALSE]  # hide CLHD
    datatable(
      df_disp[,c("Parameter","Mean","SD","CV")],
      editable = list(target="cell", columns=3),
      rownames = FALSE,
      options = list(dom='t')
    )
  })
  
  observeEvent(input$prior_table_ui_cell_edit, {
    info <- input$prior_table_ui_cell_edit
    if (info$col == 3) { # CV column in display
      df <- prior_vals()
      df_disp <- df[df$StanName != "CLHD_pop", , drop = FALSE]
      stan_name <- df_disp$StanName[info$row]
      idx <- which(df$StanName == stan_name)
      df$CV[idx] <- as.numeric(info$value)
      prior_vals(df)
    }
  })
  
  get_priors_list <- reactive({
    pv <- prior_vals()
    setNames(lapply(seq_len(nrow(pv)), function(i) list(mean=pv$Mean[i], sd=cv2sd(pv$CV[i]))),
             pv$StanName)
  })
  
  # Manual Inputs path -------------------------------------------------------
  rv <- reactiveValues(n_obs = 2)
  observeEvent(input$add_obs,    { rv$n_obs <- min(rv$n_obs+1, 10) })
  observeEvent(input$remove_obs, { rv$n_obs <- max(rv$n_obs-1, 1)  })
  
  output$obs_inputs <- renderUI({
    tagList(lapply(seq_len(rv$n_obs), function(i){
      fluidRow(
        column(6, numericInput(paste0("obs",i),  paste0("Obs Conc ",i," (mg/L):"), value = NA)),
        column(6, numericInput(paste0("time",i), paste0("Time ",i," (hr):"),       value = NA))
      )
    }))
  })
  
  observeEvent(input$simulate, {
    obs_vals  <- sapply(seq_len(rv$n_obs), function(i) input[[paste0("obs",i)]])
    time_vals <- sapply(seq_len(rv$n_obs), function(i) input[[paste0("time",i)]])
    od <- na.omit(data.frame(conc=obs_vals, time=time_vals))
    
    if (!nrow(od)) { showModal(modalDialog("Enter at least one valid observation.")); return() }
    if (input$crcl <= 0) { showModal(modalDialog("CrCL must be > 0.")); return() }
    if (input$crrt==1 && input$flow <= 0) { showModal(modalDialog("CRRT flow must be > 0 when CRRT=Yes.")); return() }
    if (input$fu <= 0 || input$fu > 1) { showModal(modalDialog("fu must be between 0 and 1.")); return() }
    
    nd <- input$n_doses
    if (nd <= 0) { showModal(modalDialog("Number of doses must be at least 1.")); return() }
    if (input$dose <= 0 || input$infusion_time < 0 || input$tau <= 0) {
      showModal(modalDialog("Dose must be >0, infusion time >=0 (0=bolus), and tau >0.")); return()
    }
    
    if (isTRUE(input$use_load)) {
      ld <- as.numeric(input$ld_dose)
      dose_vec     <- c(ld, rep(input$dose, max(nd-1,0)))
      inf_time_vec <- c(0.5, rep(input$infusion_time, max(nd-1,0)))
      start_times  <- if (nd >= 2) c(0, 0.5 + (0:(nd-2)) * input$tau) else 0
    } else {
      dose_vec     <- rep(input$dose, nd)
      inf_time_vec <- rep(input$infusion_time, nd)
      start_times  <- (0:(nd-1)) * input$tau
    }
    
    # Manual path: constant covariates, force HD off
    cov_tbl <- data.frame(TIME = 0,
                          CRCL = as.numeric(input$crcl),
                          CRRT = as.integer(input$crrt),
                          FLOW = as.numeric(input$flow),
                          HD   = 0L)
    cov_at_obs <- eval_cov_at_times(od$time, cov_tbl)
    
    res <- tryCatch(
      run_map_meropenem(
        obs_df = od,
        dosing_df = data.frame(dose = dose_vec, tinf = inf_time_vec, start = start_times),
        cov_tbl = cov_tbl, cov_at_obs = cov_at_obs,
        fu = input$fu, priors = get_priors_list(),
        mic = input$mic, tau = input$tau
      ),
      error = function(e){ showModal(modalDialog("Error: ", e$message)); NULL }
    )
    if (is.null(res)) return()
    
    output$paramTable <- renderTable({
      transform(res$map_params, Value = sprintf("%.3f", as.numeric(Value)))
    }, bordered=TRUE, striped=TRUE)
    
    output$ftmicTable <- renderTable({
      data.frame(Metric = res$ftmic$Metric,
                 Value  = ifelse(is.na(res$ftmic$Value), "NA", sprintf("%.1f%%", res$ftmic$Value)))
    }, bordered=TRUE, striped=TRUE)
    
    output$concPlot <- renderPlotly({
      df <- data.frame(time=res$times, Individual=res$indiv, Population=res$pop)
      od2 <- res$obs
      p <- ggplot(df,aes(time))+
        geom_line(aes(y=Individual,color="Individual"))+
        geom_line(aes(y=Population,color="Population"))+
        geom_hline(yintercept=input$mic, linetype="dashed", color="gray30")+
        geom_point(data=od2,aes(x=time,y=conc), color="black", size=3)+
        scale_color_manual("",values=c(Individual="indianred",Population="dodgerblue"))+
        theme_minimal() + labs(x="Time (h)",y="Conc (mg/L)")
      ggplotly(p)
    })
    
    # NEW: IWRES plot
    output$iwresPlot <- renderPlotly({
      r <- res$resid_df
      p <- ggplot(r, aes(x = time, y = IWRES)) +
        geom_hline(yintercept = 0, linewidth = 0.4) +
        geom_hline(yintercept = c(-2, 2), linetype = "dashed", linewidth = 0.4) +
        geom_point(size = 2) +
        theme_minimal() +
        labs(x = "Time (h)", y = "IWRES")
      ggplotly(p)
    })
    
    updateTabsetPanel(session,"mainTabs","Plot")
  })
  
  # Upload path --------------------------------------------------------------
  upload_state <- reactiveValues(
    raw = NULL,
    subject_ids = NULL,
    last_run_subject = NULL,
    last_run_metrics = NULL,
    cache = new.env(parent=emptyenv()),
    global_metrics = NULL,
    global_predictions = NULL,
    run_log = NULL
  )
  
  observeEvent(input$parse_btn, {
    req(input$datafile)
    df <- tryCatch(readr::read_csv(input$datafile$datapath, show_col_types = FALSE),
                   error = function(e) NULL)
    if (is.null(df) || !nrow(df)) { showModal(modalDialog("Could not read CSV or file is empty.")); return() }
    
    names(df) <- toupper(names(df))
    
    if (!"AMT"  %in% names(df) && "DOSE" %in% names(df)) df$AMT  <- df$DOSE
    if (!"TINF" %in% names(df) && "DUR"  %in% names(df)) df$TINF <- df$DUR
    if (!"DV"   %in% names(df) && "OUT"  %in% names(df)) df$DV   <- df$OUT
    if (!"FLOW" %in% names(df) && "QB"   %in% names(df)) df$FLOW <- df$QB
    
    if (!"CRRT" %in% names(df)) {
      if (all(c("CVVH", "CVVHD", "CVVHDF") %in% names(df))) {
        df$CRRT <- dplyr::coalesce(df$CVVH, df$CVVHD, df$CVVHDF)
      } else if ("CVVH" %in% names(df)) {
        df$CRRT <- df$CVVH
      } else if ("CVVHD" %in% names(df)) {
        df$CRRT <- df$CVVHD
      } else if ("CVVHDF" %in% names(df)) {
        df$CRRT <- df$CVVHDF
      }
    }
    
    df <- df %>%
      mutate(
        ID   = as.character(ID),
        EVID = suppressWarnings(as.integer(EVID)),
        AMT  = suppressWarnings(as.numeric(AMT)),
        TINF = suppressWarnings(as.numeric(TINF)),
        TIME = suppressWarnings(as.numeric(TIME)),
        DV   = suppressWarnings(as.numeric(DV)),
        CRCL = suppressWarnings(as.numeric(CRCL)),
        CRRT = suppressWarnings(as.integer(CRRT)),
        FLOW = suppressWarnings(as.numeric(FLOW)),
        HD   = if ("HD" %in% names(.)) suppressWarnings(as.integer(HD)) else NA_integer_
      )
    
    upload_state$raw <- df
    upload_state$subject_ids <- sort(unique(df$ID))
    
    output$upload_preview <- renderDT({
      datatable(head(df, 50), options = list(scrollX = TRUE, pageLength = 5))
    })
    
    output$subject_picker <- renderUI({
      selectizeInput("subject_id", "Subject ID", choices = upload_state$subject_ids, multiple = FALSE)
    })
    output$report_subject_picker <- renderUI({
      selectizeInput("subject_id_report", "Subject ID", choices = upload_state$subject_ids, multiple = FALSE)
    })
  })
  
  cohort_prepare <- reactive({
    req(upload_state$raw)
    df <- upload_state$raw
    
    obs_counts <- df %>% filter(EVID==0, !is.na(DV)) %>% count(ID, name="n_levels_total")
    
    min_req <- as.integer(input$min_levels %||% 1)
    eligible <- obs_counts %>% filter(n_levels_total >= min_req)
    
    cap <- switch(input$level_cap,
                  cap1 = 1L, cap2 = 2L, cap3 = 3L,
                  all  = Inf,
                  1L)
    cap_label <- if (is.infinite(cap)) "all" else as.character(cap)
    
    capped_obs <- df %>%
      filter(EVID==0, !is.na(DV)) %>%
      semi_join(eligible, by="ID") %>%
      arrange(ID, TIME) %>%
      group_by(ID) %>%
      { if (is.infinite(cap)) . else slice_head(., n = cap) } %>%
      ungroup()
    
    list(
      eligible_subjects = eligible$ID,
      cap = cap,
      cap_label = cap_label,
      summary = eligible %>% arrange(ID),
      capped_obs = capped_obs
    )
  })
  
  observeEvent(input$preview_cohort, {
    co <- cohort_prepare()
    out <- co$summary %>%
      transmute(`Subject ID` = ID,
                `Total available levels` = n_levels_total,
                `Levels used (cap)` = co$cap_label)
    output$cohort_summary <- renderDT({
      datatable(out, options = list(pageLength = 10), rownames = FALSE)
    })
    updateTabsetPanel(session, "uploadTabs", selected = "Cohort Preview")
  })
  
  observeEvent(input$back_to_subject, {
    updateTabsetPanel(session, "uploadTabs", selected = "Subject Viewer")
  })
  
  observeEvent(input$run_subject_btn, {
    req(upload_state$raw, input$subject_id)
    
    sdf <- upload_state$raw %>% filter(ID == input$subject_id)
    
    dosing <- sdf %>%
      filter(EVID == 1, AMT > 0) %>%
      mutate(TINF = tidyr::replace_na(TINF, 0)) %>%
      transmute(dose = AMT, tinf = TINF, start = TIME) %>%
      arrange(start)
    
    cap <- switch(input$level_cap, cap1=1L, cap2=2L, cap3=3L, all=Inf, 1L)
    
    obs <- sdf %>%
      filter(EVID == 0, !is.na(DV)) %>%
      transmute(time = TIME, conc = DV) %>%
      arrange(time) %>%
      { if (is.infinite(cap)) . else slice_head(., n = cap) }
    
    if (!nrow(dosing)) { showModal(modalDialog("No dosing rows for selected subject.")); return() }
    if (!nrow(obs))    { showModal(modalDialog("No qualifying observation rows for selected subject.")); return() }
    
    cov_tbl    <- build_cov_changes(sdf)
    cov_at_obs <- eval_cov_at_times(obs$time, cov_tbl)
    fu   <- if (isTRUE(input$upload_fu_override)) input$fu else 0.98
    
    res <- tryCatch(
      run_map_meropenem(
        obs_df = obs, dosing_df = dosing,
        cov_tbl = cov_tbl, cov_at_obs = cov_at_obs,
        fu = fu, priors = get_priors_list(),
        mic = input$mic, tau = input$tau
      ),
      error = function(e){ showModal(modalDialog("Error: ", e$message)); NULL }
    )
    if (is.null(res)) return()
    
    upload_state$last_run_subject <- input$subject_id
    upload_state$last_run_metrics <- res$map_params
    upload_state$cache[[as.character(input$subject_id)]] <- res
    
    output$upload_paramTable <- renderTable({
      transform(res$map_params, Value = sprintf("%.3f", as.numeric(Value)))
    }, bordered=TRUE, striped=TRUE)
    
    output$upload_ftmicTable <- renderTable({
      data.frame(Metric = res$ftmic$Metric,
                 Value  = ifelse(is.na(res$ftmic$Value), "NA", sprintf("%.1f%%", res$ftmic$Value)))
    }, bordered=TRUE, striped=TRUE)
    
    output$upload_concPlot <- renderPlotly({
      dfp <- data.frame(time=res$times, Individual=res$indiv, Population=res$pop)
      p <- ggplot(dfp,aes(time))+
        geom_line(aes(y=Individual,color="Individual"))+
        geom_line(aes(y=Population,color="Population"))+
        geom_hline(yintercept=input$mic, linetype="dashed", color="gray30")+
        geom_point(data=res$obs,aes(x=time,y=conc), color="black", size=3)+
        scale_color_manual("",values=c(Individual="indianred",Population="dodgerblue"))+
        theme_minimal() + labs(x="Time (h)",y="Conc (mg/L)")
      ggplotly(p)
    })
    
    # NEW: IWRES plot (upload single subject)
    output$upload_iwresPlot <- renderPlotly({
      r <- res$resid_df
      p <- ggplot(r, aes(x = time, y = IWRES)) +
        geom_hline(yintercept = 0, linewidth = 0.4) +
        geom_hline(yintercept = c(-2, 2), linetype = "dashed", linewidth = 0.4) +
        geom_point(size = 2) +
        theme_minimal() +
        labs(x = "Time (h)", y = "IWRES")
      ggplotly(p)
    })
  })
  
  observeEvent(input$run_all_btn, {
    req(upload_state$raw)
    co <- cohort_prepare()
    ids <- co$eligible_subjects
    if (!length(ids)) { showModal(modalDialog("No eligible subjects under current filters.")); return() }
    
    pri <- get_priors_list()
    fu  <- if (isTRUE(input$upload_fu_override)) input$fu else 0.98
    
    metrics_list <- list()
    preds_list   <- list()
    log_list     <- list()
    
    withProgress(message = "Running subjects", value = 0, {
      n <- length(ids)
      for (i in seq_along(ids)) {
        incProgress(1/n, detail = paste("Subject", ids[i]))
        sid <- ids[i]
        sdf <- upload_state$raw %>% filter(ID == sid)
        
        dosing <- sdf %>%
          filter(EVID == 1, AMT > 0) %>%
          mutate(TINF = tidyr::replace_na(TINF, 0)) %>%
          transmute(dose = AMT, tinf = TINF, start = TIME) %>%
          arrange(start)
        
        obs <- co$capped_obs %>% filter(ID == sid) %>%
          transmute(time = TIME, conc = DV) %>% arrange(time)
        
        if (!nrow(dosing) || !nrow(obs)) {
          log_list[[length(log_list)+1]] <- data.frame(ID=sid, status="skipped", reason="missing dose/obs")
          next
        }
        
        cov_tbl    <- build_cov_changes(sdf)
        cov_at_obs <- eval_cov_at_times(obs$time, cov_tbl)
        
        res <- tryCatch(
          run_map_meropenem(
            obs_df = obs, dosing_df = dosing,
            cov_tbl = cov_tbl, cov_at_obs = cov_at_obs,
            fu = fu, priors = pri,
            mic = input$mic, tau = input$tau
          ),
          error = function(e) e
        )
        
        if (inherits(res, "error")) {
          log_list[[length(log_list)+1]] <- data.frame(ID=sid, status="error", reason=res$message)
          next
        }
        
        met <- res$map_params %>%
          pivot_wider(names_from = Parameter, values_from = Value) %>%
          mutate(
            `fT>MIC 0–24h (%)` = res$ftmic$Value[res$ftmic$Metric=="fT>MIC 0–24h"],
            `fT>MIC 24–48h (%)`= res$ftmic$Value[res$ftmic$Metric=="fT>MIC 24–48h"],
            `fT>MIC tau (%)`   = res$ftmic$Value[res$ftmic$Metric=="fT>MIC tau"]
          ) %>%
          mutate(ID = sid, .before = 1)
        
        metrics_list[[length(metrics_list)+1]] <- met
        
        # UPDATED: store predictions/residuals at obs times (includes IPRED/SD/IWRES)
        preds <- res$resid_df %>%
          mutate(ID = sid, .before = 1) %>%
          rename(time = time, obs_conc = conc) %>%
          transmute(ID, time, obs_conc, IPRED, PRED, SD, IWRES, PRED_SD, PWRES)
        
        preds_list[[length(preds_list)+1]] <- preds
        
        upload_state$cache[[as.character(sid)]] <- res
        log_list[[length(log_list)+1]] <- data.frame(ID=sid, status="ok", reason="")
      }
    })
    
    metrics_df <- if (length(metrics_list)) bind_rows(metrics_list) else tibble()
    preds_df   <- if (length(preds_list))   bind_rows(preds_list)   else tibble()
    log_df     <- if (length(log_list))     bind_rows(log_list)     else tibble(ID=character(), status=character(), reason=character())
    
    upload_state$global_metrics     <- metrics_df
    upload_state$global_predictions <- preds_df
    upload_state$run_log            <- log_df
    
    output$global_metrics_dt <- renderDT({
      req(nrow(metrics_df) > 0)
      datatable(
        metrics_df %>% mutate(across(where(is.numeric), ~round(.x, 3))),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    })
    
    output$run_log_dt <- renderDT({
      req(nrow(log_df) > 0)
      datatable(log_df, options = list(pageLength = 10))
    })
    
    updateTabsetPanel(session,"mainTabs","Reports")
  })
  
  observeEvent(input$toReports, updateTabsetPanel(session,"mainTabs","Reports"))
  
  output$recent_metrics <- renderTable({
    req(upload_state$last_run_metrics)
    transform(upload_state$last_run_metrics, Value = sprintf("%.3f", as.numeric(Value)))
  }, bordered=TRUE, striped=TRUE)
  
  output$download_subject_metrics <- downloadHandler(
    filename = function() {
      sid <- input$subject_id_report
      paste0("meropenem_subject_", sid %||% "NA", "_metrics.csv")
    },
    content = function(file) {
      sid <- input$subject_id_report
      validate(need(!is.null(sid), "Pick a subject."))
      res <- upload_state$cache[[as.character(sid)]]
      validate(need(!is.null(res), "Run this subject first."))
      met <- res$map_params %>%
        pivot_wider(names_from = Parameter, values_from = Value) %>%
        mutate(
          `fT>MIC 0–24h (%)` = res$ftmic$Value[res$ftmic$Metric=="fT>MIC 0–24h"],
          `fT>MIC 24–48h (%)`= res$ftmic$Value[res$ftmic$Metric=="fT>MIC 24–48h"],
          `fT>MIC tau (%)`   = res$ftmic$Value[res$ftmic$Metric=="fT>MIC tau"]
        ) %>% mutate(ID = sid, .before = 1)
      readr::write_csv(met, file)
    }
  )
  
  # UPDATED: subject predictions download includes IPRED/SD/IWRES
  output$download_subject_predictions <- downloadHandler(
    filename = function() {
      sid <- input$subject_id_report
      paste0("meropenem_subject_", sid %||% "NA", "_predictions_at_obs.csv")
    },
    content = function(file) {
      sid <- input$subject_id_report
      validate(need(!is.null(sid), "Pick a subject."))
      res <- upload_state$cache[[as.character(sid)]]
      validate(need(!is.null(res), "Run this subject first."))
      df <- res$resid_df %>%
        mutate(ID = sid, .before = 1) %>%
        rename(time = time, obs_conc = conc) %>%
        transmute(ID, time, obs_conc, IPRED, PRED, SD, IWRES, PRED_SD, PWRES)
      readr::write_csv(df, file)
    }
  )
  
  output$download_all_metrics_csv <- downloadHandler(
    filename = function() "meropenem_all_subjects_metrics.csv",
    content = function(file) {
      validate(need(!is.null(upload_state$global_metrics) && nrow(upload_state$global_metrics)>0,
                    "Run all subjects first."))
      readr::write_csv(upload_state$global_metrics %>%
                         mutate(across(where(is.numeric), ~round(.x, 6))), file)
    }
  )
  
  # UPDATED: global predictions download includes IWRES/IPRED/SD
  output$download_all_predictions_csv <- downloadHandler(
    filename = function() "meropenem_all_subjects_predictions_at_obs.csv",
    content = function(file) {
      validate(need(!is.null(upload_state$global_predictions) && nrow(upload_state$global_predictions)>0,
                    "Run all subjects first."))
      readr::write_csv(upload_state$global_predictions %>%
                         mutate(across(where(is.numeric), ~round(.x, 6))), file)
    }
  )
  
  observeEvent(input$toPriors,  updateTabsetPanel(session,"mainTabs","Population Priors"))
  observeEvent(input$toPlot,    updateTabsetPanel(session,"mainTabs","Plot"))
  observeEvent(input$toParams,  updateTabsetPanel(session,"mainTabs","Parameters"))
  observeEvent(input$toSummary, updateTabsetPanel(session,"mainTabs","Summary"))
}

shinyApp(ui, server)