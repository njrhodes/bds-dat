// models/mem_plasma_model.stan  (new array syntax; bolus support; HD term)

data {
  int<lower=1> N;
  vector[N] time;
  vector[N] conc;

  int<lower=1> ndoses;
  vector[ndoses] dose_vec;
  vector[ndoses] infusion_vec;
  vector[ndoses] start_vec;

  // Per-observation, time-varying covariates
  vector<lower=0>[N] CRCL_i;
  array[N] int<lower=0, upper=1> CRRT_i;
  array[N] real<lower=0> FLOW_i;
  array[N] int<lower=0, upper=1> HD_i;

  // Priors (means & SD on log scale for log-normal IIV)
  real<lower=0> V_pop;    real<lower=0> V_sd;
  real<lower=0> CL1_pop;  real<lower=0> CL1_sd;
  real<lower=0> CL2_pop;  real<lower=0> CL2_sd;
  real<lower=0> pwr_pop;  real<lower=0> pwr_sd;
  real<lower=0> CLHD_pop; real<lower=0> CLHD_sd;
  
  // Residual error (FIXED for MAP TDM)
  real<lower=0> add_sigma;
  real<lower=0> prop_sigma;
}
parameters {
  real<lower=-3,upper=3> eta_V;
  real<lower=-3,upper=3> eta_CL1;
  real<lower=-3,upper=3> eta_CL2;
  real<lower=-3,upper=3> eta_pwr;
  real<lower=-3,upper=3> eta_CLHD;
}
transformed parameters {
  real<lower=1e-3> V    = V_pop    * exp(eta_V);
  real<lower=1e-3> CL1  = CL1_pop  * exp(eta_CL1);
  real<lower=1e-3> CL2  = CL2_pop  * exp(eta_CL2);
  real<lower=1e-3> pwr  = pwr_pop  * exp(eta_pwr);
  real<lower=1e-6> CLHD = CLHD_pop * exp(eta_CLHD);  // L/hr, active only when HD_i==1
}
model {
  // IIV priors
  eta_V    ~ normal(0, V_sd);
  eta_CL1  ~ normal(0, CL1_sd);
  eta_CL2  ~ normal(0, CL2_sd);
  eta_pwr  ~ normal(0, pwr_sd);
  eta_CLHD ~ normal(0, CLHD_sd);

  for (i in 1:N) {
    // Total clearance at observation i
    real CLT_i = fmax(
      CL1 * pow(CRCL_i[i] / 120, pwr) * (1 - CRRT_i[i]) +
      CL2 * CRRT_i[i] * (FLOW_i[i] / 1000) +
      CLHD * HD_i[i],
      1e-6
    );
    real Ke_i = fmax(CLT_i / V, 1e-6);

    real pred = 0;
    for (d in 1:ndoses) {
      real dt = time[i] - start_vec[d];
      if (dt > 0) {
        if (infusion_vec[d] > 0) {
          // zero-order infusion
          real rate = dose_vec[d] / infusion_vec[d] / CLT_i;
          if (dt <= infusion_vec[d]) {
            pred += rate * (1 - exp(-Ke_i * dt));
          } else {
            pred += rate * (1 - exp(-Ke_i * infusion_vec[d])) * exp(-Ke_i * (dt - infusion_vec[d]));
          }
        } else {
          // IV bolus (TINF == 0)
          pred += (dose_vec[d] / V) * exp(-Ke_i * dt);
        }
      }
    }
    pred = fmin(fmax(pred, 1e-8), 1e8);
    conc[i] ~ normal(pred, sqrt(square(add_sigma) + square(pred * prop_sigma)));
  }
}
