// saved as 8schools.stan
data {
  int<lower=1> I;          // number of probes
  int<lower=1> J;          // number of patients
  int<lower=1> K;          // number of replicates
  real betahat[I,J,K];         // beta of first replicate
}

parameters {
  real<lower=0> sigma[I]; // standard deviation of probe noise
  real<lower=0> eta;      // standard deviation of sigma
  real<lower=0> tau[I];   // standard deviation of beta for each probe
  real<lower = 0, upper = 1> mu[I]; // the mean of beta for each probe
  real<lower = 0, upper = 1> beta[I, J]; // the true beta for each probe, patient
}

model {
  eta ~ normal(0, 0.1);
  for (i in 1:I) {
    sigma[i] ~ normal(0, eta);
    tau[i] ~ normal(0, 0.05);
    mu[i] ~ uniform(0, 1);
    for (j in 1:J) {
      beta[i, j] ~ normal(mu[i], tau[i]);
      for (k in 1:K) {
        betahat[i, j, k] ~ normal(betaij, sigma[i]);
      }
    }
  }
}