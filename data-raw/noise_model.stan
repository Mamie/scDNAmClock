// saved as 8schools.stan
data {
  int<lower=1> I;          // number of probes
  int<lower=1> J;          // number of patients
  int<lower=1> K;          // number of replicates
  real<lower=0,upper=1> betahat[I,J,K];     // observed beta values
}

parameters {
  vector<lower=0>[I] sigma; // standard deviation of probe noise
  real<lower=0> eta;      // hyperparameter of sigma
  vector<lower=0>[I] tau;   // standard deviation of beta for each probe
  vector<lower=0, upper=1>[I] mu; // the mean of beta for each probe
  matrix<lower=0, upper=1>[I,J] beta; // the true beta for each probe, patient
}

model {
  eta ~ normal(0, 0.05);
  for (i in 1:I) {
    sigma[i] ~ normal(0, eta);
    tau[i] ~ normal(0, 0.05);
    mu[i] ~ normal(0, 1);
    for (j in 1:J) {
      beta[i,j] ~ normal(mu[i], tau[i]);
      for (k in 1:K) {
        betahat[i,j,k] ~ normal(beta[i,j], sigma[i]);
      }
    }
  }
}
