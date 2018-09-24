data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Nclusters1; // number of levels for group 1 for random intercepts
  int<lower=1> Nsizes; // number of size terms for controlling dispersion
  int<lower=1, upper=Nclusters1> NgroupMap1[Ntotal]; // mapping variable to map each observation to group 1 
  int<lower=1, upper=Nsizes> NsizeMap[Ntotal]; // mapping variable to map size to groups
  int<lower=1> Ncol; // total number of columns in model matrix
  int y[Ntotal]; // response variable
  // additional parameters
  real gammaShape; // hyperparameters for the gamma distribution 
  real gammaRate;
  real intercept;
  real intercept_sd;
}
// transformed data {
  // }
parameters {
  // parameters to estimate in the model
  vector[Ncol] betas; // regression intercept 
  real<lower=0.1> sigmaRan1; // random effect standard deviation for group 1
  real<lower=1, upper=1000> iSize[Nsizes]; // size parameter for the nb distribution
  vector[Nclusters1] rGroupsJitter1; // number of random jitters for each level of cluster/group 1
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  mu = exp(betas[1] + rGroupsJitter1[NgroupMap1]);
}
model {
  sigmaRan1 ~ gamma(gammaShape, gammaRate);
  betas ~ normal(intercept, intercept_sd);
  // random effects sample
  rGroupsJitter1 ~ normal(0, sigmaRan1);
  // likelihood function
  y ~ neg_binomial_2(mu, iSize[NsizeMap]);
}