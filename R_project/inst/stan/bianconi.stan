functions {
  real g(real a, real b) {
    return a / (a + b);
  }
}
data {
  int<lower=1> N;
  vector[N] egfr;
  vector[N] igfr;
  vector[N] sos;
  vector[N] ras;
  vector[N] pi3k;
  vector[N] akt;
  vector[N] raf;
  vector[N] mek;
  vector[N] erk;
  vector[N] p90;
  vector[N] rasgap;
  vector[N] rafpp;
  vector[N] pp2a;
}
parameters {
  real <lower=0.0, upper=8.0>  b_egfr_sos;
  real <lower=0.0, upper=8.0>  b_igfr_sos;
  real <lower=0.0, upper=8.0>  b_egfr_pi3k;
  real <lower=0.0, upper=8.0>  b_igfr_pi3k;
  real <lower=0.0, upper=8.0>  b_sos_ras;
  real <lower=0.0, upper=8.0>  b_ras_pi3k;
  real <lower=0.0, upper=8.0>  b_pi3k_akt;
  real <lower=0.0, upper=8.0>  b_akt_raf;
  real <lower=0.0, upper=8.0>  b_ras_raf;
  real <lower=0.0, upper=8.0>  b_raf_mek;
  real <lower=0.0, upper=8.0>  b_mek_erk;
}
transformed parameters{
  vector[N] sos_mu;
  vector[N] ras_mu;
  vector[N] pi3k_mu;
  vector[N] akt_mu;
  vector[N] raf_mu;
  vector[N] mek_mu;
  vector[N] erk_mu;
  for(i in 1:N){
    sos_mu[i] = fmax(120000 * g(b_egfr_sos * egfr[i] + b_igfr_sos * igfr[i], p90[i]), .0001);
    ras_mu[i] = fmax(120000 * g(b_sos_ras * sos[i], rasgap[i]), .0001);
    pi3k_mu[i] = fmax(120000 * g(b_egfr_pi3k * egfr[i] + b_ras_pi3k * ras[i] + b_igfr_pi3k * igfr[i], 1.0), .0001);
    akt_mu[i] = fmax(600000 * g(b_pi3k_akt * pi3k[i], 1.0), .0001);
    raf_mu[i] = fmax(120000 * g(b_ras_raf * ras[i], b_akt_raf * akt[i] + rafpp[i]), .0001);
    mek_mu[i] = fmax(600000 * g(b_raf_mek * raf[i], pp2a[i]), .0001);
    erk_mu[i] = fmax(600000 * g(b_mek_erk * mek[i], pp2a[i]), .0001);
  }
}
model {
  for(i in 1:N){
    sos[i] ~ normal(sos_mu[i], sqrt(sos_mu[i]));
    ras[i] ~ normal(ras_mu[i], sqrt(ras_mu[i]));
    pi3k[i] ~ normal(pi3k_mu[i], sqrt(pi3k_mu[i]));
    akt[i] ~ normal(akt_mu[i], sqrt(akt_mu[i]));
    raf[i] ~ normal(raf_mu[i], sqrt(raf_mu[i]));
    mek[i] ~ normal(mek_mu[i], sqrt(mek_mu[i]));
    erk[i] ~ normal(erk_mu[i], sqrt(erk_mu[i]));
  }
}
generated quantities{
  vector[N] sos_sim;
  vector[N] ras_sim;
  vector[N] pi3k_sim;
  vector[N] akt_sim;
  vector[N] raf_sim;
  vector[N] mek_sim;
  vector[N] erk_sim;
  for(i in 1:N){
    sos_sim[i] = poisson_rng(sos_mu[i]);
    ras_sim[i] = poisson_rng(ras_mu[i]);
    pi3k_sim[i] = poisson_rng(pi3k_mu[i]);
    akt_sim[i] = poisson_rng(akt_mu[i]);
    raf_sim[i] = poisson_rng(raf_mu[i]);
    mek_sim[i] = poisson_rng(mek_mu[i]);
    erk_sim[i] = poisson_rng(erk_mu[i]);
  }
}
