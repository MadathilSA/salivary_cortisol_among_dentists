library(brms)

# Definition of variables in the model
# cortc - cortisol levels 
# b00age - Age (years)
# b01sex - Sex (1=male, 0=female)
# ddif  - follow-up days
# time_24 - time of saliva collection in 24 hour format
# vac_stat - vaccination status (0=not vaccinated, 1= atleast 1 dose against COVID19)


# BRMS lognormal mixed effect model with random intercept and effect of vaccination
model <- brm(cortc~ b00age + b01sex + s(ddif) + s(time_24) + (1+vac_stat|a00token),
             data=data, 
             family=lognormal(), 
             control = list(adapt_delta=0.85))


## Priors used for the model are also default priors in brms 
prior_summary(model)
# prior     class       coef    group resp dpar nlpar lb ub       source
# (flat)         b                                                                default
# (flat)         b     b00age                                             (vectorized)
# (flat)         b   b01sex02                                             (vectorized)
# (flat)         b    sddif_1                                             (vectorized)
# (flat)         b stime_24_1                                             (vectorized)
# student_t(3, -2.3, 2.5) Intercept                                                default
# lkj_corr_cholesky(1)         L                                                default
# lkj_corr_cholesky(1)         L            a00token                       (vectorized)
# student_t(3, 0, 2.5)        sd                                      0         default
# student_t(3, 0, 2.5)        sd            a00token                  0    (vectorized)
# student_t(3, 0, 2.5)        sd  Intercept a00token                  0    (vectorized)
# student_t(3, 0, 2.5)        sd  vac_stat1 a00token                  0    (vectorized)
# student_t(3, 0, 2.5)       sds                                      0         default
# student_t(3, 0, 2.5)       sds    s(ddif)                           0    (vectorized)
# student_t(3, 0, 2.5)       sds s(time_24)                           0    (vectorized)
# student_t(3, 0, 2.5)     sigma                                      0         default

## The model object produced internally in brms model

// generated with brms 2.16.3
functions {
  /* compute correlated group-level effects
  * Args: 
    *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns: 
    *   matrix of scaled group-level effects
  */ 
    matrix scale_r_cor(matrix z, vector SD, matrix L) {
      // r is stored in another dimension order than z
      return transpose(diag_pre_multiply(SD, L) * z);
    }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for splines
  int Ks;  // number of linear effects
  matrix[N, Ks] Xs;  // design matrix for the linear effects
  // data for spline s(ddif)
  int nb_1;  // number of bases
  int knots_1[nb_1];  // number of knots
  // basis function matrices
  matrix[N, knots_1[1]] Zs_1_1;
  // data for spline s(time_24)
  int nb_2;  // number of bases
  int knots_2[nb_2];  // number of knots
  // basis function matrices
  matrix[N, knots_2[1]] Zs_2_1;
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  vector[N] Z_1_2;
  int<lower=1> NC_1;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  vector[Ks] bs;  // spline coefficients
  // parameters for spline s(ddif)
  // standarized spline coefficients
  vector[knots_1[1]] zs_1_1;
  real<lower=0> sds_1_1;  // standard deviations of spline coefficients
  // parameters for spline s(time_24)
  // standarized spline coefficients
  vector[knots_2[1]] zs_2_1;
  real<lower=0> sds_2_1;  // standard deviations of spline coefficients
  real<lower=0> sigma;  // dispersion parameter
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
}
transformed parameters {
  // actual spline coefficients
  vector[knots_1[1]] s_1_1;
  // actual spline coefficients
  vector[knots_2[1]] s_2_1;
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_1;
  vector[N_1] r_1_2;
  // compute actual spline coefficients
  s_1_1 = sds_1_1 * zs_1_1;
  // compute actual spline coefficients
  s_2_1 = sds_2_1 * zs_2_1;
  // compute actual group-level effects
  r_1 = scale_r_cor(z_1, sd_1, L_1);
  r_1_1 = r_1[, 1];
  r_1_2 = r_1[, 2];
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + Xc * b + Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n];
    }
    target += lognormal_lpdf(Y | mu, sigma);
  }
  // priors including constants
  target += student_t_lpdf(Intercept | 3, -2.3, 2.5);
  target += student_t_lpdf(sds_1_1 | 3, 0, 2.5)
  - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_1_1);
  target += student_t_lpdf(sds_2_1 | 3, 0, 2.5)
  - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_2_1);
  target += student_t_lpdf(sigma | 3, 0, 2.5)
  - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sd_1 | 3, 0, 2.5)
  - 2 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(to_vector(z_1));
  target += lkj_corr_cholesky_lpdf(L_1 | 1);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
}
