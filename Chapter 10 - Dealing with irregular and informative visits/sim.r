require(nlme)
require(MASS)
require(truncnorm)


logit <- function(x) { 
  log(x) - log(1 - x)
}

expit <- function(x) { 
  1/(1 + exp(-x))
}

# Function to randomize treatment allocation
treatment_alloc_randomized <- function(age) {
  0.5
}

# Function to allocate treatment according to patient age
treatment_alloc_confounding <- function(age) {
  1/(1 + exp(-(0.7 - 0.032*age - 0.0001*(age**2))))
}

# convert correlation matrix to covariance matrix 
cor2cov <- function(sd, rho) {
  if (length(sd) != nrow(rho) & length(sd) != ncol(rho)) {
    stop("Invalid dimensions of 'sigma' and/or 'rho'")
  }
  D <- diag(sd, nrow = length(sd), ncol = length(sd))
  vmat <- D %*% rho %*% D
  vmat
}

convert_to_EDSS_scale <- function(x) {
  x <- round(x*2)/2
  x[which(x < 0)] <- 0
  x[which(x > 9.5)] <- 9.5
  x
}

sim_data_EDSS <- function(npatients = 500,
                          ncenters = 20,
                          follow_up = 12*5, # Total follow-up (number of months)
                          sd_a_t = 0.5,   # DGM - Within-visit variation in EDSS scores
                          baseline_EDSS = 1.3295,    # DGM - Mean baseline EDDS score
                          sd_alpha_ij = 1.46,    # DGM - Between-subject variation in baseline EDSS
                          sd_beta1_j = 0.20,    # DGM - Between-site variation in baseline EDSS
                          mean_age = 42.41,
                          sd_age = 10.53,
                          min_age = 18,
                          beta_age = 0.05, # DGM - prognostic effect of age
                          beta_t = 0.014,  # DGM - prognostic effect of time
                          beta_t2 = 0,    # DGM - prognostic effect of time squared
                          delta_xt = -0.007, # DGM - interaction treatment time
                          delta_xt2 = 0, # 0.0005    # DGM - interaction treatment time2
                          rho = 0.8,             # DGM - autocorrelation of between alpha_tij
                          corFUN = corAR1,       # DGM - correlation structure of the latent EDSS scores
                          tx_alloc_FUN = treatment_alloc_randomized # Treatment allocation function
                          ) 
{
  
  # Identify total number of patients
  n_total <- ncenters * npatients
  
  # Create the grid
  ytimes <- seq(from = 0, to = follow_up, by = 1)
  
  # Construct Sigma
  cs1 <- corFUN(value = rho)
  cor_matrix <- corMatrix(cs1, covariate = ytimes)
  sd_alpha <- rep(sd_a_t, (length(ytimes)))
  sigma_alpha <- cor2cov(sd = sd_alpha, rho = cor_matrix)
  
  # Draw a prognostic factor
  age <- rtruncnorm(n = n_total, 
                    a = min_age, 
                    b = Inf, 
                    mean = mean_age, 
                    sd = sd_age)
  
  # Allocate treatment for all patients. If applicable, make use of baseline info
  ptreat <- tx_alloc_FUN(age = age) 
  xtreat <- rbinom(n = n_total, size = 1, prob = ptreat)
  
  # Identify the centers
  centerid <- rep(seq(ncenters), each = npatients)
  
  # Draw the patient effects
  alpha_ij <- rnorm(n = n_total, mean = 0, sd = sd_alpha_ij)
  
  # Draw the center effects
  beta_1j <- rep(rnorm(ncenters, mean = 0, sd = sd_beta1_j), each = npatients)
  
  # Draw epsilon
  epsilon_tij_x0 <-  mvrnorm(n = n_total, mu = rep(0, length(ytimes)), Sigma = sigma_alpha)
  epsilon_tij_x1 <-  mvrnorm(n = n_total, mu = rep(0, length(ytimes)), Sigma = sigma_alpha)
  
  # Generate matrices with disease trajectory for eacht patient
  delta_baseline <- matrix(baseline_EDSS, nrow = n_total, ncol = length(ytimes))
  
  # Patient-specific baseline risk (constant  over time)
  delta_patient <- matrix(alpha_ij, nrow = n_total, ncol = length(ytimes), byrow = FALSE) 
  
  # Cluster-specific basline risk (constant over time)
  delta_cluster <- matrix(beta_1j, nrow = n_total, ncol = length(ytimes), byrow = FALSE) 
  
  # Time effect is identical for all patients, but varies over time
  delta_time <- matrix(beta_t * ytimes + beta_t2 * (ytimes**2), 
                       nrow = n_total, ncol = length(ytimes), byrow = TRUE) 
  
  # Treatment effect for received treatment
  
  delta_x1 <- matrix(delta_xt * ytimes + delta_xt2 * (ytimes**2), 
                     nrow = n_total, ncol = length(ytimes), byrow = TRUE)
  
  # Age effect
  delta_age <- matrix(beta_age * age, nrow = n_total, ncol = length(ytimes), byrow = FALSE) 
  
  latent_y_x0 <- delta_baseline + delta_patient + delta_cluster + delta_time + delta_age + epsilon_tij_x0
  latent_y_x1 <- delta_baseline + delta_patient + delta_cluster + delta_time + delta_x1 + delta_age + epsilon_tij_x1
  dsx <- cbind(l_x0 = as.vector(t(latent_y_x0)), 
               l_x1 = as.vector(t(latent_y_x1)), 
               x0 = rep(xtreat, each = length(ytimes)) == 0, 
               x1 = rep(xtreat, each = length(ytimes)) == 1,
               l_dr = NA)
  dsx[,"l_dr"] <- ifelse(dsx[,"x0"] == 1, dsx[,"l_x0"], dsx[,"l_x1"])
  
  
  mat <- matrix(NA, nrow = (n_total * length(ytimes)), ncol = 6)
  colnames(mat) <- c("centerid", 
                     "patid", # Patient ID
                     "x", # Received treatment
                     "age", # Age at treatment allocation
                     "time", # Visit time (#months since treatment allocation)
                     "y") #Observed EDSS outcome under received treatment
  
  mat[, "centerid"] <- rep(centerid, each = length(ytimes))
  mat[, "patid"] <- rep(seq_along(1:n_total), each = length(ytimes))
  mat[, "x"] <- rep(xtreat, each = length(ytimes))
  mat[, "age"] <- rep(age, each = length(ytimes))
  mat[, "time"] <- rep(ytimes, n_total)
  mat[, "y"] <- convert_to_EDSS_scale(dsx[,"l_dr"] ) 
  
  return(data.frame(mat))
}