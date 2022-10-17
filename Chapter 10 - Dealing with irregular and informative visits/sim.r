require(nlme)
require(MASS)
require(truncnorm)
require(dplyr)

logit <- function(x) { 
  log(x) - log(1 - x)
}

expit <- function(x) { 
  1/(1 + exp(-x))
}

# Function to randomize treatment allocation
treatment_alloc_randomized <- function(age, sex) {
  0.5
}

# Function to allocate treatment according to patient age
treatment_alloc_confounding <- function(age, sex) {
  1/(1 + exp(-(0.7 - 0.032*age - 0.0001*(age**2))))
}

# Treatment allocation v2 age, age squared and sex

treatment_alloc_confounding_v2 <- function( age, sex ) {
  1/(1 + exp(-(0.7 - 0.032*age - 0.0001*(age**2) + 0.2*sex)))
}

treatment_alloc_confounding_v2(30,1)

# convert correlation matrix to covariance matrix 
cor2cov <- function(sd, rho) {
  if (length(sd) != nrow(rho) & length(sd) != ncol(rho)) {
    stop("Invalid dimensions of 'sigma' and/or 'rho'")
  }
  D <- diag(sd, nrow = length(sd), ncol = length(sd))
  vmat <- D %*% rho %*% D
  return(vmat)
}

convert_to_EDSS_scale <- function(x) {
  x <- round(x*2)/2
  x[which(x < 0)] <- 0
  x[which(x > 9.5)] <- 9.5
  x
}

sim_data_EDSS <- function(npatients = 500, # Number of patients per center
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
                          delta_xt = 0, # DGM - interaction treatment time
                          delta_xt2 = 0, # 0.0005    # DGM - interaction treatment time2
                          p_female = 0.75, 
                          beta_female = -0.2 ,  ## DGM - prognostic effect of male sex
                          delta_xf = 0.1,      ## DGM - interaction sex treatment       
                          rho = 0.8,             # DGM - autocorrelation of between alpha_tij
                          corFUN = corAR1,       # DGM - correlation structure of the latent EDSS scores
                          tx_alloc_FUN = treatment_alloc_confounding_v2 # Treatment allocation function
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
  
  sex <- rbinom(n = n_total, 
                size=1,  prob= p_female)
  
  # Allocate treatment for all patients. If applicable, make use of baseline info
  ptreat <- tx_alloc_FUN(age, sex) 
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
  
  delta_x1 <- matrix( 1.4+ delta_xf * c(rep(sex, length(ytimes))) + delta_xt * ytimes + delta_xt2 * (ytimes**2), 
                      nrow = n_total, ncol = length(ytimes), byrow = TRUE)
  
  # Age effect
  delta_age <- matrix(beta_age * age, nrow = n_total, ncol = length(ytimes), byrow = FALSE) 
  
  # Sex effect
  delta_sex <- matrix(beta_female * sex, nrow= n_total, ncol = length(ytimes), byrow = FALSE) 
  
  latent_y_x0 <- delta_baseline + delta_patient + delta_cluster + delta_time + delta_age + delta_sex + epsilon_tij_x0
  latent_y_x1 <- delta_baseline + delta_patient + delta_cluster + delta_time + delta_x1 + delta_age + delta_sex + epsilon_tij_x1
  dsx <- cbind.data.frame(l_x0 = as.vector(t(latent_y_x0)), 
                          l_x1 = as.vector(t(latent_y_x1)), 
                          x0 = rep(xtreat, each = length(ytimes)) == 0, 
                          x1 = rep(xtreat, each = length(ytimes)) == 1,
                          l_dr = NA)
  dsx[,"l_dr"] <- ifelse(dsx[,"x0"] == 1, dsx[,"l_x0"], dsx[,"l_x1"])
  
  mat <- data.frame(centerid = rep(centerid, each = length(ytimes)), #center ID
                    patid = rep(seq(n_total), each = length(ytimes)), # Patient ID
                    x = rep(xtreat, each = length(ytimes)), # Received treatment
                    age = rep(age, each = length(ytimes)), # Age at treatment allocation
                    ages = age^2,
                    sex= rep(sex, each=length(ytimes)),
                    time = rep(ytimes, n_total), # Visit time (#months since treatment allocation)
                    edss = NA, # Baseline EDSS
                    y = convert_to_EDSS_scale(dsx[,"l_dr"]) ,  # Observed EDSS outcome under received treatment
                    progression = NA # Observed disease progression (0=no progression from baseline, 1=progression from baseline)
  )
  
  mat <- mat  %>% group_by(patid) %>% mutate(edss = first(y),
                                             progression = case_when((edss >= 6 & (y-edss) >= 0.5) | 
                                                                       (edss >= 1 & edss < 6 & (y-edss) >= 1.0) | 
                                                                       (edss < 1 & (y-edss) >= 1.5) ~ 1,
                                                                     TRUE ~ 0))
  
  
  return(mat)
}

# Patient visits are missing according to center
censor_visits_a1 <- function(data) {
  
  data$y_obs <- data$y
  
  ncenters <- length(unique(data$centerid))
  
  # Draw the center effects for informative censoring
  u <- rnorm(ncenters, mean = 0, sd = 0.15)
  
  # Calculate probability of missing
  data$prob_yobs <- expit(-1.94 + u[data$centerid])
  
  # By default, we always have a visit for time = 0
  data$prob_yobs[data$time == 0] <- 1
  
  # Set y_obs equal to NA where missing
  data$y_obs[rbinom(nrow(data), size = 1, prob = data$prob_yobs) == 0] <- NA
  
  data
}

# Patient visits are missing according to center and treatment
censor_visits_a2 <- function(data) {
  
  data$y_obs <- data$y
  
  ncenters <- length(unique(data$centerid))
  
  # Draw the center effects for informative censoring
  u <- rnorm(ncenters, mean = 0, sd = 0.15)
  
  # Calculate probability of missing
  data$prob_yobs <- expit(-1.6 + u[data$centerid] - data$x*0.7)
  
  # By default, we always have a visit for time = 0
  data$prob_yobs[data$time == 0] <- 1
  
  # Set y_obs equal to NA where missing
  data$y_obs[rbinom(nrow(data), size = 1, prob = data$prob_yobs) == 0] <- NA
  
  data
}

# Visit schedules are regular but differ between treatment groups
censor_visits_a3 <- function(data) {
  
  data$y_obs <- data$y
  data$prob_yobs <- 0.03 #changed
  data$prob_yobs[data$x == 0 & data$time %% 3 == 0] <- 0.35
  data$prob_yobs[data$x == 1 & data$time %% 9 == 0] <- 0.55
  data$prob_yobs[data$time == 0] <- 1
  
  # Set y_obs equal to NA where missing
  data$y_obs[rbinom(nrow(data), size = 1, prob = data$prob_yobs) == 0] <- NA
  
  data
}

# Patient visits are missing according to their received treatment and current EDSS score
censor_visits_a4 <- function(data) {
  
  data$y_obs <- data$y
  
  # Calculate probability of missing
  data$prob_yobs <- expit(-0.5  -  0.5 * data$y + 0.5 * data$x)
  
  # By default, we always have a visit for time = 0
  data$prob_yobs[data$time == 0] <- 1
  
  # Set y_obs equal to NA where missing
  data$y_obs[rbinom(nrow(data), size = 1, prob = data$prob_yobs) == 0] <- NA
  
  data
}



# Patient visits are missing according to their treatment, age and sex
censor_visits_a5 <- function(data) {
  
  data$y_obs <- data$y
  
  # Calculate probability of missing
  data$prob_yobs <- expit(-0.5 + 1.6 * data$x + 0.8 * data$sex -0.2 * data$age)
  
  # By default, we always have a visit for time = 0
  data$prob_yobs[data$time == 0] <- 1
  
  # Set y_obs equal to NA where missing
  data$y_obs[rbinom(nrow(data), size = 1, prob = data$prob_yobs) == 0] <- NA
  
  data
}







