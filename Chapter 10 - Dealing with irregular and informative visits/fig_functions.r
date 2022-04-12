library(ggplot2)
library(dplyr)

plot_example_trajectory <- function(scenario) {
  logger <- flog.namespace() # Initiate Logger
  
  simpars <- load_scenario(scenario_no = scenario, logger = logger) 
  
  # Generate a dataset
  simpars$npatients <- 1
  simpars$ncenters <- 1
  sim_data <- sim_data_no_noise_V4(simpars, logger = logger)
  
  # Introduce missing values
  sim_data <- simpars$censorFUN(data = sim_data, outcome = "y_dr", logger = logger)
  
  ggplot(sim_data, aes(x=time, y=y_dr)) +
    geom_line(col="red") +
    geom_point(aes(x=time, y=y_cens), size =2) +
    xlab("Time") +
    ylab("EDSS")
  
}

plot_md_patterns <- function(scenarios, labels) {
  
  logger <- flog.namespace() # Initiate Logger
  
  ggdat <- NULL
  
  for (i in seq_along(scenarios)) {
    simpars <- load_scenario(scenario_no = scenarios[i], logger = logger) 
    
    # Generate a dataset
    simpars$npatients <- 2000
    simpars$ncenters <- 50
    sim_data <- sim_data_no_noise_V4(simpars, logger = logger)
    
    # Introduce missing values
    sim_data <- simpars$censorFUN(data = sim_data, outcome = "y_dr", logger = logger)
    
    dat_visits <- sim_data %>% group_by(time, x) %>% 
      summarise(pr_yobs = mean(prob_yobs),
                pr_yobs_iqrl = quantile(prob_yobs, probs = 0.25),
                pr_yobs_iqrm = quantile(prob_yobs, probs = 0.50),
                pr_yobs_iqru = quantile(prob_yobs, probs = 0.75),
                scenario = scenarios[i],
                label = labels[i])
    
    ggdat <- rbind(ggdat, dat_visits)
  }
  
  # Remove rows where time=0
  ggdat <- subset(ggdat, time > 0)
  
  # Generate factors for relevant variables
  ggdat$Treatment <- factor(ggdat$x, levels = c(0,1), labels = c("DMT A", "DMT B"))
  ggdat$label <- factor(ggdat$label, levels = labels, labels = labels)
  
  
  
  g <- ggplot(ggdat, aes(x = time, y = pr_yobs_iqrm, fill = Treatment, color = Treatment)) + 
    geom_errorbar(aes(group = time, ymax = pr_yobs_iqru, ymin = pr_yobs_iqrl)) +
    geom_point(aes(group = time)) + 
    ylab("Probability of observing the EDSS score") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    xlab("Time (months)") + 
    scale_x_continuous(breaks = c(1, 9, 18, 27, 36, 45, 54, 60)) +
    scale_color_manual(values = c("#2171b5", "#6baed6", "hotpink1", "hotpink4")) +
    theme(legend.position = "bottom") +
    facet_wrap(~label + Treatment)
  
  
  return(g)
}


plot_expected_disease_trajectory <- function(scenario_no, dir) {
  simpars <- load_scenario(scenario_no = scenario_no, logger = logger) 
  ytimes <- seq(from = 0, to = simpars$follow_up)
  m_given_x0 <- simpars$intercept + 
    simpars$beta_t * ytimes + simpars$beta_t2 * (ytimes**2) + 
    simpars$beta_age * simpars$mean_age 
  
  m_given_x1 <- m_given_x0  + simpars$delta_xt * ytimes + simpars$delta_xt2 * (ytimes**2)
  
  y_given_x0 <- convert_to_EDSS_scale(m_given_x0)
  y_given_x1 <- convert_to_EDSS_scale(m_given_x1)
  
  #MNAR pattern
  pobs_x0 <- expit(-0.5  -  0.5 * y_given_x0 + 0.5 * 0)
  pobs_x0[1] <- 1
  sum(pobs_x0)
  
  pobs_x1 <- expit(-0.5  -  0.5 * y_given_x1 + 0.5 * 1)
  pobs_x1[1] <- 1
  sum(pobs_x1)
  
  # Empirical visit count
  nrep <- 500
  nvisits <- array(NA, dim = c(nrep,2))
  
  for (i in 1:nrep) {
    sim_data <- sim_data_no_noise_V4(simpars, logger = logger)
    sim_data <- simpars$censorFUN(data = sim_data, outcome = "y_dr", logger = logger)
    
    dat_visits <- sim_data %>% group_by(patid) %>% 
      summarise(x = unique(x),
                nvisits_theory = sum(prob_yobs),
                nvisits_empirical = sum(!is.na(y_cens)))
    
    nvisits[i,1] <- mean((subset(dat_visits, x == 0))$nvisits_empirical)
    nvisits[i,2] <- mean((subset(dat_visits, x == 1))$nvisits_empirical)
  }
  
  nvisits
  
}

plot_distribution_age <- function(sim_data, scenario_no, dir) {
  sim_data <- subset(sim_data, time == 0)
  sim_data$Treatment <- factor(sim_data$x, levels = c(0,1), labels = c("DMT A", "DMT B"))
  
  mean_age <- data.frame(Treatment = c("DMT A", "DMT B"), 
                         wt = c(mean(subset(sim_data, x == 0)$age), 
                                mean(subset(sim_data, x == 1)$age)))
  
  g <- ggplot(sim_data, aes(x=age, color = Treatment, fill=Treatment)) + 
    geom_histogram(aes(y=..density..), position="identity", alpha=0.5, fill = "white")+
    geom_density(alpha=0.6)+
    facet_grid(Treatment ~ .) +
    geom_vline(aes(xintercept = wt), mean_age, lty = 2) +
    xlab("Age at treatment start") + 
    scale_color_brewer(palette="Paired") + 
    scale_fill_brewer(palette="Paired") + 
    theme(legend.position = "none")
  
  if (!missing(dir)) {
    ggsave(g, filename = paste("./figures/", dir, "/sc", scenario_no, "_age_distribution.pdf", sep = ""), width = 10, height = 10)
  } else {
    return(g)
  }
}

plot_distribution_edss <- function(sim_data) {
  sim_data <- subset(sim_data, time%%3 ==0)
  sim_data$Treatment <- factor(sim_data$x, levels=c(0,1), labels = c("DMT A", "DMT B"))
  sim_data$time <- as.factor(sim_data$time)
  
  # Describe the complete data
  g <- ggplot(sim_data, aes(x=time, y=y, fill=Treatment)) + 
    geom_boxplot(outlier.shape = NA, notch=TRUE) +
    ggtitle("Disease course") +
    xlab("Treatment exposure (no. months)") +
    ylab("Expanded Disability Status Scale") + 
    scale_fill_brewer(palette="Blues")
  return (g)
}


plot_distribution_lp <- function(sim_data, scenario_no, simpars, dir) {
  sim_data$lp <- simpars$intercept + simpars$beta_t * sim_data$time + (simpars$beta_t2 * sim_data$time2) + 
    (simpars$delta_xt * sim_data$time + simpars$delta_xt2 * sim_data$time2) * sim_data$x + simpars$beta_age * sim_data$age
  
  sim_data <- subset(sim_data, time %% 3 == 0)
  sim_data$Treatment <- factor(sim_data$x, levels = c(0,1), labels = c("DMT A", "DMT B"))
  sim_data$time <- as.factor(sim_data$time)
  
  # Describe the complete data
  # Histogram of y_dr, per time point
  g <- ggplot(sim_data, aes(x = time, y = lp, fill = Treatment)) + 
    geom_boxplot(outlier.shape = NA, notch = TRUE) +
    ggtitle(paste("Scenario", scenario_no)) +
    xlab("Treatment exposure (no. months)") +
    ylab("Expanded Disability Status Scale") + 
    scale_fill_brewer(palette = "Blues")
  
  ggsave(g, filename = paste("./figures/", dir, "/sc", scenario_no, "_lp_distribution.pdf", sep = ""), width = 15, height = 10)
}


plot_prob_visit <- function(sim_data, scenario_no, dir) {
  
  sim_data$Treatment <- factor(sim_data$x, levels = c(0,1), labels = c("DMT A", "DMT B"))
  
  dat_visits <- sim_data %>% group_by(patid) %>% 
    summarise(x = unique(x),
              nvisits_theory = sum(prob_yobs),
              nvisits_empirical = sum(!is.na(y_cens)))
  
  dat_text <- data.frame(
    label = paste("Simulated visits per patient: ",
                  round(c(mean(subset(dat_visits, x == 0)$nvisits_theory), 
                          mean(subset(dat_visits, x == 1)$nvisits_theory)),1)),
    Treatment = unique(sim_data$Treatment)[order(unique(sim_data$x))]
  )
  
  g <- ggplot(sim_data, aes(x = time, y = prob_yobs)) + 
    geom_boxplot(aes(group = time, fill = Treatment), outlier.shape = NA) + 
    ylab("Probability of observing the EDSS score") +
    xlab("Time (months)") + 
    facet_wrap(~Treatment) + 
    geom_text(
      data    = dat_text,
      mapping = aes(x = Inf, y = Inf, label = label),
      hjust   = 1.05,
      vjust   = 1.5
    ) +
    scale_fill_brewer(palette = "Blues") +
    theme(legend.position = "bottom")
  
  if (!missing(dir)) {
    ggsave(g, filename = paste("./figures/", dir, "/sc", scenario_no, "_prob_yobs.pdf", sep = ""), width = 20, height = 15)
  } else {
    return(g)
  }
  
}


plot_est_coverage <- function(results, dir, logger, reference_method = "Reference (RCT)") {
  if (length(unique(results$scenario_no)) != 1) {
    stop("Invalid set of scenarios")
  }
  
  scenario_no <- results$scenario_no[1]
  
  nsim <- length(unique(results$sim_id))
  
  # Load the reference model (the "truth")
  #load(paste("results/sim_", scenario_no, "_fit_ref.RData", sep = ""))
  
  
  if (results$delta_xt[1] != 0) {
    results$ref_hr <- exp(mean(subset(results, method == reference_method)$est_logHR))
  } else {
    results$ref_hr <- 1
  }
  
  results$incil <- results$est_HR_CIl <= results$ref_hr
  results$inciu <- results$est_HR_CIu >= results$ref_hr
  results$inci <- results$incil & results$inciu
  
  
  results$method <- factor(results$method, labels = unique(results$method), levels = unique(results$method))
  
  resultsci <- results %>% group_by(method) %>% summarize(coverage = mean(inci))
  
  g <- ggplot(resultsci, aes(x = method, y = coverage)) +
    geom_hline(yintercept = 0.95, lty = 2) +
    geom_text(aes(x = method, label = scales::percent(coverage), y = coverage ), vjust = -.5) +
    geom_col()
  
  ggsave(g, filename = paste("./figures/", dir, "/sc", scenario_no, "_sim_est_coverage.pdf", sep = ""), width = 15, height = 10)
  
}

plot_est_autocorr <- function(results, simpars, dir) {
  if (length(unique(results$scenario_no)) != 1) {
    stop("Invalid set of scenarios")
  }
  
  scenario_no <- results$scenario_no[1]
  
  results <- subset(results, (method == "LME-CE (PMM)"))
  
  range <- results$est_param_range
  
  ac <- data.frame(rho = numeric(), timediff = numeric(), source = numeric())
  
  # True autocorrelation from the DGM
  ac_struct <- simpars$corFUN(value = simpars$rho)
  cor_matrix <- corMatrix(ac_struct, covariate = 1:11)
  
  for (i in 1:10) {
    rho <- cor_decay_exp(0, i, range)
    ac <- rbind(ac, cbind(rho = rho, timediff = i, source = 0))
    ac <- rbind(ac, c(rho = cor_matrix[(i + 1),1], timediff = i, source = 1))
  }
  
  ac$timediff <- as.factor(ac$timediff)
  ac$source <- factor(ac$source, levels = c(0,1), labels = c("Estimated", "Data Generation Model"))
  
  g <- ggplot(ac, aes(x = timediff, y = rho)) +
    geom_boxplot(aes(color = source)) +
    xlab("Time elapsed (#months)") +
    ylab("Estimated autocorrelation")
  
  ggsave(g, filename = paste("./figures/", dir, "/sc", scenario_no, "_sim_est_autocorr.pdf", sep = ""), width = 15, height = 10)
}


plot_est_hr <- function(results, dir, reference_method = "Reference (RCT)") {
  
  if (length(unique(results$scenario_no)) != 1) {
    stop("Invalid set of scenarios")
  }
  
  scenario_no <- results$scenario_no[1]
  
  if (results$delta_xt[1] != 0) {
    ref_hr <- exp(mean(subset(results, method == reference_method)$est_logHR))
  } else {
    ref_hr <- 1
  }
  
  results$method <- factor(results$method, labels = unique(results$method), levels = unique(results$method))
  
  
  g <- ggplot(results, aes(x = method, y = est_HR)) +
    geom_hline(yintercept = ref_hr, lty = 2) +
    geom_boxplot() +
    ylab("Estimated HR") +
    xlab("Imputation method") 
  
  ## Derive MSE for each method
  results$bias_hr <- results$est_HR - ref_hr
  mse_results <- results %>% group_by(method) %>% summarize(max_est = max(est_HR),
                                                            median_bias = median(bias_hr), 
                                                            sd_bias = sd(bias_hr),
                                                            rmse = sqrt(mean((bias_hr)**2)))
  
  g <- g + geom_text(data = mse_results, aes(x = method, y = max(max_est)+0.01, 
                                             label = paste("RMSE =", format(rmse, digits=2, nsmall=2))), col='black')
  
  if (!missing(dir)) {
    ggsave(g, filename = paste("./figures/", dir, "/sc", scenario_no, "_sim_est_hr.pdf", sep = ""), width = 15, height = 10)
  } else {
    return (g)
  }
  
}

plot_bias_param_lmece <- function(results, dir, logger) {
  results <- subset(results, method == "LME-CE (PMM)")
  
  scenario_no <- results$scenario_no[1]
  simpars <- load_scenario(scenario_no = scenario_no, logger = logger) 
  
  dat <- data.frame(est = c(results$est_param_intercept, results$est_param_time, results$est_param_age,
                            results$est_param_time_x, results$est_param_range),
                    param = c(rep("Intercept", nrow(results)),
                              rep("beta_time", nrow(results)),
                              rep("beta_age", nrow(results)),
                              rep("beta_x_time", nrow(results)),
                              rep("range", nrow(results))))
  
  a_mean <- data.frame(param = c("Intercept", "beta_time", "beta_age", "beta_x_time", "range"),
                       ref_val = c(simpars$intercept, simpars$beta_t, simpars$beta_age, simpars$delta_xt, NA))
  
  
  g <- ggplot(dat, aes(y = est)) +
    geom_hline(data = a_mean, aes(yintercept = ref_val), lty = 2) +
    geom_boxplot() +
    ylab("Estimate") +
    xlab("Model parameter") +
    facet_wrap(~ param, scales = "free")
  
  ggsave(g, filename = paste("./figures/", dir, "/sc", scenario_no, "_sim_param_lmece.pdf", sep = ""), width = 15, height = 10)
}

plot_runtime <- function(results, dir, logger) {
  
  if (length(unique(results$scenario_no)) != 1) {
    stop("Invalid set of scenarios")
  }
  
  scenario_no <- results$scenario_no[1]
  
  results$method <- factor(results$method, labels = unique(results$method), levels = unique(results$method))
  
  
  g <- ggplot(results, aes(x = method, y = runtime_analysis)) +
    geom_boxplot() +
    ylab("Runtime (seconds)") +
    xlab("Imputation method") 
  
  ggsave(g, filename = paste("./figures/", dir, "/sc", scenario_no, "_sim_runtime.pdf", sep = ""), width = 15, height = 10)
}

plot_bias_hr <- function(results, dir, reference_method = "Reference (RCT)") {
  
  if (length(unique(results$scenario_no)) != 1) {
    stop("Invalid set of scenarios")
  }
  
  scenario_no <- results$scenario_no[1]
  
  methods <- unique(results$method)
  
  if (results$delta_xt[1] != 0) {
    #results$ref_loghr <- mean(subset(results, method == reference_method)$est_logHR)
    ref_hr <- exp(mean(subset(results, method == reference_method)$est_logHR))
  } else {
    ref_hr <- 1
  }
  
  #results$bias_loghr <- results$est_logHR - results$ref_loghr
  results$bias_hr <- results$est_HR - ref_hr
  
  results$method <- factor(results$method, labels = methods, levels = methods)
  
  g <- ggplot(results, aes(x = method, y = bias_hr)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_boxplot() +
    ylab("Bias of estimated HR") +
    xlab("Imputation method") 
  
  ## Derive MSE for each method
  mse_results <- results %>% group_by(method) %>% summarize(median_bias = median(bias_hr), 
                                                            sd_bias = sd(bias_hr),
                                                            rmse = sqrt(mean((bias_hr)**2)))
  
  g <- g + geom_text(data = mse_results, aes(x = method, y = median_bias+0.3*sd_bias, 
                                             label = paste("RMSE =", format(rmse, digits=2, nsmall=2))), col='black', size=5)
  
  # 
  
  ggsave(g, filename = paste("./figures/", dir, "/sc", scenario_no, "_sim_bias_hr.pdf", sep = ""), width = 15, height = 10)
}

plot_multiple_rmse_edss <- function(f1, f2, f3, title = "") {
  load(f1)
  results <- finalMatrix
  rm(finalMatrix)
  load(f2)
  results <- rbind(results, finalMatrix)
  rm(finalMatrix)
  load(f3)
  results <- rbind(results, finalMatrix)
  
  lbl_method <- c("Reference (OBS)",
                  "Reference (RCT)",
                  "LOCF",
                  "Rounding",
                  "LME-CE (EDSS conversion)",
                  "LME-CE (PMM)")
  new_lbl_method <- c("Reference (OBS)",
                      "Reference (RCT)",
                      "LOCF",
                      "Rounding",
                      "MLMI-RND",
                      "MLMI-PMM")
  
  results$method <- factor(results$method, labels = lbl_method, levels = lbl_method)
  results$Treatment <- factor(results$delta_xt)
  results$Treatment <- factor(results$Treatment, 
                              labels = c("None", "Moderate", "Strong"), 
                              levels = c(0, -0.007, -0.014))
  
  g <- ggplot(results, aes(x = method, y = rmse, fill = Treatment)) +
    geom_boxplot() +
    ylab("RMSE of imputed EDSS values") +
    xlab("Imputation method")  + 
    guides(fill = guide_legend(title = "Treatment Effect")) +
    theme(legend.position = "bottom") +
    ggtitle(title)
  g
}

plot_multiple_est_hr <- function(f1, f2, f3, title = "", reference_method = "Reference (RCT)") {
  load(f1)
  results <- finalMatrix
  rm(finalMatrix)
  load(f2)
  results <- rbind(results, finalMatrix)
  rm(finalMatrix)
  load(f3)
  results <- rbind(results, finalMatrix)
  
  lbl_method <- c("Reference (OBS)",
                  "Reference (RCT)",
                  "LOCF",
                  "Rounding",
                  "LME-CE (EDSS conversion)",
                  "LME-CE (PMM)")
  results$method <- factor(results$method, labels = lbl_method, levels = lbl_method)
  results$Treatment <- factor(results$delta_xt)
  results$Treatment <- factor(results$Treatment, 
                              labels = c("None", "Moderate", "Strong"), 
                              levels = c(0, -0.007, -0.014))
  
  ref_hr_strong <- exp(mean(subset(results, method == reference_method & delta_xt == -0.014)$est_logHR))
  ref_hr_moderate <- exp(mean(subset(results, method == reference_method & delta_xt == -0.007)$est_logHR))
  ref_hr_none <- 1
  
  ref_dat <- data.frame(Treatment = c("None", "Moderate", "Strong"),
                        ref = c(ref_hr_none, ref_hr_moderate, ref_hr_strong))
  ref_dat$Treatment <- factor(ref_dat$Treatment, 
                              labels = c("None", "Moderate", "Strong"), 
                              levels = c("None", "Moderate", "Strong"))
  
  
  g <- ggplot(results, aes(x = method, y = est_HR, fill = Treatment)) +
    geom_hline(data = ref_dat, aes(yintercept = ref, col = Treatment), lty = 2) +
    geom_boxplot() +
    ylab("Estimated HR") +
    xlab("Imputation method") + 
    guides(fill=guide_legend(title="Treatment Effect")) +
    theme(legend.position = "bottom") +
    ggtitle(title)+ guides(col = FALSE)
  #facet_wrap( ~ Treatment, ncol = 1) 
  
  
  
  g
}

plot_rmse_edss <- function(results, dir) {
  
  if (length(unique(results$scenario_no)) != 1) {
    stop("Invalid set of scenarios")
  }
  
  scenario_no <- results$scenario_no[1]
  
  results$method <- factor(results$method, labels = unique(results$method), levels = unique(results$method))
  
  results <- subset(results, !(method %in% c("LOCF (1mo grid)", "Reference")))
  
  g <- ggplot(results, aes(x = method, y = rmse)) +
    geom_boxplot() +
    ylab("RMSE of imputed EDSS values at") +
    xlab("Imputation method") 
  
  if (!missing(dir)) {
    ggsave(g, filename = paste("./figures/", dir, "/sc", scenario_no, "_sim_rmse_edss.pdf", sep = ""), width = 15, height = 10)
  } else {
    return(g)
  }
  
  
}


plot_max_fup <- function(sim_data, scenario_no, dir) {
  
  # For each patient, only keep the last observed value of y_dr
  #sim_data$centerid <- factor(sim_data$centerid)
  sim_data <- subset(sim_data, !is.na(y_cens))
  sim_data$Treatment <- factor(sim_data$x, levels = c(0,1), labels = c("DMT A", "DMT B"))
  
  plotdat <- sim_data %>%
    group_by(patid, x, centerid) %>%
    summarize(max_fup = max(time), Treatment = unique(Treatment))
  
  # Estimate median of max_fup for each center
  mfup <- plotdat %>%
    group_by(centerid) %>%
    summarize(median_mfup = median(max_fup)) %>% 
    arrange(median_mfup)
  
  plotdat$centerid <- factor(plotdat$centerid, labels = mfup$centerid, levels = mfup$centerid)
  
  
  g <- ggplot(plotdat, aes(x = centerid, y = max_fup, fill = Treatment)) + 
    geom_boxplot() + 
    ylab("Maximum follow-up (months)") +
    xlab("MS center") +
    scale_fill_brewer(palette = "Blues")
  
  ggsave(g, filename = paste("./figures/", dir, "/sc", scenario_no, "_max_fup.pdf", sep = ""), width = 15, height = 10)
  
}

print_rmse_hr <- function(results,  reference_method = "Reference (RCT)") {
  if (length(unique(results$scenario_no)) != 1) {
    stop("Invalid set of scenarios")
  }
  
  methods <- unique(results$method)
  
  if (results$delta_xt[1] != 0) {
    #results$ref_loghr <- mean(subset(results, method == reference_method)$est_logHR)
    ref_hr <- exp(mean(subset(results, method == reference_method)$est_logHR))
  } else {
    ref_hr <- 1
  }
  
  #results$bias_loghr <- results$est_logHR - results$ref_loghr
  results$bias_hr <- results$est_HR - ref_hr
  
  results$method <- factor(results$method, labels = methods, levels = methods)
  
  
  ## Derive MSE for each method
  mse_results <- results %>% group_by(method) %>% summarize(median_bias = median(bias_hr), 
                                                            sd_bias = sd(bias_hr),
                                                            rmse = sqrt(mean((bias_hr)**2)))
  
  return(mse_results)
}

print_coverage_hr <- function(results, reference_method = "Reference (RCT)") {
  if (length(unique(results$scenario_no)) != 1) {
    stop("Invalid set of scenarios")
  }
  
  
  if (results$delta_xt[1] != 0) {
    results$ref_hr <- exp(mean(subset(results, method == reference_method)$est_logHR))
  } else {
    results$ref_hr <- 1
  }
  
  results$incil <- results$est_HR_CIl <= results$ref_hr
  results$inciu <- results$est_HR_CIu >= results$ref_hr
  results$inci <- results$incil & results$inciu
  
  
  results$method <- factor(results$method, labels = unique(results$method), levels = unique(results$method))
  
  resultsci <- results %>% group_by(method) %>% summarize(coverage = mean(inci))
  
  return(resultsci)
}