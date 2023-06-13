library(ggplot2)
library(dplyr)

plot_distribution_progression <- function(sim_data) {
  sim_data <- subset(sim_data, time %% 3 == 0)
  
  
  ggdat <- sim_data %>% 
    group_by(x, time) %>% 
      dplyr::summarize(perc = sum(progression) / n()) 
  ggdat$Treatment <- factor(ggdat$x, levels = c(0,1), labels = c("DMT A", "DMT B"))
  ggdat$time <- as.factor(ggdat$time)
  
  # Describe the complete data
  g <- ggplot(ggdat, aes(x = time, y = perc, fill = Treatment)) + 
    geom_bar(stat = "identity", position = position_dodge()) +
    ggtitle("EDSS Progression from baseline") +
    xlab("Treatment exposure (no. months)") +
    ylab("Risk of EDSS progression") + 
    scale_fill_brewer(palette = "Blues")
  return(g)
}

plot_distribution_edss <- function(sim_data) {
  sim_data <- subset(sim_data, time %% 3 == 0)
  sim_data$Treatment <- factor(sim_data$x, levels = c(0,1), labels = c("DMT A", "DMT B"))
  sim_data$time <- as.factor(sim_data$time)
  
  # Describe the complete data
  g <- ggplot(sim_data, aes(x = time, y = y, fill = Treatment)) + 
    geom_boxplot(outlier.shape = NA, notch = TRUE) +
    ggtitle("EDSS Prognosis") +
    xlab("Treatment exposure (no. months)") +
    ylab("Expanded Disability Status Scale") + 
    scale_fill_brewer(palette = "Blues")
  return(g)
}



plot_prob_visit <- function(sim_data) {
  
  sim_data$Treatment <- factor(sim_data$x, levels = c(0,1), labels = c("DMT A", "DMT B"))
  
  dat_visits <- sim_data %>% group_by(patid) %>% 
    dplyr::summarize(x = unique(x),
              nvisits_theory = sum(prob_yobs),
              nvisits_empirical = sum(!is.na(y_obs)))
  
  dat_text <- data.frame(
    label = paste("Simulated visits per patient: ",
                  round(c(mean(subset(dat_visits, x == 0)$nvisits_theory), 
                          mean(subset(dat_visits, x == 1)$nvisits_theory)),1)),
    Treatment = unique(sim_data$Treatment)[order(unique(sim_data$x))]
  )
  
  g <- ggplot(sim_data, aes(x = time, y = prob_yobs)) + 
    geom_boxplot(aes(group = time, fill = Treatment), outlier.shape = NA) + 
    ylab("Probability of a patient visit") +
    xlab("Treatment exposure (no. months)") + 
    facet_wrap(~Treatment) + 
    geom_text(
      data    = dat_text,
      mapping = aes(x = Inf, y = Inf, label = label),
      hjust   = 1.05,
      vjust   = 1.5
    ) +
    scale_fill_brewer(palette = "Blues") +
    theme(legend.position = "bottom")
  
  return(g)
}
