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

ggplot_distribution_edss <- function(sim_data) {
  require(ggplot2)

  sim_data$Treatment <- factor(sim_data$x, levels = c(0,1), labels = c("DMT A", "DMT B"))
  sim_data$time <- as.factor(sim_data$time)

  simsum <- sim_data %>% group_by(time, Treatment) %>% summarize(
    median = median(y),
    mean = mean(y),
    lqd = quantile(y)["25%"],
    uqd = quantile(y)["75%"],
    miny = min(y),
    maxy = max(y),
    month = as.numeric(as.character(time))[1]
  )

  ggplot(simsum, aes(x = month, y = median, group = Treatment)) +
    geom_line() +
    geom_point(aes(y=miny)) +
    geom_point(aes(y=maxy)) +
    geom_line(aes(y=lqd), lty = 2, color = "grey") +
    geom_line(aes(y=uqd), lty = 2, color = "grey") +
    geom_ribbon(aes(ymin = lqd, ymax = uqd), fill = "grey", alpha = 0.5) +  # Shaded area
    xlab("Treatment exposure (months)") +
    ylab("Median EDSS") +
    facet_wrap(~ Treatment)  +
    theme_bw()
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
