#rm(list=ls())

library(dplyr)
library(tidyr)
library(data.table)
library(survey)
library(MASS)
library(truncnorm)
library(optmatch) #match
library(MatchIt)
library(WeightIt) #IPW
library(cobalt) # imbalance covariates
library(mice)
library(MatchThem)
library(PSweight)
library(sandwich) 
library(missForest)
library(ggplot2)

#F0. Function to simulate MS data
generate_data <- function(n, 
                         seed = NA,
                         beta = c(-0.2, -0.2, -0.2, -0.2, -0.2), # to create heterogeneity treatment effect
                         beta.x = c(-1.54, -0.01, 0.06, 0.25, 0.5, 0.13, 0.0000003), # to calculate the outcome 
                         percentiles = seq(0, 1, by = 0.2)){
  
  #' Generate simulated count data with settings based on real-world data
  #' Assume randomized treatment and independent covariates
  #'
  #' @param n sample size; integer
  #' @param seed randomization seed; integer
  #' @param beta coefficients characterizing treatment effect heterogeneity; vector of length 5
  #'             beta[1]*trt*I(high responder to DMF) +
  #'             beta[2]*trt*I(moderate responder to DMF) +
  #'             beta[3]*trt*I(neutral) +
  #'             beta[4]*trt*I(moderate responder to TERI) +
  #'             beta[5]*trt*I(high responder to TERI)
  #'             In the absence of treatment effect heterogeneity, set all beta[1:5] with the same values
  #' @param beta.x coefficients for main effects of other covariates in the rate; vector of length 7
  #'               beta.x[1] (intercept)
  #'               beta.x[2]*ageatindex_centered
  #'               beta.x[3]*female
  #'               beta.x[4]*prerelapse_num
  #'               beta.x[5]*prevDMTefficacy== hiiumhigh (reference low efficacy)
  #'               beta.x[6]*prevDMTefficacy== none (reference low efficacy)
  #'               beta.x[7]*premedicalcost
  #' @param percentiles percentiles (cumulative, monotone increasing) to define each responder subgroup; vector of floats of length 6
  #' @return data - a data frame of \code{n} rows and 8 columns of simulated data; data.frame
  
  
  set.seed(seed)
  
  if (percentiles[1] != 0 | percentiles[length(percentiles)] != 1 | length(percentiles) != 6) {
    stop("Wrong values of percentiles!")
  }
  
  # Create an empty shell
  ds <- data.table(NULL)
  
  # Define X, A, and time
  ds[, female := rbinom(n = n, size = 1, prob = 0.75)]
  ds[, ageatindex_centered := round(rtruncnorm(n, a = 18, b = 64, mean = 48, sd = 12), 0) - 48 ] # rounded to integers
  ds[, prerelapse_num := rpois(n = n, lambda = 0.44)]
  ds[, prevDMTefficacy := sample(x = c("None", "Low_efficacy", "Medium_high_efficacy"), #previous DM treatment efficacy
                               size = n, replace = TRUE, prob = c(0.45, 0.44, 0.11))]
  ds[, prevDMTefficacy := factor(prevDMTefficacy, labels = c("None", "Low_efficacy", "Medium_high_efficacy"))]
  ds[, premedicalcost := pmin(round(exp(rnorm(n = n, mean = 8.9, sd = 1.14)), 2), 600000)] # rounded to 2 decimal points
  ds[, numSymptoms := sample(x = c("0", "1", ">=2"),
                            size = n, replace = TRUE, prob = c(0.67, 0.24, 0.09))] # nuisance variable; do not include in the score
  ds[, numSymptoms := factor(numSymptoms, labels = c("0", "1", ">=2"))]
  ds[, finalpostdayscount := ceiling(rgamma(n = n, shape = 0.9, scale = 500))] # rounded up to integers
  ds[, finalpostdayscount := ifelse(finalpostdayscount > 2096, 2096, finalpostdayscount)] # truncate at the max follow up day, 2096
  ds[, finalpostdayscount := ifelse((finalpostdayscount > 2090) & (runif(1, 0, 1) < .5), 29, finalpostdayscount)] # mimic the 1 month peak;  move roughly half of the large values to 29

  # Define treatment allocation 

  XB <- model.matrix(~.,ds) %*% c(1.22,0.3,-0.1,0.2, 0.35,0.7,-0.00005,0.17,0.02,0)# ~75% people allocated in DMF arm based on (age,female,prerelapse_num,DMT efficacy,costs,numSymptoms)
  pi <- exp(XB)/(1 + exp(XB))
  ds[, trt := as.numeric(runif(n) <=pi)]
  ds[, treatment := as.factor(ifelse(trt == 1, "DMF", "TERI"))] 
  
  # Define Y (using all above PS predictors except for numSymptoms)
  xmat.score <- as.matrix(model.matrix(~ ageatindex_centered + female + prerelapse_num + prevDMTefficacy + premedicalcost, ds)) 
  gamma <- matrix(c(-0.33, # Intercept
                    -0.001, # Age
                    0.05, # female
                    -0.002, # prerelapse_num
                    0.33, 0.02, # Medium/high and none DMT efficacy
                    -0.0000005), nrow = 7) # premedicalcost
  score = exp(xmat.score %*% gamma)
  ds[, Iscore := cut(score, quantile(score, percentiles), include.lowest = T, labels = seq(1, 5))]
  
  xmat.rate <- as.matrix(model.matrix(~ ageatindex_centered + female + prerelapse_num + prevDMTefficacy + premedicalcost +
                              Iscore + trt + trt*Iscore, ds))
  
  betas <- matrix(c(beta.x[1], # Intercept
                    beta.x[2], # Age
                    beta.x[3], # female
                    beta.x[4], # prerelapse_num
                    beta.x[5], 
                    beta.x[6], # Medium/high and none DMT efficacy
                    beta.x[7], # premedicalcost
                    # Previous beta.x: -1.54 Intercept, -0.01 Age, 0.06 female, 0.47 prerelapse_num, 0.68, 0.13 Medium/high and none DMT efficacy, 0.000003 premedicalcost
                    0, 0, 0, 0, # Iscore categories 2:5 (reference group is Iscore1, high responders to DMF)
                    beta[1], beta[2] - beta[1], beta[3] - beta[1], beta[4] - beta[1], beta[5] - beta[1]))
  rate <-  exp(xmat.rate %*% betas)
  ds[, y := rnegbin(n = n, mu = rate*finalpostdayscount/365.25, theta = 3)] # post treatment number of relapses
 
  ds[, Iscore := factor(Iscore, labels = c("High A1","Moderate A1","Neutral","Moderate A0","High A0"))]
  ds[, years := finalpostdayscount / 365.25]
  ds[, age := ageatindex_centered + 48]
  data <- ds[,c("age","female", "prevDMTefficacy", "premedicalcost", "numSymptoms", "prerelapse_num", "treatment", "y", "years","Iscore")]
  return(data)
}


# F1.1 Function to transform database with dummy to categorical variables----

databack <- function(data) {
  data<-setDT(data)
  data[,pmar:=NULL]
  data[,pmart:=NULL]
  data[,pmart2:=NULL]
  data[,pmnar:=NULL]
  data[,female:=as.factor(female)]
  if (!is.factor(data$prevDMTefficacy)){
    data[, prevDMTefficacy := factor(prevDMTefficacy, labels = c("None", "Low_efficacy", "Medium_high_efficacy"))]}
  if (!is.factor(data$numSymptoms)){
    data[, numSymptoms := factor(numSymptoms, labels = c("0", "1", ">=2"))]  }  
  if (!is.factor(data$treatment)){
    data[, treatment := factor(treatment, labels = c("TERI", "DMF"))]} 
  
  data[,DMF := as.numeric(treatment == "DMF")]
  data[,pde.ind := as.factor(as.numeric(is.na(prevDMTefficacy) == FALSE))]  # indicator of previous DMT efficacy
  data[,pde.mi := as.factor(ifelse(is.na(prevDMTefficacy),"na",prevDMTefficacy))]  # missing indicator replacement of previous DMT efficacy
  data[,pmc.ind := as.factor(as.numeric(is.na(premedicalcost) == FALSE))] # indicator of premedical cost
  data[,pmc.mi := ifelse(is.na(premedicalcost),0,premedicalcost)] # missing indicator replacement of premedical cost
  data[,pns.ind := as.factor(as.numeric(is.na(numSymptoms) == FALSE))] # indicator of number of symptoms
  data[,pns.mi := as.factor(ifelse(is.na(numSymptoms),"na",numSymptoms))] # missing indicator replacement of number of symptoms
  data[,prn.ind := as.factor(as.numeric(is.na(prerelapse_num) == FALSE))] # indicator of prerelapse number
  data[,prn.mi := ifelse(is.na(prerelapse_num),0,prerelapse_num)]# missing indicator replacement of prerelapse number
  return(data)
}



#F2. Function to generate missing data ---


getmissdata <- function(data, scenario = "MAR", seed = 12345){
  set.seed(seed)
  
  # Transform categorical data to dummies
  dat_misspattern <- data %>% mutate(treatment = ifelse(treatment == "DMF", 1, 0)) %>%
    mutate(dummy=1) %>%
    spread(key=prevDMTefficacy, value=dummy, fill=0, sep ="_") %>%
    mutate(dummy=1) %>%
    spread(key=numSymptoms, value=dummy, fill=0, sep ="_") %>%
    mutate(dummy=1) %>%
    spread(key=Iscore, value=dummy, fill=0, sep ="_")
  
  # Generate missing values for premedical cost (MCAR)
  pattern <- rep(1, ncol(dat_misspattern))
  pattern[which( colnames(dat_misspattern)=="premedicalcost")] <- 0
  md1 <- ampute(dat_misspattern, patterns = pattern, prop = 0.10, mech = "MCAR")#$amp
  
  # Generate missing values for prevDMTefficacy and numSymptoms (MAR)
  cols_prevDMTefficacy <- grepl( "prevDMTefficacy", colnames(dat_misspattern), fixed = TRUE)
  cols_numSymptoms <- grepl( "numSymptoms", colnames(dat_misspattern), fixed = TRUE)
  pattern <- data.frame(matrix(1,nrow = 3,  ncol(dat_misspattern)))
  pattern[c(1,3),which(cols_prevDMTefficacy)] <- 0
  pattern[c(2,3),which(cols_numSymptoms)] <- 0
  # Alter the weights such that missingness only depends on observed values of 
  # * prevDMTefficacy (for pattern 2)
  # * numSymptoms (for pattern 1)
  # * age (for all patterns)
  # * female (for all patterns)
  md_temp <- ampute(dat_misspattern, patterns = pattern, prop = 0.3, mech = "MAR")
  weights <- data.frame(matrix(0,nrow = 3,  ncol(dat_misspattern)))
  colnames(weights) <- colnames(md_temp$weights)
  weights$age <- -0.2
  weights$female <- 0.3
  weights[2, cols_prevDMTefficacy] <- 0.2
  weights[1, cols_numSymptoms] <- 0.1
  md2 <- ampute(dat_misspattern, patterns = pattern, weights = weights, prop = 0.3, mech = "MAR")
  
  # Set missing data for prerelapsenum
  pattern <- rep(1, ncol(dat_misspattern))
  pattern[which( colnames(dat_misspattern)=="prerelapse_num")] <- 0
  
  weights <- rep(0, ncol(dat_misspattern))
  names(weights) <- colnames(md_temp$weights)
  
  if (scenario == "mcar") {
    mech <- "MCAR"
  } else if (scenario == "MAR") {
    weights["age"] <- 1/48
    weights["female"] <- 1
    mech <- "MAR"
  } else if (scenario == "MART") {
    weights["age"] <- 0.5/48
    weights["female"] <- 0.4
    weights["treatment"] <- 1
    mech <- "MAR"
  } else if (scenario == "MARTY") {
    weights["age"] <- 0.5/48
    weights["female"] <- 0.4
    weights["treatment"] <- 1
    weights["y"] <- 3
    mech <- "MAR"
  } else if (scenario == "MNAR") {
    weights["age"] <- 0.5/48
    weights["female"] <- 0.4
    weights["prerelapse_num"] <- 4
    mech <- "MNAR"
  } else {
    stop("Scenario not supported!")
  }
  md3 <- ampute(dat_misspattern, patterns = pattern, weights = weights, prop = 0.5, mech = mech)
  
  ampdata <- dat_misspattern
  ampdata$premedicalcost <- md1$amp$premedicalcost
  ampdata[,which(cols_prevDMTefficacy)] <- md2$amp[,which(cols_prevDMTefficacy)]
  ampdata[,which(cols_numSymptoms)] <- md2$amp[,which(cols_numSymptoms)]
  ampdata$prerelapse_num <- md3$amp$prerelapse_num
  
  ## Collapse dummies back into categorical variables
  ampdata$treatment <- factor(ampdata$treatment, levels = c(1,0), labels = c("DMF", "TERI"))
  ampdata <- ampdata %>% mutate(prevDMTefficacy = ifelse(prevDMTefficacy_None == 1, "None",
                                                        ifelse(prevDMTefficacy_Low_efficacy == 1, "Low_efficacy", "Medium_high_efficacy")),
                                numSymptoms = ifelse(numSymptoms_0 == 1, "0", ifelse(numSymptoms_1 == 1, "1", ">=2")))

  ampdata$prevDMTefficacy <- factor(ampdata$prevDMTefficacy,
                                   levels = c("None", "Low_efficacy", "Medium_high_efficacy"),
                                   labels = c("None", "Low_efficacy", "Medium_high_efficacy"))
  ampdata$numSymptoms <- factor(ampdata$numSymptoms,
                                levels = c("0", "1", ">=2"),
                                labels = c("0", "1", ">=2"))
  
  # Get rid of dummies
  cols_prevDMTefficacy <- which(grepl( "prevDMTefficacy_", colnames(ampdata), fixed = TRUE))
  cols_numSymptoms <- which(grepl( "numSymptoms_", colnames(ampdata), fixed = TRUE))
  
  return(ampdata %>% dplyr::select(-starts_with("prevDMTefficacy_")) %>% dplyr::select(-starts_with("numSymptoms_")))
}



#F3. Function to  get treatment effect estimands after mice imputation----

getmicest <- function(data,estimandv,CC,Tform,approachv){
  
  if (estimandv=="ATE"){ # for ATE
    methodv <- "full"
    replacev <- FALSE
    } else {             # for ATT
      methodv <- "nearest"
      replacev <- TRUE
    }

 
  # Propensity score model
  formula.full <- DMF ~ age + female + prevDMTefficacy + premedicalcost + prerelapse_num
  formula.mi   <- DMF ~ age + female + pde.ind + pde.mi + pmc.ind + pmc.mi + prn.ind + prn.mi
  formula.mic  <- DMF ~ age + female + pde.ind + prevDMTefficacy + pmc.ind + premedicalcost + prn.ind + prerelapse_num
  
  if(Tform == 1){
    formula <- formula.full
    }else if(Tform == 2){
      formula <- formula.mi
    }else {
        formula <- formula.mic}
  
   # Matching based on PS model
   matched.datasets <- matchthem(formula,
                                datasets = data,
                                approach = approachv,
                                method = methodv,
                                caliper = 0.2,
                                family = binomial,
                                estimand =estimandv,
                                distance = "glm",
                                link = "logit",
                                replace = replacev) 
  
   matched.results <- summary(pool(with(matched.datasets,
                                       svyglm(y ~ DMF + offset(log(years)), family = poisson(link = "log")),
                                       cluster = TRUE)),
                                       conf.int = TRUE)
  
   weighted.datasets <- weightthem(formula,
                                  datasets = data,
                                  approach = approachv,
                                  method = 'ps',
                                  estimand = estimandv)
  
   weighted.results <- summary(pool(with(weighted.datasets,
                                        svyglm(y ~ DMF + offset(log(years)), family = poisson(link = "log")),
                                        cluster = TRUE)),conf.int=TRUE)
  
  results <- setDT(rbind(matched.results,weighted.results))[c(2,4),c(2,3,7,8)]
  results[, method := c("Matching","IPTW")]
  results[, estimand := estimandv]
  return(results)
}


#F4. Function to estimate the treatment effect in a complete dataset----


getest <- function(data, 
                   estimandv = "ATE", # Estimate the ATE or ATT  
                   Tform, # PS model formula
                   CC = FALSE, # use the complete case dataset,
                   approachv = NULL){
  
  # Prepare output
  result <- data.frame("method" = character(),
                       "estimand" = character(),
                       "estimate" = numeric(),
                       "std.error" = numeric(),
                       "LCI" = numeric(),
                       "UCI" = numeric())
  
  if (CC) { # Get Complete Case dataset
    data <- data[complete.cases(data), ]
  }
  
    data[, DMF := as.numeric(treatment == "DMF")]

    if (estimandv == "ATE"){ # for ATE
      methodv <- "full"
      replacev <- FALSE
    } else {             # for ATT
      methodv <- "nearest"
      replacev <- TRUE
    }
  
    # Propensity score model
    if (Tform == 1) {
      formula <- DMF ~ age + female + prevDMTefficacy + premedicalcost + prerelapse_num
    } else if (Tform == 2) {
      formula <- DMF ~ age + female + pde.ind + pde.mi + pmc.ind + pmc.mi + prn.ind + prn.mi
    } else {
      formula <-DMF ~ age + female + pde.ind + prevDMTefficacy + pmc.ind + premedicalcost + prn.ind + prerelapse_num
    }
  
  # Apply Matching
  mout <- matchit(formula, 
                  data = data,
                  family = binomial,
                  method = methodv,
                  caliper = 0.2,
                  std.caliper = TRUE,
                  estimand = estimandv,
                  distance = "glm",
                  link = "logit",
                  replace = replacev) # we apply with replacement here due to small number of treated patients

  if (estimandv == "ATT") {
     mdata <- as.data.table(get_matches(object = mout, id = "matching_id"))
     assign("mdata", mdata, envir = .GlobalEnv)
     match_mod <- glm("y ~ DMF + offset(log(years))",
                      family = poisson(link = "log"),
                      data = mdata)
     # Estimate cluster-robust standard error
     tx_var <- vcovCL(match_mod, cluster = ~ subclass + matching_id, sandwich = TRUE)
     
    } else if (estimandv == "ATE") {
      mdata <- as.data.table(MatchIt::match.data(object = mout))
      assign("mdata", mdata, envir = .GlobalEnv)
      match_mod <- glm("y ~ DMF + offset(log(years))",
                       family = poisson(link = "log"),
                       data = mdata,
                       weights = weights)
      # Estimate robust variance-covariance matrix
      tx_var <- vcovCL(match_mod, cluster = ~ subclass, sandwich = TRUE) 
    }
  
  result <- result %>% add_row(method = "Matching", 
                               estimand = estimandv,
                               estimate = coef(match_mod)["DMF"],
                               std.error = sqrt(tx_var["DMF", "DMF"]),
                               LCI = coef(match_mod)["DMF"] + qnorm(0.025)*sqrt(tx_var["DMF", "DMF"]),
                               UCI = coef(match_mod)["DMF"] + qnorm(0.975)*sqrt(tx_var["DMF", "DMF"]))
  
  # Apply IPTW
  wout  <- weightit(formula, data = data, estimand = estimandv, method = "ps")
  rhcSvy  <- svydesign(ids = ~ 1, data = data, weights = ~ wout$weights)
  ipw_mod <- svyglm("y ~ DMF + offset(log(years))",
                    family  =  poisson(link = "log"),
                    design = rhcSvy)
  
  # Note: svyglm always returns 'model-robust' standard errors; 
  # the Horvitz-Thompson-type standard errors used everywhere in the survey 
  # package are a generalisation of the model-robust 'sandwich' estimators. 
  
  result <- result %>% add_row(method = "IPTW", 
                               estimand = estimandv,
                               estimate = coef(ipw_mod)["DMF"],
                               std.error =  sqrt(diag(vcov(ipw_mod))["DMF"]),
                               LCI = coef(ipw_mod)["DMF"] + qnorm(0.025)*sqrt(diag(vcov(ipw_mod))["DMF"]),
                               UCI = coef(ipw_mod)["DMF"] + qnorm(0.975)*sqrt(diag(vcov(ipw_mod))["DMF"]))
  
  result<-rename(result, `2.5 %`=LCI, `97.5 %`=UCI)
  
  return(result)
}

#F5. Function to impute mice separated by groups ----
separate_mice <- function(data, form_y, method, m = 5) {
  phr_DMF1  <- subset(data, DMF == 1)
  phr_DMF0  <- subset(data, DMF == 0)
  DMF1_imps <- mice(phr_DMF1, m = m, form = form_y, method = method)
  DMF0_imps <- mice(phr_DMF0, m = m, form = form_y, method = method)
  imps <- rbind(DMF1_imps, DMF0_imps)
  return(imps)
}


#F6. Function to calculate all estimates at once ----


allest <- function(homo, hete, functionv, Tformv, CCv, typev, approachv = NULL){
  homoATE <- rbindlist(lapply(homo, FUN = functionv, estimandv = "ATE", Tform = Tformv, CC = CCv, approach = approachv))
  homoATT <- rbindlist(lapply(homo, FUN = functionv, estimandv = "ATT", Tform = Tformv, CC = CCv, approach = approachv))
  heteATE <- rbindlist(lapply(hete, FUN = functionv, estimandv = "ATE", Tform = Tformv, CC = CCv, approach = approachv))
  heteATT <- rbindlist(lapply(hete, FUN = functionv, estimandv = "ATT", Tform = Tformv, CC = CCv, approach = approachv))
  comb <- setDT(rbind(homoATE, homoATT, heteATE, heteATT))
  comb[, Scenario  := rep(c("MCAR","MCAR","MAR","MAR","MART","MART","MART2","MART2","MNAR","MNAR"),4)]
  comb[, Treatment := c(rep("Homogeneus",20),rep("Heterogeneus",20))]
  comb[, type      := typev]
}



#F7.  Density plot by treatment
plot_den <- function(datap,var,varlab) {
  mu<-datap[,.(grp.mean=mean(get(var))),treatment]
  ggplot(datap, aes(x=get(var), color=treatment)) +
    geom_density()+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=treatment),
               linetype="dashed")+
    labs(x=varlab, y = "Density",colour="Treatment")+
    scale_color_brewer(palette="Accent") + theme_minimal()+theme(legend.position="top")
}  

#F8.  Proportion plot by treatment
plot_count <- function(datap,var,varlab) {
  counts <-datap[,.(grp.count=.N),list(treatment,get(var))]
  counts[, percentage := prop.table(grp.count), treatment] 
  ggplot( counts,aes(x = get, y = percentage, fill = treatment, colour = treatment)) + 
    geom_bar(stat = "identity", position = "dodge", alpha = 0.5)  +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_fill_brewer(palette = "Accent") +
    scale_color_brewer(palette = "Accent")+
    labs(y ="Percentage (%)",x= varlab, colour="Treatment",fill="Treatment")+theme(legend.position="top")
} 

#F9.  Annual relapse rate plot
Relapserate_plot <-function(datahom,datahet){
  datall<-setDT(rbind(datahom,datahet))
  datall[,HTE:=c(rep("No HTE",nrow(datahom)),rep("High HTE",nrow(datahet)))]
  datall<-datall[!(HTE=="No HTE"&treatment=="TERI")]
  datall[,treatment:=ifelse(treatment=="TERI","TERI",ifelse(HTE=="No HTE"&treatment=="DMF","DMF(homogeneous)","DMF(heterogeneous)"))]
  datall[,treatment:=factor(treatment,levels=c("TERI","DMF(homogeneous)","DMF(heterogeneous)"))]
  datall[,AR:=y/years]
  datall[,Iscore:=as.factor(Iscore)]
  levels(datall$Iscore)<-c("High DMF","Moderate DMF","Neutral","Moderate TERI","High TERI")
  datallgroup<-datall[,.(grp.mean=mean(AR),sd=sd(AR),n=.N),list(treatment,Iscore)]
  datallgroup[,LCI:=grp.mean-1.96*sd/sqrt(n)]
  datallgroup[,UCI:=grp.mean+1.96*sd/sqrt(n)]
  
  pd <- position_dodge(.1)
  ggplot(datallgroup,aes(x=Iscore,y=grp.mean,group=treatment,colour=treatment))+
    geom_errorbar(aes(ymin=LCI, ymax=UCI), 
                  width=.1, position=pd) +
    geom_line(position=pd) + 
    geom_point(position=pd, size=2)+
    labs(y="Annual post relapse rate",x="Effect treatment group (Iscore)",color="")+
    scale_colour_manual(values = c("#BEAED4","#7FC97F","#FDC086"))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45),legend.position = "top")
}

formatMSdata <- function(data) {
  data$female <- factor(data$female, levels = c(1,0), labels = c("Yes", "No"))
  data$prevDMTefficacy<- factor(data$prevDMTefficacy, 
                                      levels = c("None", "Low_efficacy", "Medium_high_efficacy"), 
                                      labels = c("None", "Low", "Medium or High"))
  #sim_NRS$PRMSGR <- factor(sim_NRS$PRMSGR, levels = c(1,0), labels = c("Yes", "No"))
  
  label(data$age)       <- "Age"
  label(data$female)       <- "Female Sex"
  label(data$premedicalcost)       <- "Prior medical costs"
  label(data$prevDMTefficacy)       <- "Efficacy of previous DMT"
  label(data$numSymptoms)       <- "Number of prior symptoms"
  label(data$prerelapse_num)       <- "Number of prior relapses"
  
  units(data$age)       <- "years"
  
  return(data)
}
