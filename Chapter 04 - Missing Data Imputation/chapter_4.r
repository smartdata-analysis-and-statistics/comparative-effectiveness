rm(list=ls())
library(data.table)
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
  #(Intercept) female ageatindex_centered prerelapse_num prevDMTefficacyLow prevDMTefficacyMedium_high_efficacy     premedicalcost numSymptoms0 numSymptoms1 finalpostdayscount                                         
  #c(1.38,0.5,-1.2,1.6,1.2,2.5,-0.001,0.5,1.2,0)
  XB <- model.matrix(~.,ds)%*% c(1.22,0.3,-0.1,0.2, 0.35,0.7,-0.00005,0.17,0.02,0)# ~75% people allocated in DMF arm based on (age,female,prerelapse_num,DMT efficacy,costs,numSymptoms)
  pi <- exp(XB)/(1+exp(XB))
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
                    beta.x[5], beta.x[6], # Medium/high and none DMT efficacy
                    beta.x[7], # premedicalcost
                    # Previous beta.x: -1.54 Intercept, -0.01 Age, 0.06 female, 0.47 prerelapse_num, 0.68, 0.13 Medium/high and none DMT efficacy, 0.000003 premedicalcost
                    0, 0, 0, 0, # Iscore categories 2:5 (reference group is Iscore1, high responders to DMF)
                    beta[1], beta[2] - beta[1], beta[3] - beta[1], beta[4] - beta[1], beta[5] - beta[1]))
  rate <-  exp(xmat.rate %*% betas)
  ds[, y := rnegbin(n = n, mu = rate*finalpostdayscount/365.25, theta = 3)] # post treatment number of relapses
 
  ds[, Iscore := factor (Iscore, labels= c("High A1","Moderate A1","Neutral","Moderate A0","High A0"))]
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


getmissdata <- function(data){
  set.seed(12345)
  data<-as.data.table(data)
  data[,DMF := as.numeric(treatment == "DMF")]
  # missingness on MS tracking variables (prevDMTefficacy, numSymptoms)
  data[,pmar:=-age*0.001+0.05*female]
  data[,pmart:=+0.075*DMF*(-age*0.001+0.05*female)-age*0.001+0.05*female]
  data[,pmart2:=+0.075*DMF*(-age*0.001+0.05*female)-age*0.001+0.05*female-0.005*y]
  data[,pmnar:=-age*0.001+0.05*female+1.2*prerelapse_num]
  amp_initial <- ampute(data)
  
  sick_patterns <- rbind(amp_initial$patterns[c(3,5),],c(1,1,0,1,0,1,1,1,1,1,1,1,1,1))
  sick_weights <- ampute(data, patterns = sick_patterns,prop=0.3, mech="MAR")$weights
  sick_weights[,"prevDMTefficacy"] <- c(0,0,0)
  sick_weights[,"numSymptoms"] <- c(0,0,0)
  sick_weights[,"age"] <- c(-0.2,-0.2,-0.2)
  sick_weights[,"female"] <- c(0.3,0.3,0.3)
 
  sick_mcar <- ampute(data, patterns = sick_patterns, prop = 0.3, mech = "MCAR")$amp
  sick_mar <- ampute(data, patterns = sick_patterns, weights = sick_weights, prop = 0.3, mech = "MAR")$amp
  
  # missingness on administrative variables (premedical cost)
  cost_mcar<-ampute(data, patterns = c(1,1,1,0,1,1,1,1,1,1,1,1,1,1,1), prop = 0.10, mech = "MCAR")$amp
  
  # missingness on prereleapse number
  prer_mcar <- ampute(data, patterns = c(1,1,1,1,1,0,1,1,1,1,1,1,1,1,1), prop = 0.5, mech = "MCAR")$amp
  prer_mar <- ampute(data, patterns =  c(1,1,1,1,1,0,1,1,1,1,1,1,1,1,1), weights = c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0), prop = 0.5, mech = "MAR")$amp
  prer_mart <- ampute(data, patterns = c(1,1,1,1,1,0,1,1,1,1,1,1,1,1,1), weights = c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0), prop = 0.5, mech = "MAR")$amp
  prer_mart2 <- ampute(data, patterns =c(1,1,1,1,1,0,1,1,1,1,1,1,1,1,1), weights = c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0), prop = 0.5, mech = "MAR")$amp
  prer_mnar <- ampute(data, patterns = c(1,1,1,1,1,0,1,1,1,1,1,1,1,1,1), weights = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1), prop = 0.5, mech = "MAR")$amp

  
  #MCAR
  ampdata1 <- sick_mcar
  ampdata1$premedicalcost <- cost_mcar$premedicalcost
  ampdata1$prerelapse_num <- prer_mcar$prerelapse_num
  
  
  #MAR
  ampdata2 <- sick_mar
  ampdata2$premedicalcost <- cost_mcar$premedicalcost
  ampdata2$prerelapse_num <- prer_mar$prerelapse_num
  
  
  #MAR-y(treatment)
  ampdata3 <- sick_mar
  ampdata3$premedicalcost <- cost_mcar$premedicalcost
  ampdata3$prerelapse_num <- prer_mart$prerelapse_num
  
  
  #MAR-y(MNAR)
  ampdata4 <- sick_mar
  ampdata4$premedicalcost <- cost_mcar$premedicalcost
  ampdata4$prerelapse_num <- prer_mart2$prerelapse_num
  
  #MAR-y(MNAR)
  ampdata5 <- sick_mar
  ampdata5$premedicalcost <- cost_mcar$premedicalcost
  ampdata5$prerelapse_num <- prer_mnar$prerelapse_num
  
  
  return(list(ampdata1 = databack(ampdata1),
              ampdata2 = databack(ampdata2),
              ampdata3 = databack(ampdata3),
              ampdata4 = databack(ampdata4),
              ampdata5 = databack(ampdata5)
              ))
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
  results[, method := c("Matching","IPW")]
  results[, estimand := estimandv]
  return(results)
}


#F4. Function to estimate the treatment effect in a complete dataset----


getest <- function( data, 
                    estimandv = "ATE", # Estimate the ATE or ATT  
                    Tform, # PS model formula
                    CC = FALSE, # use the complete case dataset,
                    approachv){
  
  if (CC) { # Get Complete Case dataset
    data <- data[complete.cases(data), ]
  }
  
    data[, DMF := as.numeric(treatment == "DMF")]

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
    
  
  # Apply Matching
    
  mout <- MatchIt::matchit(formula, 
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
     mdata <- as.data.table(MatchIt::get_matches(object = mout, id = "matching_id"))
     assign("mdata", mdata, envir = .GlobalEnv)
     match_mod <- glm("y ~ DMF + offset(log(years))",
                      family = poisson(link = "log"),
                      data = mdata)
     # Estimate cluster-robust standard error
     tx_var <- sandwich::vcovCL(match_mod, cluster = ~ subclass + matching_id, sandwich = TRUE)
     
    } else if (estimandv == "ATE") {
      mdata <- as.data.table(MatchIt::match.data(object = mout))
      assign("mdata", mdata, envir = .GlobalEnv)
      match_mod <- glm("y ~ DMF + offset(log(years))",
                       family = poisson(link = "log"),
                       data = mdata,
                       weights = weights)
      # Estimate robust variance-covariance matrix
      tx_var <- sandwich::vcovCL(match_mod, cluster = ~ subclass, sandwich = TRUE) 
    }
  
    match_fit <- coef(match_mod)["DMF"]
    match_se  <- sqrt(tx_var["DMF", "DMF"])
    match_res <- c(match_fit, match_se)
  
  # Apply IPW
  wout  <- weightit(formula, data = data, estimand = estimandv, method = "ps")
  rhcSvy  <- svydesign(ids = ~ 1, data = data, weights = ~ wout$weights)
  ipw_mod <- svyglm("y ~ DMF + offset(log(years))",
                    family  =  poisson(link = "log"),
                    design = rhcSvy)
  
  ipw_fit <- coef(ipw_mod)["DMF"]
  ipw_se  <- sqrt(diag(vcov(ipw_mod))["DMF"])
  ipw_res <- c(ipw_fit,ipw_se )
  
  # Combine Matching and IPW results
  res <- rbind(match_res,ipw_res)
  colnames(res) <- c("estimate","std.error")
  res <- as.data.table(res)
  za <- qnorm(1 - 0.05 / 2) # for 95% CIs
  res[, "2.5 %" := estimate-za*std.error]
  res[, "97.5 %" := estimate+za*std.error]
  res[, method := c("Matching","IPW")]
  res[, estimand := estimandv]
  
  return(res)
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

