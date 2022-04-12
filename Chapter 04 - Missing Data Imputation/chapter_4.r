## code to prepare `countExample` dataset

library(tidyverse, warn.conflicts = FALSE)
library(MASS)
library(truncnorm)
library(magrittr)
library(reshape2)
library(fastDummies)
library(dplyr)
library(optmatch) #match
library(MatchIt)
library(WeightIt) #IPW
library(cobalt) # imbalance covariates
library(mice)
library(data.table)
library(PSweight)
library(sandwich) 
library(lmtest)
library(MatchThem)
library(survey)
library(missForest)
library(ggplot2)

#F0. Function to simulate MS data
simcountdata <- function(n, 
                         seed = NA,
                         beta = c(-0.5, -0.25, 0, 0.25, 0.5),
                         beta.x = c(-1.54, -0.01, 0.06, 0.25, 0.5, 0.13, 0.0000003),
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
  
  if (percentiles[1] != 0 | percentiles[length(percentiles)] != 1 | length(percentiles) != 6){message("Wrong values of percentiles!")}
  
  # Create an empty shell
  ds <- data.frame(matrix(NA, nrow = n, ncol = 10))
  colnames(ds) <- c("treatment", "ageatindex_centered", "female", "prerelapse_num",
                    "prevDMTefficacy", "premedicalcost", "numSymptoms",
                    "postrelapse_num", "finalpostdayscount", "group")
  
  # Define X, A, and time
  ds %<>%
    mutate(trt =                 rbinom(n = n, size = 1, prob = 0.75),
           treatment =           as.factor(ifelse(trt == 1, "DMF", "TERI")),
           female =              rbinom(n = n, size = 1, prob = 0.75),
           ageatindex_centered = round(rtruncnorm(n, a = 18, b = 64, mean = 48, sd = 12), 0) - 48, # rounded to integers
           prerelapse_num =      rpois(n = n, lambda = 0.44),
           prevDMTefficacy =     sample(x = c("None", "Low efficacy", "Medium and high efficacy"),
                                        size = n, replace = TRUE, prob = c(0.45, 0.44, 0.11)),
           premedicalcost =      pmin(round(exp(rnorm(n = n, mean = 8.9, sd = 1.14)), 2), 600000), # rounded to 2 decimal points
           numSymptoms =         sample(x = c("0", "1", ">=2"),
                                        size = n, replace = TRUE, prob = c(0.67, 0.24, 0.09)), # nuisance variable; do not include in the score
           finalpostdayscount =  ceiling(rgamma(n = n, shape = 0.9, scale = 500)), # rounded up to integers
           finalpostdayscount =  ifelse(finalpostdayscount > 2096, 2096, finalpostdayscount), # truncate at the max follow up day, 2096
           finalpostdayscount =  ifelse((finalpostdayscount > 2090) & (runif(1, 0, 1) < .5), 29, finalpostdayscount), # mimic the 1 month peak;  move roughly half of the large values to 29
           group =              "Simulated")
  
  # Define Y
  xmat.score <- model.matrix(~ ageatindex_centered + female + prerelapse_num + prevDMTefficacy + premedicalcost, ds) %>% as.matrix()
  gamma <- matrix(c(-0.33, # Intercept
                    -0.001, # Age
                    0.05, # female
                    -0.002, # prerelapse_num
                    0.33, 0.02, # Medium/high and none DMT efficacy
                    -0.0000005), nrow = 7) # premedicalcost
  ds <- ds %>% mutate(score = exp(xmat.score %*% gamma),
                      Iscore = cut(score, quantile(score, percentiles), include.lowest = T, labels = seq(1, 5)))
  
  xmat.rate <- model.matrix(~ ageatindex_centered + female + prerelapse_num + prevDMTefficacy + premedicalcost +
                              Iscore + trt + trt*Iscore, ds) %>% as.matrix()
  
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
  ds <- ds %>% mutate(postrelapse_num = rnegbin(n = n, mu = rate*finalpostdayscount/365.25, theta = 3))
  
  return(list(data = ds, betas = betas, percentiles = percentiles, rate=rate))
}


#F1. Function to generate complete data ----
generate_data <- function(n = 10000, 
                          beta = c(-0.2, -0.2, -0.2, -0.2, -0.2),
                          seed = NA) {
  
  data <- simcountdata(n = n, 
                       seed = seed,
                       beta = beta,
                       #beta=c(-0.2, -0.2, -0.2, -0.2, -0.2) 
                       #beta = c(log(0.4), log(0.5), log(1), log(1.1), log(1.2)), 
                       # this beta is medium level of heterogeneity (see other levels in "Simulation design.pptx")
                       # other options:
                       #  none - beta = c(-0.2. -0.2, -0.2, -0.2, -0.2)
                       #  low - beta = c(log(0.7), log(0.75), log(1), log(1.05), log(1.1))
                       #  high - beta = c(log(0.3), log(0.5), log(1), log(1.1), log(1.5))
                       beta.x = c(-1.54, -0.01, 0.06, 0.25, 0.5, 0.13, 0.0000003)
  )$data %>% dplyr::rename(previous_treatment = prevDMTefficacy,
                           age = ageatindex_centered,
                           y = postrelapse_num,
                           years = finalpostdayscount,
                           previous_number_relapses = prerelapse_num,
                           previous_number_symptoms = numSymptoms,
                           previous_cost = premedicalcost) %>%
    mutate(previous_treatment = factor(previous_treatment, labels = c("drugA", "drugB", "drugC")),
           previous_number_symptoms = factor(previous_number_symptoms, labels = c("0", "1", ">=2")),
           years = years / 365.25,
           age = age + 48,
           DMF = as.factor(ifelse(treatment == "DMF",1,0)))%>%
    dplyr::select(age, female, previous_treatment, previous_cost, previous_number_symptoms, previous_number_relapses, DMF, y, years,Iscore)
  return(data)
}


#F1.1 Function to transform database with dummy to categorical variables----
databack<-function(data){
  data<-setDT(data)
  data[,female:=as.factor(female)]
  data[,previous_treatment:=as.factor(ifelse(previous_treatmentdrugA==1,"drugA",ifelse(previous_treatmentdrugB==1,"drugB",ifelse(is.na(previous_treatmentdrugA),NA,"drugC"))))]
  data[,previous_number_symptoms:=as.factor(ifelse(previous_number_symptoms1==1,"1",ifelse(previous_number_symptoms..2==1,">=2",ifelse(is.na(previous_number_symptoms1),NA,"0"))))]
  data[,Iscore:=as.factor(ifelse(Iscore2==1,"2",ifelse(Iscore3==1,"3",ifelse(Iscore4==1,"4",ifelse(Iscore5==1,"5","1")))))]
  data[,DMF:=as.factor(DMF1)]
  data[,DMF1:=NULL]
  data[,previous_number_symptoms:=factor(previous_number_symptoms,levels=c("0","1",">=2"))]
  data[,previous_treatmentdrugA:=NULL]
  data[,previous_treatmentdrugB:=NULL]
  data[,previous_number_symptoms1:=NULL]
  data[,previous_number_symptoms..2:=NULL]
  data[,Iscore2:=NULL]
  data[,Iscore3:=NULL]
  data[,Iscore4:=NULL]
  data[,Iscore5:=NULL]
  data[,prer.ind:= as.factor(as.numeric(is.na(previous_number_relapses)==FALSE))] 
  data[,prer.mi:=ifelse(is.na(previous_number_relapses),0,previous_number_relapses)]
  data[,prco.ind:= as.factor(as.numeric(is.na(previous_cost)==FALSE))] 
  data[,prco.mi:=ifelse(is.na(previous_cost),0,previous_cost)]
  data[,prns.ind:= as.factor(as.numeric(is.na(previous_number_symptoms)==FALSE))] 
  data[,prns.mi:=ifelse(is.na(previous_number_symptoms),"na",previous_number_symptoms)]
  return(as.data.frame(data))
}

#F2.Function to generate missing data ---
getmissdata <- function(data){
  data<-data.frame(model.matrix(~.-1,data))
  data$previous_treatmentdrugC<-NULL
  prer_mcar<-ampute(data, patterns=c(1,1,1,1,1,1,1,0,1,1,1,1,1,1,1),prop=0.5, mech="MCAR")$amp
  prer_mar<-ampute(data, patterns = c(1,1,1,1,1,1,1,0,1,1,1,1,1,1,1),weights=c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0), prop = 0.5, mech = "MAR")$amp
  prer_mart<-ampute(data, patterns = c(1,1,1,1,1,1,1,0,1,1,1,1,1,1,1),weights=c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0), prop = 0.5, mech = "MAR",type="LEFT")$amp
  prer_mnar<-ampute(data, patterns = c(1,1,1,1,1,1,1,0,1,1,1,1,1,1,1),weights=c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0), prop = 0.5, mech = "MNAR",type="LEFT")$amp
  cost_mcar<-ampute(data, patterns=c(1,1,1,1,0,1,1,1,1,1,1,1,1,1,1),prop=0.3, mech="MCAR")$amp
  presym_mcar<-ampute(data, patterns=c(1,1,1,1,1,0,0,1,1,1,1,1,1,1,1),prop=0.15, mech="MCAR")$amp

  
  #MCAR
  ampdata1<-prer_mcar
  ampdata1[,c("previous_number_symptoms1","previous_number_symptoms..2")]<-presym_mcar[,c("previous_number_symptoms1","previous_number_symptoms..2")]
  ampdata1$previous_cost<-cost_mcar$previous_cost

  #MAR
  ampdata2<-prer_mar
  ampdata2[,c("previous_number_symptoms1","previous_number_symptoms..2")]<-presym_mcar[,c("previous_number_symptoms1","previous_number_symptoms..2")]
  ampdata2$previous_cost<-cost_mcar$previous_cost
  
  
  #MAR-y(treatment)
  ampdata3<-prer_mart
  ampdata3[,c("previous_number_symptoms1","previous_number_symptoms..2")]<-presym_mcar[,c("previous_number_symptoms1","previous_number_symptoms..2")]
  ampdata3$previous_cost<-cost_mcar$previous_cost

  
  #MAR-y(MNAR)
  ampdata4<-prer_mnar
  ampdata4[,c("previous_number_symptoms1","previous_number_symptoms..2")]<-presym_mcar[,c("previous_number_symptoms1","previous_number_symptoms..2")]
  ampdata4$previous_cost<-cost_mcar$previous_cost

  
  return(list(ampdata1=databack(ampdata1),
              ampdata2=databack(ampdata2),
              ampdata3=databack(ampdata3),
              ampdata4=databack(ampdata4)))
}


#F3. Function to  get treatment effect estimands after mice imputation----
getmicest<-function(data,estimandv,CC,Tform,approachv,replacev=FALSE){

  methodv <- ifelse(estimandv=="ATE","full","nearest")
  
  formula.full <-DMF ~ age + female + previous_treatment + previous_cost+  previous_number_symptoms + previous_number_relapses
  formula.mi <- DMF ~ age + female + previous_treatment + previous_cost+ previous_number_symptoms + prer.ind + prer.mi+ prco.ind + prco.mi + prns.ind + prns.mi
  formula.mic <-DMF ~ age + female + previous_treatment + previous_cost+ previous_number_symptoms + prer.ind + previous_number_relapses + prco.ind + previous_cost + prns.ind + previous_number_symptoms
  
  if(Tform==1){
    formula=formula.full
    }else if(Tform==2){
      formula=formula.mi
      }else {
        formula=formula.mic}
  
  matched.datasets <- matchthem(formula,
                                datasets = data,
                                approach = approachv,
                                method = methodv,
                                caliper = 0.2,
                                family = binomial,
                                estimand =estimandv,
                                distance = "glm",
                                link = "logit",
                                replace = replacev) #change that!!! cluster standard errors check !!!
  
  matched.results <- summary(pool(with(matched.datasets,
                                       svyglm(y ~ DMF+offset(log(years)), family = poisson(link = "log")),
                                       cluster = TRUE)),conf.int = TRUE)
  
  weighted.datasets <- weightthem(formula,
                                  datasets = data,
                                  approach = approachv,
                                  method = 'ps',
                                  estimand = estimandv)
  
  weighted.results <- summary(pool(with(weighted.datasets,
                                        svyglm(y ~ DMF+offset(log(years)), family = quasipoisson(link = "log")),
                                        cluster=TRUE)),conf.int=TRUE)
  
  results<-setDT(rbind(matched.results,weighted.results))[c(2,4),c(2,3,7,8)]
  results[,method:=c("Matching","IPW")]
  results[,estimand:=estimandv]
  return(results)
}




#F4. Function to estimate the treatment effect in a complete dataset----
getest <- function(data, 
                   estimandv = "ATE", # Estimate the ATE or ATT ? 
                   Tform, 
                   CC = FALSE, # Omit incomplete cases?
                   approachv, 
                   replacev = FALSE # perform matching with or without replacement?
                   ){
  if (CC) {
    data <- data[complete.cases(data), ]
  }
  
  formula.full <- DMF ~ age + female + previous_treatment + previous_cost+  previous_number_symptoms + previous_number_relapses
  formula.mi <- DMF ~ age + female + previous_treatment + prer.ind + prer.mi+ prco.ind + prco.mi + prns.ind + prns.mi
  formula.mic <-DMF ~ age + female + previous_treatment + prer.ind + previous_number_relapses + prco.ind + previous_cost + prns.ind + previous_number_symptoms
  
  if(Tform==1){
    formula=formula.full
  }else if(Tform==2){
    formula=formula.mi
  }else {
    formula=formula.mic}
  
  methodv <- ifelse(estimandv == "ATE", "full", "nearest")

  # Apply matching
  mout <- matchit(formula, 
                  data = data,
                  family = binomial,
                  method = methodv,
                  caliper = 0.2,
                  std.caliper = TRUE,
                  estimand = estimandv,
                  distance = "glm",
                  link = "logit",
                  replace = replacev)
  
  # Step 2: retrieve matched sample
  if (estimandv == "ATE"){
    mdata <- match.data(mout)
  } else{
    mdata <- get_matches(mout)
  }
  
  match_mod <- glm("y ~ DMF + offset(log(years))", 
                   family = poisson(link = "log"),
                   data = mdata)
  
  ## I made some changes here to ensure you obtain a (cluster-)robust standard error
  # TODO
  
  match_fit <- summary(match_mod)$coefficients[, 1]["DMF1"]
  match_se <- summary(match_mod)$coefficients[, 2]["DMF1"]
  match_res<-c(match_fit,match_se )
  
  # Apply IPW
  wout <- weightit(formula,data = data, estimand = estimandv, method = "ps")
  data$ipw<-wout$weights
  ipw_mod <- glm("y ~ DMF + offset(log(years))", 
                 family = poisson(link = "log"),
                 data = data,
                 weights = ipw)
  ipw_fit <- summary(ipw_mod)$coefficients[, 1]["DMF1"]
  ipw_se <- coeftest(ipw_mod, vcov = vcovHC)["DMF1","Std. Error"] #Robust estimate 
  ipw_res<-c(ipw_fit,ipw_se )
  res<-rbind(match_res,ipw_res)
  colnames(res)<-c("estimate","std.error")
  res<-as.data.table(res)
  za <- qnorm(1 - 0.05 / 2) # for 95% CIs
  res[,"2.5 %":= estimate-za*std.error]
  res[,"97.5 %":= estimate+za*std.error]
  res[,method:= c("Matching","IPW")]
  res[,estimand:= estimandv]
}


#F5.Function to impute mice separated by groups ----
separate_mice <- function(data, form_y, method) {
  phr_DMF1  <- subset(data, DMF == 1)
  phr_DMF0  <- subset(data, DMF == 0)
  DMF1_imps <- mice(phr_DMF1,m = 5, form = form_y, method = method)
  DMF0_imps <- mice(phr_DMF0, m = 5, form = form_y, method = method)
  imps <- rbind(DMF1_imps, DMF0_imps)
  return(imps)
}

#F6. Function to calculate all estimates at once ----

allest<-function(homo=mdata_noHTE,hete=mdata_hiHTE,functionv,Tformv,CCv,typev,approachv=NULL,replacev=FALSE){
  homoATE<-rbindlist(lapply(homo,FUN=functionv,estimandv="ATE",Tform=Tformv,CC=CCv,approach=approachv,replacev=replacev))
  homoATT<-rbindlist(lapply(homo,FUN=functionv,estimandv="ATT",Tform=Tformv,CC=CCv,approach=approachv,replacev=replacev))
  heteATE<-rbindlist(lapply(hete,FUN=functionv,estimandv="ATE",Tform=Tformv,CC=CCv,approach=approachv,replacev=replacev))
  heteATT<-rbindlist(lapply(hete,FUN=functionv,estimandv="ATT",Tform=Tformv,CC=CCv,approach=approachv,replacev=replacev))
  comb<-setDT(rbind(homoATE,homoATT,heteATE,heteATT))
  comb[,Scenario:=rep(c("MCAR","MCAR","MAR","MAR","MART","MART","MNAR","MNAR"),4)]
  comb[,Treatment:=c(rep("Homogeneus",16),rep("Heterogeneus",16))]
  comb[,type:=typev]
}
