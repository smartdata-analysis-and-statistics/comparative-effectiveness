# Load packages ----
# Package names
packages <- c("data.table","dplyr","tidyr","kableExtra","table1", # data formatting and plotting
              "MASS", "truncnorm", # "data simulation
              "mixgb", # modelling
              "ggpubr", "naniar", # data visualization
              "optmatch", "MatchIt", "MatchThem", "WeightIt","PSweight", "sandwich", "cobalt","survey", #PS estimation
              "marginaleffects", #
              "mice", "missForest","ggmice", # Multiple imputation
              "ranger", "mixgb")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Data generation ----
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
#'              # the treatment effect options are:
#                  none - beta = c(-0.2. -0.2, -0.2, -0.2, -0.2)
#                  medium - beta= c(log(0.4), log(0.5), log(1), log(1.1), log(1.2))
#                  low - beta = c(log(0.7), log(0.75), log(1), log(1.05), log(1.1))
#                  high - beta = c(log(0.3), log(0.5), log(1), log(1.1), log(1.5))
#' @param beta.x coefficients for main effects of other covariates in the rate; vector of length 7
#'               beta.x[1] (intercept)
#'               beta.x[2]*ageatindex_centered
#'               beta.x[3]*female
#'               beta.x[4]*prerelapseNum
#'               beta.x[5]*prevDMTefficacy== hiiumhigh (reference low efficacy)
#'               beta.x[6]*prevDMTefficacy== none (reference low efficacy)
#'               beta.x[7]*premedicalcost
#' @param percentiles percentiles (cumulative, monotone increasing) to define each responder subgroup; vector of floats of length 6
#' @return data - a data frame of \code{n} rows and 8 columns of simulated data; data.frame


generate_data <- function( n,
                          seed = NA,
                          beta = c(-0.2, -0.2, -0.2, -0.2, -0.2), # to create a homogeneous level of treatment effect
                          beta.x = c(-1.54, -0.01, 0.06, 0.25, 0.5, 0.13, 0.0000003), # to calculate the outcome
                          percentiles = seq(0, 1, by = 0.2)){


  set.seed(seed)

  if (percentiles[1] != 0 | percentiles[length(percentiles)] != 1 | length(percentiles) != 6) {
    stop("Wrong values of percentiles!")
  }

  # Create an empty shell
  ds <- data.table(NULL)

  # Define X, A, and time
  ds[, gender := rbinom(n = n, size = 1, prob = 0.75)]
  ds[, ageatindex_centered := round(rtruncnorm(n, a = 18, b = 64, mean = 48, sd = 12), 0) - 48 ] # rounded to integers
  ds[, prerelapseNum := rpois(n = n, lambda = 0.44)]
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

  XB <- model.matrix(~.,ds) %*% c(1.22,0.3,-0.1,0.2, 0.35,0.7,-0.00005,0.17,0.02,0)# ~75% people allocated in DMF arm based on (age,gender,prerelapseNum,DMT efficacy,costs,numSymptoms)
  pi <- exp(XB)/(1 + exp(XB))
  ds[, trt := as.numeric(runif(n) <=pi)]
  ds[, treatment := as.factor(ifelse(trt == 1, "DMF", "TERI"))]

  # Define Y (using all above PS predictors except for numSymptoms)
  xmat.score <- as.matrix(model.matrix(~ ageatindex_centered + gender + prerelapseNum + prevDMTefficacy + premedicalcost, ds))
  gamma <- matrix(c(-0.33, # Intercept
                    -0.001, # Age
                    0.05, # female
                    -0.002, # prerelapseNum
                    0.33, 0.02, # Medium/high and none DMT efficacy
                    -0.0000005), nrow = 7) # premedicalcost
  score = exp(xmat.score %*% gamma)
  ds[, Iscore := cut(score, quantile(score, percentiles), include.lowest = T, labels = seq(1, 5))]

  xmat.rate <- as.matrix(model.matrix(~ ageatindex_centered + gender + prerelapseNum +
                                        prevDMTefficacy + premedicalcost +
                                        Iscore + trt + trt*Iscore, ds))

  betas <- matrix(c(beta.x[1], # Intercept
                    beta.x[2], # Age
                    beta.x[3], # female
                    beta.x[4], # prerelapseNum
                    beta.x[5],
                    beta.x[6], # Medium/high and none DMT efficacy
                    beta.x[7], # premedicalcost
                    # Previous beta.x: -1.54 Intercept, -0.01 Age, 0.06 female, 0.47 prerelapseNum, 0.68, 0.13 Medium/high and none DMT efficacy, 0.000003 premedicalcost
                    0, 0, 0, 0, # Iscore categories 2:5 (reference group is Iscore1, high responders to DMF)
                    beta[1], beta[2] - beta[1], beta[3] - beta[1], beta[4] - beta[1], beta[5] - beta[1]))
  rate <- exp(xmat.rate %*% betas)

  ds[, y := rnegbin(n = n, mu = rate*finalpostdayscount/365.25, theta = 3)] # post treatment number of relapses

  ds[, years := finalpostdayscount / 365.25]
  ds[, age := ageatindex_centered + 48]
  ds[, gender:=factor(gender,labels=c("Male","Female"))]
  ds[, logPremedicalcost:=log(premedicalcost)]
  data <- ds[,c("age","gender", "prevDMTefficacy", "logPremedicalcost", "numSymptoms", "prerelapseNum", "treatment", "y","years","Iscore")]

  return(data)
}


# Missing generation ----
#' Function to generate missing data ----
#'
#' @param data Completed dataset
#' @param scenario Options of missing generation on logPremedicalcost covariate:
#'                  MCAR(Missing complete at random),
#'                  MAR(P(R)=f(age,gender,treatment)),
#'                  MART(P(R)=f(age,gender,treatment)),
#'                  MARTY(P(R)=f(age,gender,treatment,Y)),
#'                  MNAR(P(R)=f(age,gender,prerelapseNum) )
#' @param seed
#'
#' @return incomplete dataset

get_missdata <- function(data, scenario = "MAR", seed = 1234){

  data0<-data
  data <- dplyr::select(data, -Iscore)


  # Transform categorical data to dummies
  dat_misspattern <- data %>% mutate(gender = ifelse(gender == "Female", 1, 0)) %>%
    mutate(dummy=1) %>%
    mutate(treatment = ifelse(treatment == "DMF", 1, 0)) %>%
    mutate(dummy=1) %>%
    spread(key=prevDMTefficacy, value=dummy, fill=0, sep ="_") %>%
    mutate(dummy=1) %>%
    spread(key=numSymptoms, value=dummy, fill=0, sep ="_")


  md_temp <- ampute(dat_misspattern, mech = "MAR")

  # Generate missing data patterns for age and prerelapseNum (MAR)
  pattern <- data.frame(matrix(1,nrow = 3,  ncol(dat_misspattern)))
  colnames(pattern) <- colnames(md_temp$weights)
  pattern[c(1,3),"age"] <- 0
  pattern[c(2,3),"prerelapseNum"] <- 0
  weights <- data.frame(matrix(0,nrow = 3,  ncol(dat_misspattern)))
  colnames(weights) <- colnames(md_temp$weights)
  weights$gender <- 0.7
  weights$y <- 1
  weights$age[2] <- 1.25
  weights$prerelapseNum[1] <- 0.95
  freq<-c(0.3, 0.4, 0.2)

  md2 <- ampute(dat_misspattern, freq=freq, patterns = pattern, prop = 0.20, mech = "MAR")


  # Generate missing data patterns for prevDMTefficacy and numSymptoms (MAR)
  cols_prevDMTefficacy <- grepl( "prevDMTefficacy", colnames(dat_misspattern), fixed = TRUE)
  cols_numSymptoms <- grepl( "numSymptoms", colnames(dat_misspattern), fixed = TRUE)
  pattern <- data.frame(matrix(1,nrow = 3,  ncol(dat_misspattern)))
  colnames(pattern) <- colnames(md_temp$weights)
  pattern[c(1,3),which(cols_prevDMTefficacy)] <- 0
  pattern[c(2,3),which(cols_numSymptoms)] <- 0
  # Alter the weights such that missingness only depends on observed values of
  # * prevDMTefficacy (for pattern 2)
  # * numSymptoms (for pattern 1)
  # * age (for all patterns)
  # * gender (for all patterns)
  weights <- data.frame(matrix(0,nrow = 3,  ncol(dat_misspattern)))
  colnames(weights) <- colnames(md_temp$weights)
  weights$age <- -1
  weights$gender <- 1
  weights[2, cols_prevDMTefficacy] <- 1.2
  weights[1, cols_numSymptoms] <- 0.9
  freq<-c(0.35, 0.35, 0.3)
  md3 <- ampute(dat_misspattern, patterns = pattern, freq=freq,weights = weights, prop = 0.35, mech = "MAR")

  # Set missing data for logPremedicalcost
  pattern <- rep(1, ncol(dat_misspattern))
  pattern[which( colnames(dat_misspattern)=="logPremedicalcost")] <- 0

  weights <- rep(0, ncol(dat_misspattern))
  names(weights) <- colnames(md_temp$weights)

  if (scenario == "MCAR") {
    mech <- "MCAR"
  } else if (scenario == "MAR") {
    weights["age"] <- -1
    weights["gender"] <- 1
    mech <- "MAR"
  } else if (scenario == "MART") {
    weights["age"] <- -1
    weights["gender"] <- 1
    weights["treatment"] <- 1.5
    mech <- "MAR"
  } else if (scenario == "MARTY") {
    weights["age"] <- -1
    weights["gender"] <- 1
    weights["treatment"] <- 1.5
    weights["y"] <- 1.5
    mech <- "MAR"
  } else if (scenario == "MNAR") {
    weights["age"] <- -1
    weights["gender"] <- 1
    weights["logPremedicalcost"] <- 1.5
    mech <- "MNAR"
  } else {
    stop("Scenario not supported!")
  }
  md4 <- ampute(dat_misspattern, patterns = pattern, weights = weights, prop = 0.40, mech = mech)

  ampdata <- dat_misspattern
  ampdata$prerelapseNum <- md2$amp$prerelapseNum
  ampdata$age <- md2$amp$age
  ampdata[,which(cols_prevDMTefficacy)] <- md3$amp[,which(cols_prevDMTefficacy)]
  ampdata[,which(cols_numSymptoms)] <- md3$amp[,which(cols_numSymptoms)]
  ampdata$logPremedicalcost <- md4$amp$logPremedicalcost
  colnames(ampdata)

  ## Collapse dummies back into categorical variables
  ampdata$gender <-
  ampdata$gender <- factor(ampdata$gender, levels = c(1,0), labels = c("Female", "Male"))
  ampdata$treatment <- factor(ampdata$treatment, levels = c(1,0), labels = c("DMF", "TERI"))
  ampdata <- ampdata %>% mutate(prevDMTefficacy = ifelse(is.na(prevDMTefficacy_None),NA,
                                                         ifelse(prevDMTefficacy_None == 1, "None",
                                                         ifelse(prevDMTefficacy_Low_efficacy == 1, "Low_efficacy", "Medium_high_efficacy"))),
                                numSymptoms = ifelse(is.na(numSymptoms_0),NA,ifelse(numSymptoms_0 == 1, "0", ifelse(numSymptoms_1 == 1, "1", ">=2"))))

  ampdata$prevDMTefficacy <- factor(ampdata$prevDMTefficacy,
                                    levels = c("None", "Low_efficacy", "Medium_high_efficacy"),
                                    labels = c("None", "Low_efficacy", "Medium_high_efficacy"))
  ampdata$numSymptoms <- factor(ampdata$numSymptoms,
                                levels = c("0", "1", ">=2"),
                                labels = c("0", "1", ">=2"))

  # Get rid of dummies
  cols_prevDMTefficacy <- which(grepl( "prevDMTefficacy_", colnames(ampdata), fixed = TRUE))
  cols_numSymptoms <- which(grepl( "numSymptoms_", colnames(ampdata), fixed = TRUE))

  out <- ampdata %>% dplyr::select(-starts_with("prevDMTefficacy_")) %>%
         dplyr::select(-starts_with("numSymptoms_"))%>%
         mutate(Iscore=data0$Iscore)

  return(out)
}



ATE_estimation <- function(data, # mice object or dataset
                           estimand = "ATE", # Estimate the ATE or ATT
                           PSform = "Pred", # Propensity score model form: "Pred" : only specify the confounders,"Mind": specify missing indicator and deterministically imputed variable,"MIMind": specify missing indicator and stochastically imputed variable.
                           approach = NULL, # within, across
                           model = "Homogeneous", # model to estimate: ,"y ~DMF*prevDMTefficacy + offset(log(years))"
                           analysis ="Undefined",
                           variable= NA){ # name of the analysis

  main.model <- "y ~ treatment*(gender + age + logPremedicalcost + prerelapseNum + prevDMTefficacy + numSymptoms) + offset(log(years))"

    if (model == "Homogeneous") {
      main.model <<- "y ~ treatment + gender + age + logPremedicalcost + prerelapseNum + prevDMTefficacy + numSymptoms + offset(log(years))"
    }

    type_data <- "data.frame"

    # Propensity score estimation ----

    # Propensity score models

    if (PSform == "Pred") { # only specify the confounders
      PS.formula <- treatment ~ gender + age + logPremedicalcost + prerelapseNum + prevDMTefficacy
    } else { # Mind pr MIMind
      PS.formula <- treatment ~  gender + age.mind + age + lpmc.mind + logPremedicalcost + prn.mind + prerelapseNum + pde.mind + prevDMTefficacy
    }

    # Include additional variables required for PS score

    ## Transform mice object to data.frame
    if(!inherits(data, "data.frame")){
      type_data ="mice"
      data <- complete(data, 'long', include = TRUE)
    }

    # Add missing indicators for Mind and MIMind approaches
    if (PSform != "Pred"){
    data <- data %>%
      mutate(age.mind = as.numeric(is.na(age) == FALSE), #  missing indicator of age
             lpmc.mind = as.numeric(is.na(logPremedicalcost) == FALSE), # missing indicator of premedical cost
             prn.mind = as.numeric(is.na(prerelapseNum) == FALSE), # missing indicator of prerelapse number
             pde.mind = as.numeric(is.na(prevDMTefficacy) == FALSE), # missing indicator of previous DMT efficacy
             pns.mind = as.numeric(is.na(numSymptoms) == FALSE)) #missing indicator of number of symptoms

     if (PSform == "Mind"){
       data <- data %>%
         mutate(age = ifelse(is.na(age),0,age), # deterministic imputation of age
                logPremedicalcost = ifelse(is.na(logPremedicalcost),0,logPremedicalcost), # # deterministic imputation of lo premedical cost
                prerelapseNum = ifelse(is.na(prerelapseNum),0,prerelapseNum), #  deterministic imputation of prerelapse number
                prevDMTefficacy = as.factor(ifelse(is.na(prevDMTefficacy),"na",prevDMTefficacy)), # deterministic imputation of previous DMT efficacy
                numSymptoms = as.factor(ifelse(is.na(numSymptoms),"na",numSymptoms))) # deterministic imputation of number of symptom
        }
    }

    ## Transform back to mice object
    if(type_data =="mice"){
      type_data ="mice"
      data <- as.mids(data)
    }

    # Get balanced data ----

    ## Set parameters according to estimator ----
    if (estimand == "ATE"){ # for ATE
      methodv <- "full"
      replacev <- FALSE
    } else {             # for ATT
      methodv <- "nearest"
      replacev <- TRUE
    }

    ATE_var <- NA
    if (type_data == "data.frame"){
        ### data.frame object ----
        #### data matching ----

        mout <- matchit(PS.formula,
                        data = data,
                        family = binomial,
                        method = methodv,
                        caliper = 0.2,
                        std.caliper = TRUE,
                        estimand = estimand,
                        distance = "glm",
                        link = "logit",
                        replace = replacev) # we apply with replacement here due to small number ofDMFed patients


          mdata <- match.data(mout)
          assign("mdata", mdata, envir = .GlobalEnv)

      #### glm model ----

          fit_mod <- glm(formula = as.formula(main.model),
                         family = poisson(link = "log"),
                         data = mdata,
                         weights = weights)


            ATE <- avg_comparisons(fit_mod,
                                   variables = list(treatment = c("TERI","DMF")),
                                   vcov = "HC3", # Var-cov correction
                                   wts = "weights",
                                   newdata = if (estimand == "ATT"){
                                     subset(mdata, treatment == "DMF")}else{mdata}
                                   )%>%data.frame()
           if(!is.na(variable)){
            ATE_var <-  comparisons( fit_mod,
                                 vcov = "HC3",
                                 variables = list(treatment = c("TERI","DMF")),
                                 wts = "weights",
                                 by = variable,
                                 newdata = if (estimand == "ATT"){
                                   subset(mdata, treatment == "DMF")}else{mdata}
                                  )%>%data.frame()
            }

      }else{ # Imputed datasets

        #### balancing  ----
        mdata <- matchthem(formula = PS.formula,
                                      datasets = data,
                                      approach = approach,
                                      method = methodv,
                                      caliper = 0.1,
                                      family = binomial,
                                      estimand =estimand,
                                      distance = "glm",
                                      link = "logit",
                                      replace = replacev)


        #### glm model ----
        dat <- complete( mdata, action = "all")

        mod <- lapply(dat, \(i)
                      glm(formula = as.formula(main.model),
                          family = poisson(link = "log"),
                          weights = weights,
                          data = i))
        # ATE estimation

        ATE_i <-lapply(seq_along(mod), \(i)
                     avg_comparisons(mod[[i]],
                                     variables = list(treatment = c("TERI","DMF")),
                                     vcov = "HC3",
                                     wts = "weights",
                                     newdata = if (estimand == "ATT"){
                                                subset(dat[[i]], treatment == "DMF")}else{dat[[i]]}))
        ATE<- data.frame(summary(mice::pool(ATE_i), conf.int = TRUE))%>%
              rename('conf.low'="X2.5..",'conf.high'="X97.5..")

        if(!is.na(variable)){

          ATE_var_i <-lapply(seq_along(mod), \(i)
                         comparisons( mod[[i]],
                                   vcov = "HC3",
                                   variables = list(treatment = c("TERI","DMF")),
                                   wts = "weights",
                                   by = variable,
                                   newdata = if (estimand == "ATT"){
                                     subset(dat[[i]], treatment == "DMF")}else{dat[[i]]}))


          vec_val <- as.vector(unlist(unique(setDT(dat[[1]])[,..variable])))

          ATE_var <- data.frame()
          for( val in vec_val ){
            res_val <- data.frame(summary(mice::pool(lapply(seq_along( ATE_var_i ), \(i)
                                                      ATE_var_i[[i]]%>%
                                                        filter( get({{variable}}) == val))),
                                          conf.int = TRUE))%>%
              mutate(!!variable:= val)%>%
              rename('conf.low'="X2.5..",'conf.high'="X97.5..")

            ATE_var <- rbind(ATE_var,res_val)
          }
        }
      }

    ATE$analysis <- analysis
    if(!is.na(variable)){ATE_var$analysis <- analysis}

    return(list( ATE = ATE, ATE_var = ATE_var ))}



















# Function to estimate the treatment effect ----

get_est <- function(data, # mice object or dataset
                    estimand = "ATE", # Estimate the ATE or ATT
                    PSform = "Pred", # Propensity score model form: "Pred" : only specify the confounders,"Mind": specify missing indicator and deterministically imputed variable,"MIMind": specify missing indicator and stochastically imputed variable.
                    method = "Matching", #Balance method:"Matching": Propensity score matching,"IPTW": Inverse propensity treatment weighting)
                    CC = FALSE, # use the complete case dataset,
                    approach = NULL, # within, across
                    model = "Homogeneous", # model to estimate: ,"y ~DMF*prevDMTefficacy + offset(log(years))"
                    analysis ="Undefined"){ # name of the analysis


  if( model == "Homogeneous"){
    main.model <<- "y ~ DMF + offset(log(years))"
    }else{
      main.model <<- "y ~ DMF+ DMF*prevDMTefficacy + offset(log(years))"
      }

  type_data ="data.frame"

  # Propensity score estimation ----

  # Propensity score models

  if (PSform == "Pred") { # only specify the confounders
    PS.formula <- treatment ~ gender + age + logPremedicalcost + prerelapseNum + prevDMTefficacy
  } else if (PSform == "Mind") { # Mind
    PS.formula <- treatment ~ gender + age.mind + age.dimp + lpmc.mind + lpmc.dimp + prn.mind + prn.dimp + pde.mind + pde.dimp
  } else if (PSform == "MIMind") {  # "MIMind"
    PS.formula <- treatment ~  gender + age.mind + age + lpmc.mind + logPremedicalcost + prn.mind + prerelapseNum + pde.mind + prevDMTefficacy
  } else if (PSform == "MIMind") {  # "MIMind"
    PS.formula <- treatment ~  gender + age.mind + age + lpmc.mind + logPremedicalcost + prn.mind + prerelapseNum + pde.mind + prevDMTefficacy
  }

  # Include additional variables required for PS score

   ## Transform mice object to data.frame
  if(!inherits(data, "data.frame")){
    type_data ="mice"
    data <- complete(data, 'long', include = TRUE)
  }

  data <- data %>%
            mutate(age.mind = as.numeric(is.na(age) == FALSE), #  missing indicator of age
                   age.dimp = ifelse(is.na(age),0,age), # deterministic imputation of age
                   lpmc.mind = as.numeric(is.na(logPremedicalcost) == FALSE), # missing indicator of premedical cost
                   lpmc.dimp = ifelse(is.na(logPremedicalcost),0,logPremedicalcost), # # deterministic imputation of lo premedical cost
                   prn.mind = as.numeric(is.na(prerelapseNum) == FALSE), # missing indicator of prerelapse number
                   prn.dimp = ifelse(is.na(prerelapseNum),0,prerelapseNum), #  deterministic imputation of prerelapse number
                   pde.mind = as.numeric(is.na(prevDMTefficacy) == FALSE), # missing indicator of previous DMT efficacy
                   pde.dimp = as.factor(ifelse(is.na(prevDMTefficacy),"na",prevDMTefficacy)), # deterministic imputation of previous DMT efficacy
                   pns.mind = as.numeric(is.na(numSymptoms) == FALSE), #missing indicator of number of symptoms
                   pns.dimp = as.factor(ifelse(is.na(numSymptoms),"na",numSymptoms))) # deterministic imputation of number of symptom

   ## Transform back to mice object
   if(type_data =="mice"){
      type_data ="mice"
      data <- as.mids(data)
      }

   ## Get Complete Case dataset
   if(CC){
     data <- data %>% filter(complete.cases(.))
     }


  # Get balanced data ----

  ## Set parameters according to estimator ----
  if (estimand == "ATE"){ # for ATE
    methodv <- "full"
    replacev <- FALSE
  } else {             # for ATT
    methodv <- "nearest"
    replacev <- TRUE
  }


  ## PS Matching ----
  if(method=="Matching"){

    if (type_data == "data.frame"){
      ### data.frame object ----
      #### data matching ----

      mout <- matchit(PS.formula,
                      data = data,
                      family = binomial,
                      method = methodv,
                      caliper = 0.2,
                      std.caliper = TRUE,
                      estimand = estimand,
                      distance = "glm",
                      link = "logit",
                      replace = replacev) # we apply with replacement here due to small number ofDMFed patients

      #### glm model ----

      if (estimand == "ATT") {
        mdata <- as.data.table(get_matches(object = mout, id = "matching_id"))
        assign("mdata", mdata, envir = .GlobalEnv)
        match_mod <- glm(formula = as.formula(main.model),
                         family = poisson(link = "log"),
                         data = mdata)
        # Estimate cluster-robust standard error
        tx_var <- vcovCL(match_mod, cluster = ~ subclass + matching_id, sandwich = TRUE)

        TE <-avg_predictions(fit1, variables = "treatment",
                        vcov = "HC3",
                        newdata = subset(md, A == 1),
                        wts = "weights")

      } else if (estimand == "ATE") {
        mdata <- as.data.table(MatchIt::match.data(object = mout))
        assign("mdata", mdata, envir = .GlobalEnv)
        match_mod <- glm(formula =as.formula(main.model),
                         family = poisson(link = "log"),
                         data = mdata,
                         weights = weights)
        # Estimate robust variance-covariance matrix
        tx_var <- vcovCL(match_mod, cluster = ~ subclass, sandwich = TRUE)

      } else {
        stop ("Estimand not supported!")
      }

    # Result object
    estimate <- coefficients(match_mod)
    std.error <- sqrt(diag(tx_var))
    result <- data.frame(analysis = analysis, method = method, estimand=estimand, term=names(estimate), estimate=estimate,std.error=std.error)
    result$LCI <- result$estimate+qnorm(0.025)*std.error
    result$UCI <- result$estimate+qnorm(0.975)*std.error



      ### MICE object ----
    }else{

      ####data matching ----
      matched.datasets <- matchthem(formula = PS.formula,
                                    datasets = data,
                                    approach = approach,
                                    method = methodv,
                                    caliper = 0.1,
                                    family = binomial,
                                    estimand =estimand,
                                    distance = "glm",
                                    link = "logit",
                                    replace = replacev)
      #### glm model ----

      matched.results<- summary(mice::pool(with(matched.datasets,
                                           svyglm(formula = as.formula(main.model),
                                                  family = poisson(link = "log")),
                                           cluster = TRUE)),
                                 conf.int = TRUE)


      # Result object
      result <- data.frame(analysis = analysis, method = method, estimand=estimand,
                           term = matched.results$term,
                           estimate=matched.results$estimate,
                           std.error=matched.results$std.error,
                           LCI= matched.results$`2.5 %`,
                           UCI= matched.results$`97.5 %`)

    }

  }

  # IPWT ----
  if(method=="IPTW"){
    if (inherits(data, "data.frame")){
      ### data.frame object ----

      #### data weghting ----
      wout  <- weightit(PS.formula, data = data, estimand = estimand, method = "ps")
      rhcSvy  <- survey::svydesign(ids = ~ 1, data = data, weights = ~ wout$weights)

      #### glm model ----
      ipw_mod <- survey::svyglm(formula = as.formula(main.model),
                                family  =  poisson(link = "log"),
                                design = rhcSvy)

      # Note: svyglm always returns 'model-robust' standard errors;
      # the Horvitz-Thompson-type standard errors used everywhere in the survey
      # package are a generalisation of the model-robust 'sandwich' estimators.
      estimate <- coef(ipw_mod)
      std.error <- sqrt(diag(vcov(ipw_mod)))
      result <- data.frame(analysis = analysis, method = method, estimand=estimand, term=names(estimate), estimate=estimate,std.error=std.error)
      result$LCI <- result$estimate+qnorm(0.025)*std.error
      result$UCI <- result$estimate+qnorm(0.975)*std.error

    }else{
      ### MICE object ----

      #### data weghting ----
      weighted.datasets <- weightthem(PS.formula,
                                      datasets = data,
                                      approach = approach,
                                      method = 'ps',
                                      estimand = estimand)
      #### glm model ----
      weighted.results <- summary(mice::pool(with(weighted.datasets,
                                            svyglm(formula =as.formula(main.model), family = poisson(link = "log")),
                                            cluster = TRUE)),conf.int=TRUE)

      result <- data.frame(analysis = analysis, method = method, estimand=estimand,
                           term = weighted.results$term,
                           estimate=weighted.results$estimate,
                           std.error=weighted.results$std.error,
                           LCI= weighted.results$`2.5 %`,
                           UCI= weighted.results$`97.5 %`)


    }
  }

  # Calculate Treatment effect
  if( model == "Homogeneous"){
    TE <- data.frame(analysis = analysis,
                     method = method,
                     estimand=estimand,
                     term="Treatment",
                     exp(c(1,1)%*%as.matrix(result[,c("estimate","LCI","UCI")]))
                    -exp(c(1,0)%*%as.matrix(result[,c("estimate","LCI","UCI")])))
  }else{
    TE_none <- exp(c(1,1,0,0,0,0)%*%as.matrix(result[,c("estimate","LCI","UCI")])) - exp(c(1,0,0,0,0,0)%*%as.matrix(result[,c("estimate","LCI","UCI")]))
    TE_low <- exp(c(1,1,1,0,1,0)%*%as.matrix(result[,c("estimate","LCI","UCI")])) - exp(c(1,0,1,0,0,0)%*%as.matrix(result[,c("estimate","LCI","UCI")]))
    TE_mhig <- exp(c(1,1,0,1,0,1)%*%as.matrix(result[,c("estimate","LCI","UCI")])) - exp(c(1,1,0,1,0,0)%*%as.matrix(result[,c("estimate","LCI","UCI")]))
    TE <- data.frame(analysis = analysis,
                     method = method,
                     estimand=estimand,
                     term=c("None","Low","Medium_High"),
                     rbind(TE_none,TE_low,TE_mhig))}

  return(TE)
}


# No output in the prediction model
form_nt <- list(prevDMTefficacy ~ age + gender + logPremedicalcost + numSymptoms + prerelapseNum + logyears,
                logPremedicalcost ~ age + gender + prevDMTefficacy + numSymptoms + prerelapseNum + logyears,
                numSymptoms ~ age + gender + logPremedicalcost + prevDMTefficacy + prerelapseNum + logyears,
                prerelapseNum ~ age + gender + logPremedicalcost + prevDMTefficacy + numSymptoms + logyears,
                age ~ prerelapseNum + gender + logPremedicalcost + prevDMTefficacy + numSymptoms + logyears)
form_nt <- name.formulas(form_nt)

# Output in the prediction model
form_y <- list(prevDMTefficacy ~ age + gender + logPremedicalcost + numSymptoms + treatment + prerelapseNum + logyears + y,
               logPremedicalcost ~ age + gender + prevDMTefficacy + numSymptoms + treatment + prerelapseNum + logyears + y,
               numSymptoms ~ age + gender + logPremedicalcost + prevDMTefficacy + prerelapseNum + treatment + logyears + y,
               prerelapseNum ~ age + gender + logPremedicalcost + prevDMTefficacy + numSymptoms + treatment + logyears + y,
               age ~ prerelapseNum + gender + logPremedicalcost + prevDMTefficacy + numSymptoms + treatment + logyears + y)
form_y <- name.formulas(form_y)

# Prediction model with interaction T*X
form_i <- list(prevDMTefficacy ~ treatment*(age + gender + logPremedicalcost + numSymptoms + prerelapseNum + logyears),
               logPremedicalcost ~ treatment*(age + gender + prevDMTefficacy + numSymptoms + prerelapseNum + logyears),
               numSymptoms ~ treatment*(age + gender + logPremedicalcost + prevDMTefficacy + prerelapseNum + logyears),
               prerelapseNum ~ treatment*(age + gender + logPremedicalcost + prevDMTefficacy + numSymptoms + logyears),
               age ~ treatment*(prerelapseNum + gender + logPremedicalcost + prevDMTefficacy + numSymptoms + logyears))
form_i <- name.formulas(form_i)

# Prediction model with interaction T*X, T*Y with X*Y!!
form_iy <- list(prevDMTefficacy ~ (y + treatment)*(age + gender + logPremedicalcost + numSymptoms + prerelapseNum + logyears) + y*treatment,
                logPremedicalcost ~ (y + treatment)*(age + gender + prevDMTefficacy + numSymptoms + prerelapseNum + logyears) + y*treatment,
                numSymptoms ~ (y + treatment)*(age + gender + logPremedicalcost + prevDMTefficacy + prerelapseNum + logyears) + y*treatment,
                prerelapseNum ~ (y + treatment)*(age + gender + logPremedicalcost + prevDMTefficacy + numSymptoms + logyears) + y*treatment,
                age ~ (y + treatment)*(prerelapseNum + gender + logPremedicalcost + prevDMTefficacy + numSymptoms + logyears) + y*treatment)
form_iy <- name.formulas(form_iy)

# Prediction model with interaction T*X, T*Y
form_iytt <- list(prevDMTefficacy ~ treatment*(age + gender + logPremedicalcost + numSymptoms + prerelapseNum + logyears + y),
                 logPremedicalcost ~ treatment*(age + gender + prevDMTefficacy + numSymptoms + prerelapseNum + logyears + y),
                 numSymptoms ~ treatment*(age + gender + logPremedicalcost + prevDMTefficacy + prerelapseNum + logyears + y),
                 prerelapseNum ~ treatment*(age + gender + logPremedicalcost + prevDMTefficacy + numSymptoms + logyears + y),
                 age ~ treatment*(prerelapseNum + gender + logPremedicalcost + prevDMTefficacy + numSymptoms + logyears + y))
form_iytt <- name.formulas(form_iytt)

form_yt  <- list(prevDMTefficacy ~ age + gender + logPremedicalcost + numSymptoms + treatment + prerelapseNum + logyears + y,
               logPremedicalcost ~ age + gender + prevDMTefficacy + numSymptoms + treatment + prevDMTefficacy*treatment + prerelapseNum + logyears + y,
               numSymptoms ~ age + gender + logPremedicalcost + prevDMTefficacy + prerelapseNum + treatment + prevDMTefficacy*treatment + logyears + y,
               prerelapseNum ~ age + gender + logPremedicalcost + prevDMTefficacy + numSymptoms + treatment + prevDMTefficacy*treatment + logyears + y,
               age ~ prerelapseNum + gender + logPremedicalcost + prevDMTefficacy + numSymptoms + treatment + prevDMTefficacy*treatment + logyears + y)
form_yt <- name.formulas(form_yt)






library(ggplot2)

#' Plot distribution full dataset
#' @param data full dataset


# Density plot by treatment
plot_den <- function(datap,var,varlab) {
  mu<-datap[,.(grp.mean=mean(get(var))),treatment]
  ggplot(datap, aes(x=get(var), color=treatment)) +
    geom_density()+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=treatment),
               linetype="dashed")+
    labs(x=varlab, y = "Density",colour="Treatment")+
    scale_color_brewer(palette="Accent") + theme_minimal()+theme(legend.position="top")
}

# Proportion ment

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

# Annual relapse rate plot
Relapserate_plot <-function(data){
  data[,AR:=y/years]
  data[,prevDMTefficacy:=as.factor(prevDMTefficacy)]
  levels(data$prevDMTefficacy)<-c("None","Low_efficacy","Medium_high_efficacy")
  datallgroup<-data[,.(grp.mean=mean(AR),sd=sd(AR),n=.N),list(treatment,prevDMTefficacy)]
  datallgroup[,LCI:=grp.mean-1.96*sd/sqrt(n)]
  datallgroup[,UCI:=grp.mean+1.96*sd/sqrt(n)]
  pd <- position_dodge(.1)
  ggplot(datallgroup,aes(x=prevDMTefficacy,y=grp.mean,group=treatment,colour=treatment))+
    geom_errorbar(aes(ymin=LCI, ymax=UCI),
                  width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size=2)+
    labs(y="Annual post relapse rate",x="Previous DMT treatment efficacy ",color="")+
    scale_colour_manual(values = c("#BEAED4","#7FC97F"))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45),legend.position = "top")
}





Relapserate_plot0 <-function(data){
  data[,AR:=y/years]
  data[,Iscore:=as.factor(Iscore)]
  levels(data$Iscore)<-c("High DMF","Moderate DMF","Neutral","Moderate TERI","High TERI")
  datallgroup<-data[,.(grp.mean=mean(AR),sd=sd(AR),n=.N),list(treatment,Iscore)]
  datallgroup[,LCI:=grp.mean-1.96*sd/sqrt(n)]
  datallgroup[,UCI:=grp.mean+1.96*sd/sqrt(n)]
  pd <- position_dodge(.1)
  ggplot(datallgroup,aes(x=Iscore,y=grp.mean,group=treatment,colour=treatment))+
    geom_errorbar(aes(ymin=LCI, ymax=UCI),
                  width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size=2)+
    labs(y="Annual post relapse rate",x="Effect treatment group",color="")+
    scale_colour_manual(values = c("#BEAED4","#7FC97F"))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45),legend.position = "top")
}


plot_dis <-function(data){
  p1<-plot_den(datap=data,var="age",varlab="Age (years)")  #Age density
  p2<-plot_count(datap=data,var="gender",varlab="Gender") # Gender proportion
  levels(data$prevDMTefficacy)<-c("None","Low","Med-High")
  p3<-plot_count(datap=data,var="prevDMTefficacy",varlab="Previous DMT efficacy") # Previuos treatment proportion
  p4<-plot_den(datap=data,var="logPremedicalcost",varlab="Log previous medical cost")  #log cost density
  p5<-plot_count(datap=data,var="prerelapseNum",varlab="Previous relapse number") # Previous relapses proportion
  p6<-Relapserate_plot0(data=data)
  fplot<-ggpubr::ggarrange(p2,p1,p3,p4,p5,p6, ncol = 3, nrow = 2,common.legend = TRUE, legend="bottom")# All plots together
  return(fplot)
  }



formatMSdata <- function(data) {
  label(data$age)       <- "Age"
  label(data$gender)       <- "Gender"
  label(data$logPremedicalcost)       <- "Log prior medical costs"
  label(data$prevDMTefficacy)       <- "Efficacy of previous DMT"
  label(data$numSymptoms)       <- "Number of prior symptoms"
  label(data$prerelapseNum)       <- "Number of prior relapses"

  units(data$age)       <- "years"

  return(data)
}

