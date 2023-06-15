# ------------------------------------------------------------------
# Project: Chapter 18 of the RWE book
#
# Purpose: To provide supporting functions for chapter16.R
#
# Platform: Windows
# R Version: 4.2.2
#
#   Modifications:
#
#   Date			By			Description
# --------		--------	-----------------------------
#  16FEB2022  pj      Start the script
#  04MAR2022  pj      Add cvmodel functions and getValue
#  08MAR2022  pj      Modify y for list and dwols, modify true d calculation, add accuracy and agreement plots
#  13JUN2022  pj      Change from A to T for treatment
#  17AUG2022  gs      Added comments + minor edits
#  26AUG2022  pj      Go through GS' comments and make edits
#  31AUG2022  pj      Suppress startup message when loading library
# ------------------------------------------------------------------


simcountdata <- function(n, seed = 999,
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
  #'               beta.x[5]*prevDMTefficacy== mediumhigh (reference low efficacy)
  #'               beta.x[6]*prevDMTefficacy== none (reference low efficacy)
  #'               beta.x[7]*premedicalcost
  #' @param percentiles percentiles (cumulative, monotone increasing) to define each responder subgroup; vector of floats of length 6
  #' @return data - a data frame of \code{n} rows and 8 columns of simulated data; data.frame

  library(truncnorm)
  library(magrittr)
  suppressPackageStartupMessages(library(reshape2))
  library(fastDummies)

  set.seed(seed)

  if (percentiles[1] != 0 | percentiles[length(percentiles)] != 1 | length(percentiles) != 6){message("Wrong values of percentiles!")}

  # Create an empty shell
  ds <- data.frame(matrix(NA, nrow = n, ncol = 10))
  colnames(ds) <- c("trt", "ageatindex_centered", "female", "prerelapse_num",
                    "prevDMTefficacy", "premedicalcost", "numSymptoms",
                    "postrelapse_num", "finalpostdayscount", "group")

  # Define X, A, and time
  ds %<>%
    mutate(trt =                 rbinom(n = n, size = 1, prob = 0.75),
           # treatment =           ifelse(trt == 1, "drug1", "drug0"),
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
                      Iscore = cut(score, quantile(score, percentiles), include.lowest = TRUE, labels = seq(1, 5)))

  xmat.rate <- model.matrix(~ ageatindex_centered + female + prerelapse_num + prevDMTefficacy + premedicalcost +
                              Iscore + trt + trt*Iscore, ds) %>% as.matrix()

  betas <- matrix(c(beta.x[1], # Intercept
                    beta.x[2], # Age
                    beta.x[3], # female
                    beta.x[4], # prerelapse_num
                    beta.x[5], beta.x[6], # Medium/high and none DMT efficacy
                    beta.x[7], # premedicalcost
                    0, 0, 0, 0, # Iscore categories 2:5 (reference group is Iscore1, high responders to DMF)
                    beta[1], beta[2] - beta[1], beta[3] - beta[1], beta[4] - beta[1], beta[5] - beta[1]))
  rate <-  exp(xmat.rate %*% betas)
  ds <- ds %>% mutate(postrelapse_num = rnegbin(n = n, mu = rate*finalpostdayscount/365.25, theta = 3))

  return(list(data = ds, betas = betas, percentiles = percentiles))
}

IPTWfun <- function(PSmodel, data, newdata = NULL){
  #' Estimate PS and append PS and IPTW to dataset
  #'
  #' @param PSmodel A formula with treatment variable of the left-hand side and covariates separated by "+" on the right-hand side
  #' @param data Dataframe from which to fetch all variables in PSmodel
  #' @param newdata Data for which we want to predict the PS. If NULL, dataframe supplied in \code{data} is used
  #'
  #' @return Same dataframe as supplied by \code{newdata}, with 2 additional columns ps and itpw

  if(class(PSmodel) != "formula") stop("PSmodel must be a formula.")
  if(is.null(newdata)) newdata <- data

  ps <- glm(PSmodel, family = "binomial", data = data)
  newdata$ps <- predict(ps, newdata, type = "response")
  newdata <- newdata %>% mutate(iptw = ifelse(trt == 1, 1/ps, 1/(1 - ps)))
  # Note: if one category absent in training data, PS coefficient will be NA, need to make sure to keep only the columns in testdata for which coefficient is not NA
  return(newdata)
}

cvvalue <- function(data, xvar, method, n.fold, n.cv, seed, RCT = FALSE){
  #' n.fold cross-validation to estimate the value under the estimated optional ITR
  #'
  #' @param data Input data; data.frame
  #' @param method PM method to be implemented; string
  #'               Possible methods: 'all1','all0','dWOLS','listDTR2', 'poisson', 'tworeg', 'contrastReg'
  #' @param xvar X variables to be included in the PS model and PM model, must be numeric or dummy variables; vector
  #' @param n.fold Number of CV folds; integer
  #' @param n.cv Number of CV iterations; integer
  #' @param RCT Whether treatment is randomized (TRUE) or not. If RCT=TRUE, the PS is taken as the proportion of patients treated with 1
  #' @param seed randomization seed
  #'
  #' @return cvresult - A list of length n.cv, each with sublist of length n.fold, containing the following elements:
  #'              - dhat: the estimated optimal ITR d; a vector of length equal to the number of rows in this fold
  #'              - vhat.dhat: a list of two elements to calculate the value
  #'                    - U: the numerator of the value
  #'                    - W: the denominator of the value

  # Create empty shell for output
  cvresult <- vector("list", n.cv)
  names(cvresult) <- paste0("cv.i", 1:n.cv)

  for (cv.i in 1:n.cv){
    cat("\nCV iteration =", cv.i, "out of", n.cv)
    this.seed = seed*100 + cv.i
    set.seed(this.seed)

    # Determine outcome, which depends on the method
    candidates <- list(all1 = "y",
                       all0 = "y",
                       poisson = "y",
                       dWOLS = "mlogarr0001",
                       listDTR2 = "mlogarr0001",
                       contrastReg = "y",
                       twoReg = "y")
    yvar <- candidates[[method]]
    cat("\n  Outcome is:", yvar, "for method", method, "\n")

    # Format data
    input <- data.frame(y = data[[yvar]], trt = data$drug1, time = log(data$years), data[xvar])

    # Create CV folds
    folds <- createFolds(input$trt, k = n.fold, list = TRUE) # stratified CV

    # Go through each CV
    cvresult[[paste0("cv.i", cv.i)]] <- eachCV(data = input, method = method, xvar = xvar, outcome = "y", folds = folds, n.fold = n.fold, seed = seed, RCT = RCT)
  }# end of loop for cv.i

  return(cvresult)
}

eachCV <- function(data, method, xvar, outcome, folds, n.fold, seed, RCT){
  #' One iteration of cross-validation
  #' Estimate PS in each training and testing dataset after each splitting
  #'
  #' @param data Input data; data.frame
  #' @param method PM method to be implemented; string
  #'               Possible methods: 'all1','all0','dWOLS','listDTR2', 'poisson', 'twoReg', 'contrastReg'
  #' @param xvar X variables to be included in the PS model and PM model, must be numeric or dummy variables; vector
  #' @param outcome Response variable (whether higher outcome is better depends on method) specified by user; string
  #' @param folds Row indices of the input data in the same CV fold are saved in each element of list; list of length n.fold
  #' @param n.fold Number of CV folds; integer
  #' @param RCT Whether treatment is randomized (TRUE) or not. If RCT=TRUE, the PS is taken as the proportion of patients treated with 1
  #' @param seed randomization seed
  #'
  #' @return eachcv - Estimated value for each CV fold; vector of size n.fold

  # Create empty shell for output
  eachcv <- vector("list", n.fold)
  names(eachcv) <- paste0("fold", 1:n.fold)

  sim.big = NULL

  # Go through each CV fold
  for (fold.i in 1:n.fold){
    cat("\n CV fold =", fold.i, "out of", n.fold)
    traindata <- data[-folds[[fold.i]],]
    testdata <- data[folds[[fold.i]],]

    # Calculate PS/IPTW
    if (RCT){
      trainps <- mean(traindata$trt)
      traindata <- traindata %>% mutate(ps = trainps, iptw = ifelse(trt == 1, 1/ps, 1/(1-ps)))
      testdata <- testdata %>% mutate(ps = trainps, iptw = ifelse(trt = 1, 1/ps, 1/(1-ps)))
    } else {
      fpstr <- as.formula(paste("trt ~", paste0(xvar, collapse = "+")))
      traindata <- IPTWfun(data = traindata, PSmodel = fpstr)
      testdata <- IPTWfun(data = traindata, PSmodel = fpstr, newdata = testdata)
    }

    # Apply PM method
    if (method %in% c("all1", "all0")){
      if (method == "all1") d <- 1
      if (method == 'all0') d <- 0
      eachcv[[paste0("fold", fold.i)]]$dhat <- rep(d, nrow(testdata))
      eachcv[[paste0("fold", fold.i)]]$vhat.dhat <- getValue(y = testdata[["y"]],
                                                             t = testdata[["trt"]],
                                                             d.hat = rep(d, nrow(testdata)),
                                                             p.hat = testdata$ps,
                                                             fu = exp(testdata[['time']]))
    } else if (method == "poisson"){
      traindata1 <- traindata %>% filter(trt == 1) # T = 1
      traindata0 <- traindata %>% filter(trt == 0) # T = 0

      eachcv[[paste0("fold", fold.i)]] <- tryCatch({itrLuPoisson(traindata1 = traindata1,
                                                                 traindata0 = traindata0,
                                                                 testdata = testdata,
                                                                 xvar = xvar,
                                                                 sim.big = sim.big)
      }, error = function(error_message) {
        cat("\n")
        message(error_message)
        cat("\nSolution: Set U and W to NA and move on.\n")
        list(dhat = NA, vhat.dhat = list(U = NA, W = NA), dhat.big = NA, v.dhat= NA)
      }
      )
    } else if (method == "dWOLS"){
      eachcv[[paste0("fold", fold.i)]] <- itrDWOLS(traindata = traindata,
                                                   testdata = testdata,
                                                   outcome = outcome,
                                                   Xoutcome = c(xvar, "time"),
                                                   Xinteraction = xvar,
                                                   dWOLSweight = "IPTW",
                                                   sim.big = sim.big)
    } else if (method == "listDTR2"){
      eachcv[[paste0("fold", fold.i)]] <- tryCatch({ itrLIST(traindata = traindata,
                                                             testdata = testdata,
                                                             outcome = outcome, # listdtr assumes higher better
                                                             treatment = "trt",
                                                             xvar = xvar,
                                                             maxlen = 2L,
                                                             seed = seed + fold.i,
                                                             sim.big = sim.big)
      }, error = function(error_message) {
        cat("\n")
        message(error_message)
        cat("\nSolution: Set U and W to NA and move on.\n")
        list(dhat = NA, vhat.dhat = list(U = NA, W = NA), dhat.big = NA, v.dhat = NA)
      }
      )
    } else if (method %in% c("twoReg", "contrastReg")){
      eachcv[[paste0("fold", fold.i)]] <- tryCatch({itrLuDR(traindata = traindata,
                                                            testdata = testdata,
                                                            xvar = xvar,
                                                            RCT = RCT,
                                                            tree.depth = 2,
                                                            n.trees = 100,
                                                            Kfold = 5,
                                                            B = 3,
                                                            seed.cf = 3,
                                                            plot.gbmperf = FALSE,
                                                            sim.big = sim.big)$valueContrastReg
      }, error = function(error_message) {
        cat("\n")
        message(error_message)
        cat("\nSolution: Set U and W to NA and move on.\n")
        list(dhat = NA, vhat.dhat = list(U = NA, W = NA), dhat.big = NA, v.dhat = NA)
      }
      )
    } else {
      stop("Method not recognized! Pick one from ['all1', 'all0', 'dWOLS', 'listDTR2', 'poisson', 'twoReg', 'contrastReg'].")
    }# end of if-else statements for method
  }# end of CV fold loop for fold.i

  return(eachcv)
}

itrDWOLS <- function(traindata, testdata, outcome, Xoutcome, Xinteraction, dWOLSweight = "IPTW", sim.big = NULL){
  #' ITR based on doubly robust dWOLS
  #'
  #' @param traindata - training data with both arms; data.frame
  #' @param testdata - test data with both arms (assume two arms) and same variables as traindata; data.frame
  #' @param outcome - assumed to be favorable if lower, continuous for dWOLS; string
  #' @param Xoutcome - covariates to be included as main effects in the outcome model; a vector of strings
  #' @param Xinteraction - covariates with an interaction term with treatment, must be a subset of Xoutcome; a vector of strings
  #' @param weight - weights used in dWOLS; "IPTW" or "overlap"
  #' @param sim.big a large test set object (independent from data) to get true value function; list object with elements `data` and `betas`
  #'
  #' @return if testdata != NULL: a list with:
  #'              - dhat: the estimated ITR of every observation in `testdata` (1 means arm 1 and 0 means arm 0)
  #'              - vhat.dhat: estimated value of the estimated ITR of every observation in `testdata`
  #'     OR if testdata == NULL: a list with:
  #'              - fit: dataframe with coefficient of ITR with dWOLS (coef) and SE based on 200 bootstrap
  #'              - score: dataframe with ID, score (score_dwols), optimal ITR (itr_dwols)
  #' @return In addition, if sim.big != NULL: the list also contains:
  #'              - dhat.big: the estimated ITR of every observation in `sim.big` (1 means arm 1 and 0 means arm 0)
  #'              - v.dhat: true value of the estimated ITR of every observation in `sim.big`

  output <- list()

  # Save outcome input
  text.outcome <- outcome

  # Specify subset variables
  X.beta <- list(as.formula(paste0("~ ", paste0(Xoutcome, collapse = " + "))))
  X.psi <- list(as.formula(paste0("~ ", paste0(Xinteraction, collapse = " + "))))
  outcome <- as.numeric(traindata[,which(colnames(traindata) == outcome)])

  # Calculate SE if testdata == NULL
  var = "none"
  if(is.null(testdata)) var = "bootstrap"

  # Fit dWOLS
  if(dWOLSweight == "IPTW"){
    weightiptw <- function(w) 1/w
    mod <- DTRreg(outcome = outcome, blip.mod = X.psi, tf.mod = X.beta, treat.mod = list(trt ~ 1), method = "dwols", data = traindata, treat.mod.man = list(traindata$ps), weight = weightiptw, var.estim = var)
  } else {
    mod <- DTRreg(outcome = outcome, blip.mod = X.psi, tf.mod = X.beta, treat.mod = list(trt ~ 1), method = "dwols", data = traindata, treat.mod.man = list(traindata$ps), var.estim = var)
  }

  # Derive optimal treatment in test
  if(is.null(testdata)){
    scoredwols <- model.matrix(X.psi[[1]], data = traindata) %*% mod$psi[[1]]
    itr <- factor(ifelse(scoredwols > 0, 1, 0))

    # output score = -scoredwols to match other methods i.e. score > 0 means treatment 1 is better, < 0 means treatment 0 is better
    output <- list(fit = data.frame(coef = mod$psi[[1]], SE = sqrt(diag(mod$covmat[[1]]))), score = data.frame(ID = traindata$ID, score_dwols = -scoredwols, itr_dwols = itr))
    colnames(output$score)[2:3] <- c(paste0("score_dwols_", text.outcome), paste0("itr_dwols_", text.outcome))
  } else {
    itr <- factor(ifelse(model.matrix(X.psi[[1]], data = testdata) %*% mod$psi[[1]] > 0, 1, 0))
    output$dhat <- itr
    output$vhat.dhat <- getValue(y = (exp(-testdata[["y"]]))*exp(testdata[["time"]]), t = testdata[["trt"]], d.hat = itr, p.hat = testdata[['ps']], fu = exp(testdata[['time']]))
     # y is arr here, so we need to convert it back to post-index relapse number.
  }

  # If large independent data are provided
  if (!is.null(sim.big)){
    temp.big <- format.countdata(data = sim.big$data, yvar = "mlogarr0001", timevar = "finalpostdayscount", trtvar = "trt",
                                 xcontinuousvars = c("ageatindex_centered", "prerelapse_num", "premedicalcost", "postrelapse_num", "offset", "FUweight"),
                                 xcategoricalvars = c("female", "prevDMTefficacy"), imputation.method = NULL)
    testdata.big <- data.frame(y = temp.big$y, trt = temp.big$trt, time = log(temp.big$time), temp.big$x)
    itr.big <- ifelse(model.matrix(X.psi[[1]], data = testdata.big) %*% mod$psi[[1]] > 0, 1, 0)
    output$dhat.big <- itr.big
    output$v.dhat <- getTrueValue(ds = sim.big$data, d.hat = itr.big, betas = sim.big$betas)
  }

  return(output)
}

itrLIST <- function(traindata, testdata, outcome, treatment, xvar,
                    maxlen = 2L, seed = seed + fold.i, sim.big = NULL){
  #' ITR with list-based DTR
  #' Listdtr assumes higher outcome is better but our outcome is assumed to be better if lower
  #'
  #' @param traindata - training data with both arms; data.frame
  #' @param testdata - test data with both arms (assume two arms) and same variables as traindata; data.frame
  #' @param outcome - outcome variable, assumed to be favorable if lower; string
  #' @param treatment - treatment variable (assume binary); string
  #' @param xvar - X variables; vector of strings
  #' @param maxlen - maximum number of nodes (i.e. if-else statements); integer (ends with "L")
  #' @param seed - randomization seed; integer
  #' @param sim.big a large test set object (independent from data) to get true value function; list object with elements `data` and `betas`
  #'
  #' @return if testdata != NULL: a list with:
  #'              - dhat: the estimated ITR of every observation in `testdata` (1 means arm 1 and 0 means arm 0)
  #'              - vhat.dhat: estimated value of the estimated ITR of every observation in `testdata`
  #'      OR if testdata == NULL: a list with:
  #'              - score: dataframe with ID and optimal ITR (itr_listDTR)
  #' @return In addition, if sim.big != NULL: the list also contains:
  #'              - dhat.big: the estimated ITR of every observation in `sim.big` (1 means arm 1 and 0 means arm 0)
  #'              - v.dhat: true value of the estimated ITR of every observation in `sim.big`

  output <- list()

  # Prepare X variables into a matrix (with dummy categorical variables)
  Xtrain <- as.matrix(traindata[, xvar])
  if(is.null(testdata)){
    Xtest <- Xtrain
  } else {
    Xtest <- as.matrix(testdata[, xvar])
  }

  # Fit the list-DTR model
  listmod <- listdtr(y = traindata[[outcome]],
                     a = traindata[[treatment]],
                     x = Xtrain,
                     stage.x = rep(1, ncol(Xtrain)),
                     seed = seed, maxlen = maxlen)

  # Optimal ITR in test data
  itr <- predict(listmod, xnew = Xtest, stage = 1)

  if(is.null(testdata)){
    output <- list(coef = listmod, score = data.frame(ID = traindata$ID, itr_listDTR = itr))
    colnames(output$score)[2] <- paste0("itr_listDTR", maxlen, "_", outcome)
  } else {
    output$dhat <- itr
    output$vhat.dhat <- getValue(y = (exp(-testdata[["y"]]))*exp(testdata[["time"]]), t = testdata[["trt"]], d.hat = itr, p.hat = testdata[['ps']], fu = exp(testdata[['time']]))
  } # y is arr here, so we need to convert it back to post-index relapse number.

  # If large independent data are provided
  if (!is.null(sim.big)){
    temp.big <- format.countdata(data = sim.big$data, yvar = "mlogarr0001", timevar = "finalpostdayscount", trtvar = "trt",
                                 xcontinuousvars = c("ageatindex_centered", "prerelapse_num", "premedicalcost", "postrelapse_num", "offset", "FUweight"),
                                 xcategoricalvars = c("female", "prevDMTefficacy"), imputation.method = NULL)
    testdata.big <- data.frame(y = temp.big$y, trt = temp.big$trt, time = log(temp.big$time), temp.big$x)
    Xtest.big <- as.matrix(testdata.big[, c(categoricalvars, continuousvars)])
    itr.big <- predict(listmod, xnew = Xtest.big, stage = 1)
    output$dhat.big <- itr.big
    output$v.dhat <- getTrueValue(ds = sim.big$data, d.hat = itr.big, betas = sim.big$betas)
  }

  return(output)
}

itrLuPoisson <- function(traindata1, traindata0, testdata, xvar, sim.big = NULL){

  #' ITR based on Poisson regression
  #'
  #' @param traindata1 - training data but with only arm = 1; data.frame
  #' @param traindata0 - training data but with only arm = 0; data.frame
  #' @param testdata - test data with both arms (assume two arms) and same variables as traindata; data.frame
  #' @param xvar - X variables to be included in the zero-inflated regression; vector of strings
  #' @param sim.big a large test set object (independent from data) to get true value function; list object with elements `data` and `betas`
  #'
  #' @return if testdata != NULL: a list with:
  #'              - dhat: the estimated ITR of every observation in `testdata` (1 means arm 1 and 0 means arm 0)
  #'              - vhat.dhat: estimated value of the estimated ITR of every observation in `testdata`
  #'     OR if testdata == NULL: a list with:
  #'              - fit: dataframe with coefficient of Poisson fit trt=1 (coef1) and trt=0 (coef0), and corresponding SE (SE1 and SE0)
  #'              - score: dataframe with ID, score (score_pois), optimal ITR (itr_pois)
  #' @return In addition, if sim.big != NULL: the list also contains:
  #'              - dhat.big: the estimated ITR of every observation in `sim.big` (1 means arm 1 and 0 means arm 0)
  #'              - v.dhat: true value of the estimated ITR of every observation in `sim.big`

  output <- list()

  # Save traindata ID
  ID <- c(traindata0$ID, traindata1$ID)

  # Separate into training trt = 1, trt = 0, remove extra variables
  traindata1 <- traindata1 %>% dplyr::select(all_of(xvar),y, time)
  traindata0 <- traindata0 %>% dplyr::select(all_of(xvar),y, time)

  # Fit Poisson by treatment arm
  fit1 <- glm(y ~. - time + offset(time), family = "poisson", data = traindata1)
  fit0 <- glm(y ~. - time + offset(time), family = "poisson", data = traindata0)
  beta1.ini <- fit1$coef
  beta0.ini <- fit0$coef
  delta2 <- beta1.ini - beta0.ini

  # Predict score
  if(is.null(testdata)){
    traindata <- rbind(traindata0, traindata1)
    xtot <- traindata %>% dplyr::select(all_of(xvar))
  } else{
    xtot <- testdata %>% dplyr::select(all_of(xvar))
  }
  xtot <- as.matrix(cbind(1, xtot))
  scorepois <- as.numeric(as.matrix(xtot) %*% delta2)

  # Optimal treatment: scorepois > 0 means treatment 0 is preferred, pois < 0 means treatment 1 is prefered
  itr <- factor(ifelse(scorepois >= 0, 0, 1))

  if(is.null(testdata)){
    output <- list(fit = data.frame(coef1 = beta1.ini, SE1 = summary(fit1)$coef[,2], coef0 = beta0.ini, SE0 = summary(fit0)$coef[,2]), score = data.frame(ID = ID, score_lupois = scorepois, itr_lupois = itr))
  } else {
    output$dhat <- itr
    output$vhat.dhat <- getValue(y = testdata[["y"]], t = testdata[["trt"]], d.hat = itr, p.hat = testdata[['ps']], fu = exp(testdata[['time']]))
  }

  # If large independent data are provided
  if (!is.null(sim.big)){
    temp.big <- format.countdata(data = sim.big$data, yvar = "mlogarr0001", timevar = "finalpostdayscount", trtvar = "trt",
                                 xcontinuousvars = c("ageatindex_centered", "prerelapse_num", "premedicalcost", "postrelapse_num", "offset", "FUweight"),
                                 xcategoricalvars = c("female", "prevDMTefficacy"), imputation.method = NULL)
    testdata.big <- data.frame(y = temp.big$y, trt = temp.big$trt, time = log(temp.big$time), temp.big$x)
    xtot.big <- testdata.big %>% dplyr::select(all_of(categoricalvars), all_of(continuousvars))
    xtot.big <- as.matrix(cbind(1, xtot.big))
    scorepois.big <- as.numeric(as.matrix(xtot.big) %*% delta2)
    itr.big <- ifelse(scorepois.big >= 0, 0, 1)
    output$dhat.big <- itr.big
    output$v.dhat <- getTrueValue(ds = sim.big$data, d.hat = itr.big, betas = sim.big$betas)
  }

  return(output)
}

twoarmglmcount.dr <- function(y, x, time, trt, ps, f1.predictor, f0.predictor, error.control = 1e-3, max.iter = 150, tune = c(0.5, 2.0)){
  y <- y*exp(-time)
  x.aug <- cbind(1, x)
  p.aug <- length(x.aug[1, ])

  # Estimate coefficients with NR algorithm
  beta <- rep(0, p.aug)
  mae <- Inf
  iter <- 0
  epsilon.min <- Inf
  while(mae > error.control && iter <= max.iter && epsilon.min > 0){
    eta <- as.numeric(exp(x.aug %*% beta))
    error <- (trt*(y - eta*f0.predictor/2 - f1.predictor/2)*(1 - ps) - (1 - trt)*(eta*y - eta*f0.predictor/2 - f1.predictor/2)*ps)/(eta*ps + (1 - ps))
    score <- colSums(x.aug*error)
    slopewt <- (y + f0.predictor*(trt/ps - 1)/2 + f1.predictor*((1 - trt)/(1 - ps) - 1)/2)*eta*ps*(1 - ps)/(eta*ps + (1 - ps))^2
    slope <- t(x.aug*slopewt) %*% x.aug
    epsilon.min <- eigen(slope)$value[p.aug]
    if(iter == 0) epsilon0 <- epsilon.min + epsilon.min*(epsilon.min < 0) # fixed to all iterations
    beta <- beta + solve(slope + diag(tune[2]*abs(epsilon0), p.aug, p.aug)) %*% score*tune[1] # adding the diagonal matrix to slove potential singualrity issues
    mae <- sum(abs(score))
    iter <- iter + 1
  }

  converge1 <- 1*(mae <= error.control)
  converge2 <- 0

  # If NP did not converge, solve for minimizing the L2-norm (sum of squares) of the score function to avoid inversing the slope inverse (generally slower than NP)
  if(converge1 == 0){
    lossf <- function(beta){
      eta <- as.numeric(exp(x.aug %*% beta))
      eta.max <- max(eta)
      error <- (trt*(y - eta*f0.predictor/2 - f1.predictor/2)*(1 - ps) - (1 - trt)*(eta*y - eta*f0.predictor/2 - f1.predictor/2)*ps)/(eta*ps + (1 - ps))
      score <- colSums(x.aug*error)
      return(sum(score^2))
    }

    initial.value <- lossf(rep(0, p.aug))
    fit <- optim(rep(0, p.aug), fn = lossf) #, control=list(trace=T))
    beta <- fit$par
    converge2 <- 1*(abs(fit$value) < initial.value/100)
  }

  beta <- as.vector(beta)
  eta <- as.numeric(exp(x.aug %*% beta))
  error <- (trt*(y - eta*f0.predictor/2 - f1.predictor/2)*(1 - ps) - (1 - trt)*(eta*y - eta*f0.predictor/2 - f1.predictor/2)*ps)/(eta*ps + (1 - ps))
  score <- colSums(x.aug*error)
  slopewt <- (y + f0.predictor*(trt/ps - 1)/2 + f1.predictor*((1 - trt)/(1 - ps) - 1)/2)*eta*ps*(1 - ps)/(eta*ps + (1 - ps))^2
  slope <- t(x.aug*slopewt) %*% x.aug
  sigma <- solve(slope) %*% (t(x.aug*error^2) %*% x.aug) %*% solve(slope)

  return(list(coef = beta, vcov = sigma, converge = (converge1 + converge2 > 0)))
}

itrLuDR <- function(traindata, testdata, xvar, RCT = FALSE, tree.depth = 2, n.trees = 200, Kfold = 6, B = 3, seed.cf = 3, plot.gbmperf = FALSE, sim.big = NULL){

  #' ITR based on two regressions and contrast regression
  #'
  #' @param traindata - training data with both arms; data.frame
  #' @param testdata - test data with both arms (assume two arms) and same variables as traindata; data.frame
  #' @param xvar - Categorical X variables to be included in the model; vector of strings
  #' @param RCT Whether treatment is randomized. If RCT=T, the PS is the proportion of patients treated with DMF
  #' @param tree.depth Depth of individual trees in boosting (usually 2-3); integer
  #' @param n.trees Maximum number of trees in boosting (usually 100-1000); integer
  #' @param Kfold Number of folds (parts) used in cross-fitting to partition the data; integer
  #' @param B Number of time cross-fitting is repeated to reduce Monte Carlo variability; integer
  #' @param seed.cf Randomization seed for cross-fitting partitions; integer
  #' @param plot.gbmperf Plot the performance measures in the GBM method; boolean
  #' @param sim.big a large test set object (independent from data) to get true value function; list object with elements `data` and `betas`
  #'
  #' @return if testdata != NULL: a list with two elements, tworeg and contrast reg, each with:
  #'              - dhat: the estimated ITR of every observation in `testdata` (1 means arm 1 and 0 means arm 0)
  #'              - vhat.dhat: estimated value of the estimated ITR of every observation in `testdata`
  #'     OR if testdata == NULL: a list with two elements, tworeg and contrast reg, each with:
  #'               - fit: coefficient of log(CATE) (coef) and SE (for contrast reg only)
  #'               - score: dataframe with ID, score (score_tworeg or score_contrastreg), optimal ITR (itr_tworeg or itr_contrastreg)
  #' @return In addition, if sim.big != NULL: a list with two elements, tworeg and contrast reg, each with:
  #'              - dhat.big: the estimated ITR of every observation in `sim.big` (1 means arm 1 and 0 means arm 0)
  #'              - v.dhat: true value of the estimated ITR of every observation in `sim.big`

  output <- list()

  # Save traindata ID
  ID <- traindata$ID

  # Prepare data
  if(RCT == FALSE){
    traindata_ps <- traindata %>%
      dplyr::select(trt, all_of(xvar))
    # Formula for PS model in cross fitting
    fps <- as.formula(paste("trt ~ ", paste0(xvar, collapse = "+")))
  }
  traindata <- traindata %>% dplyr::select(all_of(xvar), y, trt, time)

  # Prepare for cross-fitting
  N1 <- ifelse(is.numeric(traindata$trt), sum(traindata$trt), sum(as.numeric(traindata$trt) - 1))
  N0 <- nrow(traindata) - N1
  N <- N1 + N0
  p.aug <- ncol(traindata %>% dplyr::select(-y, -trt, -time)) + 1

  index1 <- rep(1:Kfold, floor(N1/Kfold))
  if(N1 > Kfold*floor(N1/Kfold)) index1 <- c(index1, 1:(N1 - Kfold*floor(N1/Kfold)))

  index0 <- rep(1:Kfold, floor(N0/Kfold))
  if(N0 > Kfold*floor(N0/Kfold)) index0 <- c(index0, Kfold + 1 - 1:(N0 - Kfold*floor(N0/Kfold)))

  delta3.mat <- delta4.mat <- matrix(NA, B, p.aug)
  sigma4.mat <- matrix(0, p.aug, p.aug)

  # Cross fitting
  converge <- rep(NA, B)
  for(bb in 1:B){
    cat("\n   Bootstrap:", bb, "out of", B, "\n")
    set.seed(bb + seed.cf)
    index1cv <- sample(index1, N1, F)
    index0cv <- sample(index0, N0, F)
    index <- rep(NA, N)
    index[traindata$trt == 1] <- index1cv
    index[traindata$trt == 0] <- index0cv

    f1.predictcv <- f0.predictcv <- pscv <- rep(NA, N)
    for(k in 1:Kfold){
      datatot_train <- traindata[index != k, ]
      if(RCT == FALSE) datatot_train_ps <- traindata_ps[index != k, ]
      trt_train <- datatot_train$trt[index != k]

      datatot_test <- traindata[index == k, ]
      if(RCT == FALSE) datatot_test_ps <- traindata_ps[index == k, ]
      x_test <- datatot_test %>% dplyr::select(all_of(xvar))

      data1 <- datatot_train %>% filter(trt == 1) %>% dplyr::select(-trt)
      set.seed(100)
      fit1.gbm <- gbm::gbm(y ~. - time + offset(time), data = data1,
                           distribution = "poisson", interaction.depth = tree.depth, n.trees = n.trees, cv.folds = 5)
      best1.iter <- max(10, gbm.perf(fit1.gbm, method = "cv", plot.it = TRUE))
      withCallingHandlers({
        f1.predictcv[index == k] <- predict(object = fit1.gbm, newdata = datatot_test, n.trees = best1.iter, type = "response")
      },
      warning=function(w) {
        if (grepl("does not add the offset", conditionMessage(w)))
          invokeRestart("muffleWarning")  # suppress warning: "predict.gbm does not add the offset to the predicted values."
      })

      data0 <- datatot_train %>% filter(trt == 0) %>% dplyr::select(-trt)
      set.seed(100)
      fit0.gbm <- gbm(y ~. - time + offset(time), data = data0, distribution = "poisson", interaction.depth = tree.depth, n.trees = n.trees, cv.folds = 5)
      best0.iter <- max(10, gbm.perf(fit0.gbm, method = "cv", plot.it = TRUE))
      withCallingHandlers({
        f0.predictcv[index == k] <- predict(object = fit0.gbm, newdata = datatot_test, n.trees = best0.iter, type = "response")
      },
      warning=function(w) {
        if (grepl("does not add the offset", conditionMessage(w)))
          invokeRestart("muffleWarning")  # suppress warning: "predict.gbm does not add the offset to the predicted values."
      })

      if(RCT == FALSE){
        pscv[index == k] <- IPTWfun(PSmodel = fps, data = datatot_train_ps, newdata = datatot_test_ps)$ps
      } else {
        pscv[index == k] <- sum(datatot_train$trt == 1)/nrow(datatot_train)
      }
    }

    ## bb-th cross fitting two regression estimator
    xb <- as.matrix(traindata %>% dplyr::select(all_of(xvar)))

    ## bb-th cross fitting contrast regression estimator
    fit_two <- twoarmglmcount.dr(y = traindata$y, x = xb, time = traindata$time, trt = traindata$trt, ps = pscv, f1.predictor = f1.predictcv, f0.predictor = f0.predictcv)
    delta4.mat[bb, ] <- fit_two$coef
    converge[bb] <- fit_two$converge
    if(converge[bb] == TRUE) sigma4.mat <- sigma4.mat + fit_two$vcov
  }
  # Final contrast regression estimator
  converge4 <- (sum(converge) > 0)
  if(converge4 == TRUE){
    delta4 <- colMeans(delta4.mat[converge == TRUE, , drop = FALSE])
    sigma4 <- sigma4.mat/sum(converge)
  } else {
    delta4 <- colMeans(delta4.mat)
    sigma4 <- sigma4.mat
  }

  # Predict score in test data
  if(is.null(testdata)){
    xtot <- traindata %>% dplyr::select(all_of(xvar))
  } else{
    xtot <- testdata %>% dplyr::select(all_of(xvar))
  }
  xtot <- as.matrix(cbind(1, xtot))
  scorecontrastreg <- as.numeric(as.matrix(xtot) %*% delta4)

  # Optimal treatment: scorepois > 0 means treatment 0 is preferred, pois < 0 means treatment 1 is prefered
  itr_contrastreg <- factor(ifelse(scorecontrastreg >= 0, 0, 1))

  if(is.null(testdata)){
    output <- list(contrastreg = list(fit = data.frame(coef = delta4, SE = sqrt(diag(sigma4.mat))), score = data.frame(ID = ID, score_contrastreg = scorecontrastreg, itr_contrastreg = itr_contrastreg)))
  } else {
    output$valueContrastReg$dhat <- itr_contrastreg
    output$valueContrastReg$vhat.dhat <- getValue(y = testdata[['y']], t = testdata[['trt']], d.hat = itr_contrastreg, p.hat = testdata[['ps']], fu = exp(testdata[['time']]))
  }

  # If large independent data are provided
  if (!is.null(sim.big)){
    temp.big <- format.countdata(data = sim.big$data, yvar = "mlogarr0001", timevar = "finalpostdayscount", trtvar = "trt",
                                 xcontinuousvars = c("ageatindex_centered", "prerelapse_num", "premedicalcost", "postrelapse_num", "offset", "FUweight"),
                                 xcategoricalvars = c("female", "prevDMTefficacy"), imputation.method = NULL)
    testdata.big <- data.frame(y = temp.big$y, trt = temp.big$trt, time = log(temp.big$time), temp.big$x)
    xtot.big <- testdata.big %>% dplyr::select(all_of(categoricalvars), all_of(continuousvars))
    xtot.big <- as.matrix(cbind(1, xtot.big))
    scorecontrastreg.big <- as.numeric(as.matrix(xtot.big) %*% delta4)
    itr_contrastreg.big <- ifelse(scorecontrastreg.big >= 0, 0, 1)
    output$valueContrastReg$dhat.big <- itr_contrastreg.big
    output$valueContrastReg$v.dhat <- getTrueValue(ds = sim.big$data, d.hat = itr_contrastreg.big, betas = sim.big$betas)
  }

  return(output)
}

getValue <- function(y, t, d.hat, p.hat, fu){

  #' Estimated value function of an estimated ITR, V.hat(d.hat), the weighted approach
  #'
  #' @param y - Outcome
  #' @param t - Treatment received 0/1
  #' @param d.hat - Estimated optimal treatment 0/1
  #' @param p.hat - Estimated propensity score, P(T=1|X)
  #' @param fu - Follow-up time

  pt <- t*p.hat + (1-t)*(1-p.hat) # P(T=t|X)
  denom <- sum((t == d.hat)*fu/pt)
  num <- sum(y*(t == d.hat)/pt)
  return(list(U = num, W = denom))
}

getTrueOptimalValue <- function(n, seed){

  #' True value function of the true optimal ITR, V(d)
  #' An empirical way to calculate the true optimal value given the setup in simdata() function
  #'
  #' @param n sample size (choose something large to get close to the theoretical truth); scalar
  #' @param seed randomization seed to generate the simulated data; scalar
  #' @return true optimal value function, $V = E[Y^(optT)]$

  sim <- simcountdata(n = n,
                      seed = seed,
                      beta = c(log(0.4), log(0.5), log(1), log(1.1), log(1.2)),
                      beta.x = c(-1.54, -0.01, 0.06, 0.25, 0.5, 0.13, 0.0000003)
  )

  ds <- sim$data %>%
    mutate(trueT = ifelse(as.numeric(Iscore) < 3, 1, 0),
           trueT = ifelse(as.numeric(Iscore) == 3, trt, trueT)) # true optimal T
    # If all 5 groups
    # optimal T is 1 if in the score group 1 and 2 (high and moderate responder to 1);
    # optimal T is either 1 or 0 if in score group 3 (neutral group);
    # optimal T is 0 if in score group 4 & 5 (high and moderate responder to 0)

  xmat.rate <- model.matrix(~ ageatindex_centered + female + prerelapse_num + prevDMTefficacy + premedicalcost +
                              Iscore + trueT + trueT*Iscore, ds) %>% as.matrix()
  rate <- exp(xmat.rate %*% sim$betas) # Not FU, this is E[Y|time=1]

  return(mean(rate)) # empirical estimate over a large number of samples
}

getTrueWorstValue <- function(n, seed){

  #' True value function of the true worse ITR
  #' An empirical way to calculate the true worse value given the setup in simdata() function
  #' The worse value is the opposite of the true optimal value
  #'
  #' @param n sample size (choose something large to get close to the theoretical truth); scalar
  #' @param seed randomization seed to generate the simulated data; scalar
  #' @return true optimal value function, $V = E[Y^(optT)]$

  sim <- simcountdata(n = n,
                      seed = seed,
                      beta = c(log(0.4), log(0.5), log(1), log(1.1), log(1.2)),
                      beta.x = c(-1.54, -0.01, 0.06, 0.25, 0.5, 0.13, 0.0000003)
  )

  ds <- sim$data %>%
      mutate(trueWorstT = ifelse(as.numeric(Iscore) < 3, 0, 1),
             trueWorstT = ifelse(as.numeric(Iscore) == 3, trt, trueWorstT)) # true worst T
    # If all 5 groups
    # worst T is 0 if in the score group 1 and 2 (high and moderate responder to T=1);
    # worst T is either 0 or 1 if in score group 3 (neutral group);
    # worst T is 1 if in score group 4 & 5 (high and moderate responder to T=0)

  xmat.rate <- model.matrix(~ ageatindex_centered + female + prerelapse_num + prevDMTefficacy + premedicalcost +
                              Iscore + trueWorstT + trueWorstT*Iscore, ds) %>% as.matrix()
  rate <- exp(xmat.rate %*% sim$betas) # Not FU, this is E[Y|time=1]

  return(mean(rate)) # empirical estimate over a large number of samples
}
