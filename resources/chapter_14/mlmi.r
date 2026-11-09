mice.impute.mlmi <- function(y, ry, x, type, wy = NULL, intercept = TRUE,
                             donors = 5L,
                             conversion = "pmm",
                             matchtype = 1L, ...)
{

  wy <- !ry
  x <- cbind(1, as.matrix(x))
  type <- c(2, type)
  names(type)[1] <- colnames(x)[1] <- "(Intercept)"

  clust <- names(type[type == -3]) # Clustering of patients
  rande <- names(type[type == 2]) # Random slopes - these are ignored
  group <- names(type[type == -2]) # Clustering of centers
  time <- names(type[type == 6]) # Time will be treated as a fixed effect and also be used to define the autocorrelation

  fixe <- names(type[type > 0])
  lev <- unique(x[, clust])
  X <- x[, fixe, drop = FALSE]
  Z <- x[, rande, drop = FALSE]
  xobs <- x[ry, , drop = FALSE]
  yobs <- y[ry]
  Xobs <- X[ry, , drop = FALSE]
  Zobs <- Z[ry, , drop = FALSE]

  fixmodel <- paste("yobs ~ ", paste(fixe[-1L], collapse = "+"))

  suppressWarnings(fit <- try(lme(formula(fixmodel),
                                  random =  formula(paste("~1|", group, "/", clust)),
                                  correlation = corExp(form = formula(paste("~", time, "|", group, "/", clust))),
                                  data = data.frame(yobs, xobs),
                                  control = list(returnObject = TRUE))))


  if (("try-error" %in% class(fit))) {
    stop("Estimation of multilevel model failed!")
  }

  # Extract the treatment effect
  #print(fixef(fit)["x"])

  vnames <- names(fixef(fit))

  # Store the observed predictions
  yhatobs <- predict(fit, level = 2)
  yimp_raw <- yimp_pmm <- y

  sigmahat <- fit$sigma # fit_intervals$sigma["est."]
  df <- nrow(fit$data) - length(fixef(fit))

  rancoef_group <- as.matrix(ranef(fit)[[group]])
  lambda_group <- t(rancoef_group) %*% rancoef_group
  psi_group_star <- lambda_group/rchisq(1, df = nrow(rancoef_group))

  rancoef_clust <- as.matrix(ranef(fit)[[clust]])
  lambda_clust <- t(rancoef_clust) %*% rancoef_clust
  psi_clust_star <- lambda_clust/rchisq(1, df = nrow(rancoef_clust))

  cs <- fit$modelStruct$corStruct
  range_star <- as.numeric(coef(cs, unconstrained = FALSE))


  fit_intervals <- try(intervals(fit), silent = TRUE)

  if (!("try-error" %in% class(fit_intervals))) {
    # Draw random sample for the range of the autocorrelation matrix
    se_log_range = (log(fit_intervals$corStruct[, "upper"]) - log(fit_intervals$corStruct[, "lower"])) / (2*qnorm(0.975))
    range_star <- exp(rnorm(1, mean = log(fit_intervals$corStruct[, "est."]), sd = se_log_range))

    se_log_tau_center <- (log(((fit_intervals$reStruct[[group]])["sd((Intercept))","upper"])) - log(((fit_intervals$reStruct[[group]])["sd((Intercept))","lower"])))/(2*qnorm(0.975))
    se_log_tau_patient <- (log(((fit_intervals$reStruct[[clust]])["sd((Intercept))","upper"])) - log(((fit_intervals$reStruct[[clust]])["sd((Intercept))","lower"])))/(2*qnorm(0.975))

    # Derive SE for the between-study SD directly from lme4
    tau_patient_star <- exp(rnorm(n = 1, mean = log((fit_intervals$reStruct[[clust]])["sd((Intercept))","est."]), sd = se_log_tau_patient))
    psi_clust_star <- tau_patient_star**2
  }

  psi_star <- as.numeric(psi_group_star + psi_clust_star)


  # Draw a random sample for the residual error
  sigma2star <- df * sigmahat^2/rchisq(n = 1, df = df)

  # Rescale the covariance matrix to the new draw of sigma
  covmat <- sigma2star * (vcov(fit)/sigmahat^2)
  rv <- t(chol(covmat))

  # Draw random sample for the beta coefficients
  beta_star <- fixef(fit) + rv %*% rnorm(ncol(rv))
  rownames(beta_star) <- vnames



  # Iterate over all patients
  for (jj in seq_along(lev)) {
    Xi <- as.matrix(Xobs[xobs[, clust] == lev[jj], ])

    if (sum(xobs[,clust] == jj)  <= 1) {
      # If we only have one observation, Xi is a single-column matrix but should be a single-row
      Xi <- t(Xi)
    }

    # Check if we need to impute anything
    if (sum(is.na(y[x[,clust] == jj])) == 0) {
      next
    }

    # Identify all relevant time points
    spatDat <- data.frame(time = sort(unique(x[x[,clust] == jj, time])))

    cs1Exp <- corExp(value = range_star, form = ~ time)
    cs1Exp <- Initialize(cs1Exp, spatDat)
    sigma_full <- cor2cov(sd = rep(sqrt(sigma2star), nrow(spatDat)), rho = corMatrix(cs1Exp))

    # Population-level predictions
    blup_pop_beta <- x[x[,clust] == jj, vnames] %*% beta_star[vnames,]

    yi <- yobs[xobs[, clust] == jj]
    ti <- as.matrix(Xobs[xobs[, clust] == jj, time])
    Zi <- as.matrix(Zobs[xobs[, clust] == jj, ])
    idx <- spatDat$time %in% ti

    sigma_all_pat <- sigma_full[idx, idx]

    psistar_zi <- psi_star %*% t(Zi)
    zi_psistar_zi <- Zi %*% psistar_zi

    resid_pop_blup <- (yi - Xi[,vnames] %*% beta_star[vnames,])

    Mi <- psistar_zi %*% chol2inv(chol(zi_psistar_zi + sigma_all_pat))
    myi <- Mi %*% resid_pop_blup # Mean of the conditional random effect for each cluster
    vyi <- psi_star - Mi %*% Zi * psi_star

    # Draw random draws for the random effects
    bi_star <- myi + sqrt(vyi) * rnorm(length(myi))

    # Derive the linear predictor on the cluster level
    blup_clust_beta <- blup_pop_beta + rep(bi_star, each = nrow(spatDat))

    y_pred_resid <- y[x[,clust] == jj] - blup_clust_beta
    y_imp_resid <- y_pred_resid

    # Identify observed and missing residuals
    dependent.ind <- which(is.na(y_pred_resid))
    given.ind <- which(!is.na(y_pred_resid))

    B <- sigma_full[dependent.ind, dependent.ind]
    C <- sigma_full[dependent.ind, given.ind, drop = FALSE]
    D <- sigma_full[given.ind, given.ind]

    CDinv <-  C %*% chol2inv(chol(D))

    # Derive the conditional residual errors
    cMu <- as.vector(CDinv %*% y_pred_resid[given.ind])
    cVar <- B - CDinv %*% t(C)

    y_pred_resid[dependent.ind] <- cMu
    y_imp_resid[dependent.ind] <-  MASS::mvrnorm(n = 1, mu = cMu, Sigma = cVar)

    ## Apply predictive mean matching
    yobsi <- yobs[xobs[, clust] == jj]
    yhatobsi <- yhatobs[xobs[, clust] == jj]
    yhatmisi <- blup_clust_beta[wy[x[, clust] == jj]] + y_pred_resid[wy[x[, clust] == jj]] #Note: we here use the expected conditional residuals
    idx_pmm <- mice:::matchindex(d = yhatobsi, t = yhatmisi, k = 5L)


    # Identify which rows in yimp_raw need to be replaced
    index_pat <- which(x[,clust] == jj)
    index_missing <- which(is.na(y))

    ## Store the raw imputations
    #yimp_raw[index_pat[index_pat %in% index_missing]] <-  blup_clust_beta[wy[x[, clust] == jj]] + y_imp_resid[dependent.ind] # We add "random" conditional residuals
    yimp_pmm[index_pat[index_pat %in% index_missing]] <- yobsi[idx_pmm]
  }


  # Return the imputed values
  yimp_pmm[wy]
}


### Fast implementation of MICE algorithm
## Residual errors can be sampled using two appraoches:
# If sampling_resid is set to "condMVN", imputed values are generated by accounting
# for autocorrelation in the observed responses
# If sampling_resid is set to "MVN", imputed values are generated independently of
# the observed responses
impute_y_mice_3l <- function(data,
                             n.imp = 10,
                             seed,
                             sampling_resid = "condMVN" # Method for sampling the residual errors
) {
  require(Matrix)

  set.seed(seed)

  ####################### PREP
  y <- data$y_obs
  yimp_raw <- yimp_pmm <- replicate(n.imp, y)
  ry <- !(is.na(data$y_obs))
  x <- data[,c("patid", "centerid", "time", "time_x", "age", "sex","edss")]
  type <- c("patid" = -3, "centerid" = -2, "time"  = 6, "time_x" = 1,
            "age" = 1,
            "sex" = 1,
            "edss" = 1)


  ####################### START
  wy <- !ry
  x <- cbind(1, as.matrix(x))
  type <- c(2, type)
  names(type)[1] <- colnames(x)[1] <- "(Intercept)"

  clust <- names(type[type == -3]) # Patient
  rande <- names(type[type == 2])
  group <- names(type[type == -2]) # Center
  time <- names(type[type == 6]) # Time

  fixe <- names(type[type > 0])
  lev <- unique(x[, clust])
  X <- x[, fixe, drop = FALSE]
  Z <- x[, rande, drop = FALSE]
  xobs <- x[ry, , drop = FALSE]
  yobs <- y[ry]
  Xobs <- X[ry, , drop = FALSE]
  Zobs <- Z[ry, , drop = FALSE]

  fixmodel <- paste("yobs ~ ", paste(fixe[-1L], collapse = "+"))

  fit <- try(lme(formula(fixmodel),
                 random =  formula(paste("~1|", group, "/", clust)),
                 correlation = corExp(form = formula(paste("~", time, "|", group, "/", clust))),
                 data = data.frame(yobs, xobs),
                 control = list(returnObject = TRUE)))


  if (("try-error" %in% class(fit))) {
    warning("Estimation of multilevel model failed!")
    return(list(data = data, fit = NULL))
  }

  vnames <- names(fixef(fit))

  # Store the observed predictions
  yhatobs <- predict(fit, level = 2)

  # Derive standard error for key model parameters
  fit_intervals <- try(intervals(fit))
  sigmahat <- fit$sigma # fit_intervals$sigma["est."]
  df <- nrow(fit$data) - length(fixef(fit))

  if (!("try-error" %in% class(fit_intervals))) {
    se_log_tau_center <- (log(((fit_intervals$reStruct[[group]])["sd((Intercept))","upper"])) - log(((fit_intervals$reStruct[[group]])["sd((Intercept))","lower"])))/(2*qnorm(0.975))
    se_log_tau_patient <- (log(((fit_intervals$reStruct[[clust]])["sd((Intercept))","upper"])) - log(((fit_intervals$reStruct[[clust]])["sd((Intercept))","lower"])))/(2*qnorm(0.975))

    # Draw random sample for the cluster effects
    tau_center_star <- exp(rnorm(n = n.imp, mean = log((fit_intervals$reStruct[[group]])["sd((Intercept))","est."]), sd = se_log_tau_center))

    # Draw random sample for the patient effects
    tau_patient_star <- exp(rnorm(n = n.imp, mean = log((fit_intervals$reStruct[[clust]])["sd((Intercept))","est."]), sd = se_log_tau_patient))

    psi_star <- (tau_center_star**2) + (tau_patient_star**2)

    # Draw random sample for the range of the autocorrelation matrix
    se_log_range = (log(fit_intervals$corStruct[, "upper"]) - log(fit_intervals$corStruct[, "lower"])) / (2*qnorm(0.975))
    range_star <- exp(rnorm(n.imp, mean = log(fit_intervals$corStruct[, "est."]), sd = se_log_range))

  } else {
    warning("Error when estimating the confidence intervals for the multilevel model")

    # Sometimes, LME is not able to determine the variance of random effect parameters
    # We then derive an approximation using methods described by Jolani
    rancoef_group <- as.matrix(ranef(fit)[[group]])
    lambda_group <- t(rancoef_group) %*% rancoef_group
    psi_group_star <- rep(lambda_group, n.imp)/rchisq(n.imp, df = nrow(rancoef_group))

    rancoef_clust <- as.matrix(ranef(fit)[[clust]])
    lambda_clust <- t(rancoef_clust) %*% rancoef_clust
    psi_clust_star <- rep(lambda_clust, n.imp)/rchisq(n.imp, df = nrow(rancoef_clust))

    psi_star <- psi_group_star + psi_clust_star

    # For the estimated range parameter, we ignore uncertainty if the variance could not be estimated
    cs <- fit$modelStruct$corStruct
    range_est <- as.numeric(coef(cs, unconstrained = FALSE))
    range_star <- rep(range_est, n.imp)
  }

  spatDat <- data.frame(time = 0:max(x[,time]))



  pbimp <- txtProgressBar(min = 0, max = n.imp, style = 3)
  for (imp in 1:n.imp) {
    # Draw a random sample for the residual error
    sigma2star <- df * sigmahat^2/rchisq(n = 1, df = df)

    # Rescale the covariance matrix to the new draw of sigma
    covmat <- sigma2star * (vcov(fit)/sigmahat^2)
    rv <- t(chol(covmat))

    # Draw random sample for the beta coefficients
    beta_star <- fixef(fit) + rv %*% rnorm(ncol(rv))
    rownames(beta_star) <- vnames


    cs1Exp <- corExp(value = range_star[imp], form = ~ time)
    cs1Exp <- Initialize(cs1Exp, spatDat)
    sigma_full <- cor2cov(sd = rep(sqrt(sigma2star), nrow(spatDat)), rho = corMatrix(cs1Exp))

    # Population-level predictions
    blup_pop_beta <- X[,vnames] %*% beta_star[vnames,]

    # Patient-level predictions
    blup_clust_beta <- y_pred_resid <- rep(0, nrow(blup_pop_beta))

    # Set up sigma for all clusters
    sigma_all_pat <- list()
    zi_psistar_zi <- list()
    psistar_zi <- list()
    resid_pop_blup <- list()
    zi_all_pat <- list()

    for (jj in lev) {
      Xi <- as.matrix(Xobs[xobs[, clust] == jj, ])
      if (sum(xobs[,clust] == jj)  <= 1) {
        # If we only have one observation, Xi is a single-column matrix but should be a single-row
        Xi <- t(Xi)
      }

      yi <- yobs[xobs[, clust] == jj]
      ti <- as.matrix(Xobs[xobs[, clust] == jj, time])
      Zi <- as.matrix(Zobs[xobs[, clust] == jj, ])
      idx <- spatDat$time %in% ti

      zi_all_pat[[jj]] <- Zi
      sigma_all_pat[[jj]] <- sigma_full[idx, idx]

      psistar_zi[[jj]] <- psi_star[imp] %*% t(Zi)
      zi_psistar_zi[[jj]] <- Zi %*% psistar_zi[[jj]]

      resid_pop_blup[[jj]] <- (yi - Xi[,vnames] %*% beta_star[vnames,])
    }

    sigma_all_pat <- do.call(bdiag, sigma_all_pat)
    psistar_zi <- do.call(c, psistar_zi)
    zi_psistar_zi <- do.call(bdiag, zi_psistar_zi)
    resid_pop_blup <- do.call(bdiag, resid_pop_blup)
    zi_all_pat <- do.call(bdiag, zi_all_pat)

    Mi <- psistar_zi %*% chol2inv(chol(zi_psistar_zi + sigma_all_pat))
    myi <- Mi %*% resid_pop_blup # Mean of the conditional random effect for each cluster
    vyi <- psi_star[imp] - Mi %*% zi_all_pat * psi_star[imp]

    # Draw random draws for the random effects
    bi_star <- myi + sqrt(vyi) * rnorm(length(myi))

    # Derive the linear predictor on the cluster level
    blup_clust_beta <- blup_pop_beta + rep(bi_star, each = nrow(spatDat))

    y_pred_resid <- y - blup_clust_beta
    y_imp_resid <- y_pred_resid

    # Identify observed and missing residuals
    dependent.ind <- which(is.na(y_pred_resid))
    given.ind <- which(!is.na(y_pred_resid))

    if (sampling_resid == "condMVN") {
      sigma_all_pat <- do.call(bdiag, replicate(length(lev), sigma_full, simplify = F))

      B <- sigma_all_pat[dependent.ind, dependent.ind]
      C <- sigma_all_pat[dependent.ind, given.ind, drop = FALSE]
      D <- sigma_all_pat[given.ind, given.ind]

      CDinv <-  C %*% chol2inv(chol(D))

      # Derive the conditional residual errors
      cMu <- as.vector(CDinv %*% y_pred_resid[given.ind])
      CH_cVar <- Cholesky(B - CDinv %*% t(C))

      y_pred_resid[dependent.ind] <- cMu
      y_imp_resid[dependent.ind] <-  sparseMVN::rmvn.sparse(1, mu = cMu, CH = CH_cVar, prec = FALSE)[1,]
    } else if (sampling_resid == "MVN") {
      # Generate a sample for each patient
      cMu <- rep(0, length(spatDat$time))
      random_resid <- rmvnorm(n = length(lev), mean = cMu, sigma = sigma_full)

      y_pred_resid[dependent.ind] <- rep(0, length(dependent.ind))
      y_imp_resid[dependent.ind] <-  as.vector(t(random_resid))[dependent.ind]
    } else {
      y_pred_resid[dependent.ind] <- rep(0, length(dependent.ind))
      y_imp_resid[dependent.ind] <-  rnorm(length(dependent.ind)) * sqrt(sigma2star)
    }

    ## Store the raw imputations
    #yimp_raw[wy, imp] <-  blup_clust_beta[wy] + y_imp_resid[dependent.ind] # We add "random" conditional residuals

    ## Apply predictive mean matching
    yhatmis <- blup_clust_beta[wy] + y_pred_resid[wy] #Note: we here use the expected conditional residuals
    idx_pmm <- mice:::matchindex(d = yhatobs, t = yhatmis, k = 5L)
    yimp_pmm[wy, imp] <- y[ry][idx_pmm]


    # Update progress bar
    setTxtProgressBar(pbimp, imp)
  }
  close(pbimp)

  #data[, paste("y_imp_lme3l_mi_conv_", seq_along(1:n.imp), sep = "")] <- convert_to_EDSS_scale(yimp_raw[,seq(n.imp)])
  return(yimp_pmm[,seq(n.imp)])
}
