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
    yimp_raw[index_pat[index_pat %in% index_missing]] <-  blup_clust_beta[wy[x[, clust] == jj]] + y_imp_resid[dependent.ind] # We add "random" conditional residuals
    yimp_pmm[index_pat[index_pat %in% index_missing]] <- yobsi[idx_pmm]
  }
  
  
  # Return the imputed values
  yimp_pmm[wy]
}
