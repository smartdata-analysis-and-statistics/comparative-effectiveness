rm(list=ls(all=TRUE))

## Simulate data
source('Chapter 10 - Dealing with irregular and informative visits/sim.r')
set.seed(9843626)

## Simulate a dataset in which age and sex affect the chances of being observed for the outcome
## We also need confounding
## A center effect
## An individual effect

## I created the new function censor_visits_a5 where informative visits depend on age, sex and treatment x

## 1. first, create a simulated dataset in which sex has been added as a variable

dataset  <- sim_data_EDSS(npatients = 500, #10000, # 2000, # previous 50000
                          ncenters = 10,
                          follow_up = 12*5, # Total follow-up (number of months)
                          sd_a_t = 0.5,   # DGM - Within-visit variation in EDSS scores
                          baseline_EDSS = 1.3295,    # DGM - Mean baseline EDDS score
                          sd_alpha_ij = 1.46,    # DGM - Between-subject variation in baseline EDSS
                          sd_beta1_j = 0.20,    # DGM - Between-site variation in baseline EDSS
                          mean_age = 42.41,
                          sd_age = 10.53,
                          min_age = 18,
                          beta_age = 0.05, # DGM - prognostic effect of age
                          beta_t = 0.014,  # DGM - prognostic effect of time
                          beta_t2 = 0,    # DGM - prognostic effect of time squared
                          delta_xt = 0, # DGM - interaction treatment time
                          delta_xt2 = 0, # 0.0005    # DGM - interaction treatment time2
                          p_female = 0.75,
                          beta_female = -0.2 ,  ## DGM - prognostic effect of male sex
                          delta_xf = 0.1,      ## DGM - interaction sex treatment
                          rho = 0.8,             # DGM - autocorrelation of between alpha_tij
                          corFUN = corAR1,       # DGM - correlation structure of the latent EDSS scores
                          tx_alloc_FUN = treatment_alloc_confounding_v2 ) ## or treatment_alloc_randomized

summary(dataset)
head(dataset)
dataset$ages<- (dataset$age)^2
ps<- predict( glm( dataset$x ~ dataset$age + dataset$ages + dataset$sex, family='binomial'), type='response')

dataset$ipt<- ifelse( dataset$x==1, 1/ps, 1/(1-ps))

## 2. then remove Y according to the informative visit process

dataset_visit <- censor_visits_a5(dataset)
summary(dataset_visit)
head(dataset_visit)

# 1.46 conditional effect adjusted for confounding (if no interaction trmt and age, sex)
lm( dataset$y ~ dataset$x + dataset$age + dataset$ages + dataset$sex )

#  1.46  marginal effect (this confirms that conditional effect is rather a good estimate here)
lm(dataset$y ~ dataset$x, weight=dataset$ipt)

# 1.26
# estimate not adjusted for confounding
lm( dataset$y ~ dataset$x )

ps2<- predict( glm( dataset_visit$x ~ dataset_visit$age +  dataset_visit$ages + dataset_visit$sex, family='binomial'), type='response')
dataset_visit$ipt2<- ifelse( dataset_visit$x==1, 1/ps2, 1/(1-ps2))

# 1.39
# estimate not accounting for the visit process (but acc. for confounding)
lm( dataset_visit$y_obs ~ dataset_visit$x, weight= dataset_visit$ipt2) ## different estimate in the observed dataset

# 1.19
# not adjusted for anything:
lm( dataset_visit$y_obs ~ dataset_visit$x) ## different estimate in the observed dataset




#######################################################################################################

## Estimation of the marginal effect of treatment

#######################
## doubly -weighted  ##
#######################

library(survival)


## visit process

dataset_visit$t1<-dataset_visit$time-1
dataset_visit$t2<-dataset_visit$time
dataset_visit$visit<- 1
dataset_visit[is.na(dataset_visit$y_obs),]$visit<- 0
#gamma<-glm( dataset_visit$visit  ~dataset_visit$x + dataset_visit$sex + dataset_visit$age , family='binomial',data=dataset_visit)$coef
#gamma<-coxph(Surv(dataset_visit$t1, dataset_visit$t2, dataset_visit$visit )~dataset_visit$x + dataset_visit$sex + dataset_visit$age , data=dataset_visit)$coef
#dataset_visit$rho_i<- 1/exp( gamma[1]*dataset_visit$x + gamma[2]*dataset_visit$sex + gamma[3]*dataset_visit$age  )

# try fitting gamma when removing the time=0 that is deterministic
dataset_visit2<- dataset_visit[dataset_visit$time!=0,]
gamma<-glm( dataset_visit2$visit  ~dataset_visit2$x + dataset_visit2$sex + dataset_visit2$age , family='binomial',data=dataset_visit2)$coef
dataset_visit$rho_i<- 1/ exp(gamma[1]+ gamma[2]*dataset_visit$x + gamma[3]*dataset_visit$sex + gamma[4]*dataset_visit$age  )
boxplot(dataset_visit$rho_i, dataset_visit$visit)

## confounding
memory.limit(10000000)

ps<-predict(glm( dataset_visit$x ~ dataset_visit$age +  dataset_visit$ages + dataset$sex, family='binomial'),type='response')
png('C:\\Users\\Janie\\Documents\\University\\McGill\\7-Other publications\\9- Book Chapter Debray Simoneau\\resultsorfigures\\boxplotipt.png')
boxplot(ps~ dataset_visit$x)
dev.off()
dataset_visit$ipt<- ifelse( dataset_visit$x==1, 1/ps, 1/(1-ps))


lm( dataset_visit$y_obs ~ dataset_visit$x, weights= ipt*rho_i, data=dataset_visit)
# 1.47


#######################
## Imputation method ##
#######################

library(mice)

source("mlmi.r")

impute <- function(dat,
                   times = seq(0, 60, by = 3),
                   n.imp = 10,
                   maxit = 5,
                   outcome_var = "edss",
                   treat_var = "trt",
                   time_var = "time") {

  dat$time <- dat[,time_var]
  dat$trt <- dat[,treat_var]

  impdat <- data.frame(patid = rep(unique(dat$patid), each = length(times)),
                       time = rep(times, length(unique(dat$patid))))

  impdat <- bind_rows(dat, impdat)

  impdat <- impdat %>% group_by(patid) %>% arrange(time, .by_group = TRUE) %>% summarize(
    centerid = first(na.omit(centerid)),
    patid = first(patid),
    time = unique(time),
    trt = first(na.omit(trt)),
    age = first(na.omit(age)),
    sex = first(na.omit(sex)))


  # add EDSS scores
  impdat <- full_join(impdat, dat, by = c("patid" = "patid", "time" = "time")) %>%
    dplyr::select(centerid.x, patid, trt.x, age.x, sex.x, time, edss)
  colnames(impdat) <- c("centerid", "patid", "trt", "age", "sex", "time", "edss")

  # Add interaction term between time and treatment
  impdat$trttime <- impdat$time * impdat$trt

  # Add interaction term between sex and treatment
  impdat$trtsex <- impdat$sex * impdat$trt

  # Prepare imputation method
  imp0 <- mice(impdat, maxit = 0)

  #imputation Method
  imeth <- imp0$method
  imeth[outcome_var] <- "mlmi"

  predM <- imp0$predictorMatrix
  predM[outcome_var, "centerid"] <- -2
  predM[outcome_var, "patid"] <- -3
  predM[outcome_var, time_var] <- 6
  diag(predM) <- 0

  # Run the imputation
  imp <-  mice(impdat, predictorMatrix = predM, method = imeth, maxit = maxit, m = n.imp)

  return(imp)
}

dataa<- data.frame(dataset_visit)
dataa<- dataa[,-c(9:18)]
colnames(dataa)[3]<- 'trt'
dataa<- dataa[,-c(8)]
head(dataa)
data_imputed<- impute(             dataa ,
                                   times = seq(0, 60, by = 3),
                                   n.imp = 10,
                                   maxit = 5,
                                   outcome_var = "edss",
                                   treat_var = "trt",
                                   time_var = "time")

head(data_imputed)
data_imp1<- complete(data_imputed, action=1)
head(data_imp1)
data_imp1$ages<- (data_imp1$age)^2
ps<- predict( glm( data_imp1$trt~ data_imp1$age + data_imp1$ages + data_imp1$sex,family='binomial'), type='response')
data_imp1$ipt<- ifelse( data_imp1$trt ==1, 1/ps, 1/(1-ps))
lm(data_imp1$edss ~ data_imp1$trt, weight=ipt, data=data_imp1)
# 1.46
