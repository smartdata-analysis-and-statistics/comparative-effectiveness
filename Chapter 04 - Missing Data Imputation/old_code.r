

```{r, fig.align="center", fig.width=6, fig.height=6, fig.cap="Figure 1:Distribution of confounders and outcome variable"}
datap <- copy(data_noHTE)
p1<-plot_den(datap=datap,var="age",varlab="Age (years)")  #Age density
p2<-plot_count(datap=datap[, gender:= as.factor(ifelse(female==1,"Female","Male"))],var="gender",varlab="Gender") # Gender proportion
p3<-plot_count(datap=datap[, efficacy := factor(ifelse(prevDMTefficacy=="Low_efficacy","Low",ifelse(prevDMTefficacy=="Medium_high_efficacy","Medium & high","None")),levels=c("None","Low","Medium & high"))],
               var="efficacy",varlab="Gender") # Previuos treatment proportion
p4<-plot_den(datap=datap[,logcost:=log(premedicalcost)],var="logcost",varlab="Log previous cost")  #log cost density
p5<-plot_count(datap=datap,var="prerelapse_num",varlab="Previous relapse number") # Previous relapses proportion
p6<-Relapserate_plot(datahom=data_noHTE,datahet=data_hiHTE)
ggpubr::ggarrange(p1,p2,p3,p4,p5,p6, ncol = 3, nrow = 2)# All plots together
```

## 2.  Estimation of treatment effect in the absence of missing data

```{r}
#1.1 Estimation of treatment effects in the complete data
true_noHTE_ATE <- getest(data = data_noHTE, estimandv = "ATE", Tform = 1, CC = TRUE)
true_noHTE_ATT <- getest(data = data_noHTE, estimandv = "ATT", Tform = 1, CC = TRUE)
true_hiHTE_ATE <- getest(data = data_hiHTE, estimandv = "ATE", Tform = 1, CC = TRUE)
true_hiHTE_ATT <- getest(data = data_hiHTE, estimandv = "ATT", Tform = 1, CC = TRUE)

true<-setDT(rbind(true_noHTE_ATE,true_noHTE_ATT,true_hiHTE_ATE,true_hiHTE_ATT))
true[,Scenario:="Full"]
true[,Treatment:=c(rep("Homogeneus",4),rep("Heterogeneus",4))]
true[,type:= "True"]
```

## 3. Missing data generation 

For generate missing data in covariates we considered the following:
  - Variables collected at patient registration are complete (age, female)
- Variables that makes a reference on MS progression were supossed MAR dependable on the age and female. The missingness on this variables makes that missing process favors the absence of this covariates on patients with moderate and high TERI treatment effect.
- Medical cost variable was missing complete at random (MCAR).
- The previous number of relapses (prerelapse_num) was simulated under the following scenarios, where the missingness of prerelapse_num depends on:
  1. MAR: age and sex favoring absence on patients with moderate and high treament effect
2. MART:  age and sex and also the treatment variable
3. MARTY: age, sex, treatment and Y the outcome variable 
4. MNAR:  age, sex and prerelapse_num 

```{r}
mdata_noHTE <- getmissdata(data_noHTE)
mdata_hiHTE <- getmissdata(data_hiHTE)
```
#4. Data imputation

```{r}
# No output in the prediction model 
form_ny <- list(prevDMTefficacy ~ age+female + years+premedicalcost+numSymptoms+treatment+prerelapse_num,
                premedicalcost ~ age+female + years+prevDMTefficacy+numSymptoms+treatment+prerelapse_num,
                numSymptoms ~ age+female + premedicalcost+years+prevDMTefficacy+prerelapse_num+treatment,
                prerelapse_num ~age+female + premedicalcost+years+prevDMTefficacy+numSymptoms+treatment)
form_ny <- name.formulas(form_ny)
# Output in the prediction model 
form_y <- list(prevDMTefficacy ~ age+female + years+premedicalcost+numSymptoms+treatment+prerelapse_num+y,
               premedicalcost~age+female + years+prevDMTefficacy+numSymptoms+treatment+prerelapse_num+y,
               numSymptoms ~age+female + premedicalcost+years+prevDMTefficacy+prerelapse_num+treatment+y,
               prerelapse_num ~age+female + premedicalcost+years+prevDMTefficacy+numSymptoms+treatment+y)
form_y <- name.formulas(form_y)
# Prediction model with interaction T*X
form_i <- list(prevDMTefficacy ~ treatment*(age+female + years+premedicalcost+numSymptoms+prerelapse_num),
               premedicalcost~ treatment*(age+female + years+prevDMTefficacy+numSymptoms+prerelapse_num),
               numSymptoms ~ treatment*(age+female + premedicalcost+years+prevDMTefficacy+prerelapse_num),
               prerelapse_num ~ treatment*(age+female + premedicalcost+years+prevDMTefficacy+numSymptoms))
form_i <- name.formulas(form_i)
# Prediction model with interaction T*X, T*Y with X*Y!!
form_iy <- list(prevDMTefficacy ~ treatment*(age+female +years+premedicalcost+numSymptoms+prerelapse_num+y)+y*(age+female + years+premedicalcost+numSymptoms+prerelapse_num),
                premedicalcost~ treatment*(age+female + years+prevDMTefficacy+numSymptoms+prerelapse_num+y)+y*(age+female + years+prevDMTefficacy+numSymptoms+prerelapse_num),
                numSymptoms ~ treatment*(age+female + premedicalcost+years+prevDMTefficacy+prerelapse_num+y)+y*(age+female + premedicalcost+years+prevDMTefficacy+prerelapse_num),
                prerelapse_num ~ treatment*(age+female + premedicalcost+years+prevDMTefficacy+numSymptoms+y)+y*(age+female + premedicalcost+years+prevDMTefficacy+numSymptoms) )
form_iy <- name.formulas(form_iy)

# Prediction model with interaction T*X, T*Y 
form_iyt <- list(prevDMTefficacy ~ treatment*(age+female +years+premedicalcost+numSymptoms+prerelapse_num+y),
                 premedicalcost~ treatment*(age+female + years+prevDMTefficacy+numSymptoms+prerelapse_num+y),
                 numSymptoms ~ treatment*(age+female + premedicalcost+years+prevDMTefficacy+prerelapse_num+y),
                 prerelapse_num ~ treatment*(age+female + premedicalcost+years+prevDMTefficacy+numSymptoms+y))
form_iyt <- name.formulas(form_iyt)


imp <- mice(mdata_noHTE[[1]],form=form_ny,m=1,iter=1)
method <- imp$method
method["numSymptoms"] <- "pmm"
method["prevDMTefficacy"] <- "pmm"

#3.1. MICE all groups together---

miceout_noHTEny <- lapply(mdata_noHTE,mice,m=5,form=form_ny,method=method)
miceout_noHTEy <- lapply(mdata_noHTE,mice,m=5,form=form_y,method=method)
miceout_noHTEi <- lapply(mdata_noHTE,mice,m=5,form=form_i,method=method)
miceout_noHTEiy <- lapply(mdata_noHTE,mice,m=5,form=form_iy,method=method)
miceout_noHTEiyt <- lapply(mdata_noHTE,mice,m=5,form=form_iyt,method=method)
miceout_hiHTEny <- lapply(mdata_hiHTE,mice,m=5,form=form_ny,method=method)
miceout_hiHTEy <- lapply(mdata_hiHTE,mice,m=5,form=form_y,method=method)
miceout_hiHTEi <- lapply(mdata_hiHTE,mice,m=5,form=form_i,method=method)
miceout_hiHTEiy <- lapply(mdata_hiHTE,mice,m=5,form=form_iy,method=method)
miceout_hiHTEiyt <- lapply(mdata_hiHTE,mice,m=5,form=form_iyt,method=method)

#3.2. MICE separated by Treatment groups ----
miceout_noHTEs <- lapply(mdata_noHTE,separate_mice,form_y=form_y,method=method)
miceout_hiHTEs <- lapply(mdata_hiHTE,separate_mice,form_y=form_y,method=method)

#3.3. Random forest in mice ----
methodrf<-method
methodrf["premedicalcost"] <- "rf"
methodrf["prerelapse_num"] <- "rf"
methodrf["numSymptoms"] <- "rf"
methodrf["prevDMTefficacy"] <- "rf"
miceout_noHTEyrf <- lapply(mdata_hiHTE,mice,m=5,form=form_y,method=methodrf)
miceout_hiHTEyrf <- lapply(mdata_hiHTE,mice,m=5,form=form_y,method=methodrf)


# 3.4. Classification and regression tree in mice ---
methodcart<-method
methodcart["premedicalcost"] <- "cart"
methodcart["prerelapse_num"] <- "cart"
methodcart["numSymptoms"] <- "cart"
methodcart["prevDMTefficacy"] <- "cart"
miceout_noHTEycart <- lapply(mdata_hiHTE,mice,m=5,form=form_y,method=methodcart)
miceout_hiHTEycart <- lapply(mdata_hiHTE,mice,m=5,form=form_y,method=methodcart)

# 3.5. Missforest ----
for_mdata_noHTE<-list(missForest(mdata_noHTE[[1]][,-c(11:16)])$ximp,
                      missForest(mdata_noHTE[[2]][,-c(11:16)])$ximp,
                      missForest(mdata_noHTE[[3]][,-c(11:16)])$ximp,
                      missForest(mdata_noHTE[[4]][,-c(11:16)])$ximp,
                      missForest(mdata_noHTE[[5]][,-c(11:16)])$ximp)
for_mdata_hiHTE<-list(missForest(mdata_hiHTE[[1]][,-c(11:16)])$ximp,
                      missForest(mdata_hiHTE[[2]][,-c(11:16)])$ximp,
                      missForest(mdata_hiHTE[[3]][,-c(11:16)])$ximp,
                      missForest(mdata_hiHTE[[4]][,-c(11:16)])$ximp,
                      missForest(mdata_noHTE[[5]][,-c(11:16)])$ximp)

#4. Apply missing imputation methods ---
ccest  <- allest(homo = mdata_noHTE, hete = mdata_hiHTE, functionv = getest, Tformv = 1, CCv = TRUE, typev = "CC" , approachv=NULL) # Complete case
minest <- allest(homo = mdata_noHTE, hete = mdata_hiHTE, functionv = getest, Tformv = 2, CCv = FALSE, typev = "MInd", approachv=NULL) # Missing indicator
micestny<-allest(homo = miceout_noHTEny, hete = miceout_hiHTEny, functionv = getmicest, Tformv = 1, CCv = FALSE, typev = "MICE.ny", approachv="within") #MICE no outcome within 
micesty<- allest(homo = miceout_noHTEy,  hete = miceout_hiHTEy,  functionv = getmicest, Tformv = 1, CCv = FALSE, typev = "MICE.y",   approachv="within") #MICE with outcome within  
micestiy<-allest(homo = miceout_noHTEiy, hete = miceout_hiHTEiy, functionv = getmicest, Tformv = 1, CCv = FALSE, typev = "MICE.inty",approachv="within") #MICE with outcome and interaction within  
micestyc<-allest(homo = miceout_noHTEy,  hete = miceout_hiHTEy,  functionv = getmicest, Tformv = 3, CCv = FALSE, typev = "MICE.y.mind",approachv="within") #MICE with outcome and MI within 
micests<- allest(homo = miceout_noHTEs,  hete = miceout_hiHTEs,  functionv = getmicest, Tformv = 1, CCv = FALSE, typev = "MICE.sep",approachv="within") #MICE with outcome and MI within 
micestiya<-allest(homo = miceout_noHTEiy, hete = miceout_hiHTEiy, functionv = getmicest, Tformv = 1, CCv = FALSE, typev = "MICE.inty.ac",approachv="across") #MICE with outcome and interaction within  
micestiyt<-allest(homo = miceout_noHTEiyt, hete = miceout_hiHTEiyt, functionv = getmicest, Tformv = 1, CCv = FALSE, typev = "MICE.intyt",approachv="within") #MICE with outcome and interaction within  
micestiyat<-allest(homo = miceout_noHTEiyt, hete = miceout_hiHTEiyt, functionv = getmicest, Tformv = 1, CCv = FALSE, typev = "MICE.intyt.ac",approachv="across") #MICE with outcome and interaction within  
micestyrf<- allest(homo = miceout_noHTEyrf,  hete = miceout_hiHTEyrf,  functionv = getmicest, Tformv = 1, CCv = FALSE, typev = "MICE.yrf",   approachv="within") #MICE with random forest  
micestycart<-allest(homo = miceout_noHTEycart,  hete = miceout_hiHTEycart,  functionv = getmicest, Tformv = 1, CCv = FALSE, typev = "MICE.ycart",   approachv="within") #MICE with cart
missF  <- allest(homo = for_mdata_noHTE, hete = for_mdata_hiHTE, functionv = getest, Tformv = 1, CCv = FALSE, typev = "MissForest" , approachv=NULL) # Random Forest

save.image(file = "/Users/jmunozav/Desktop/Chapter_book/Book_Missing_Causality15.RData")  
```

## Presence of HTE


```{r}
#1. Generate full datasets according to heterogeneity in treatment specification ----
#1. Generate full datasets according to heterogeneity in treatment effect specification ----
data_noHTE <- generate_data(n = 3000, beta = c(-0.2, -0.2, -0.2, -0.2, -0.2), seed = 1234) # no treatment effect
data_hiHTE <- generate_data(n = 3000, beta = c(log(0.4), log(0.5), log(1), log(1.1), log(1.2)), seed = 1234) # strong treatment effect
```


#load(file = "/Users/jmunozav/Desktop/Chapter_book/Book_Missing_Causality5.RData")

# Results under different scenarios
```{r}
#7 Make the plot ----
alldata<-rbind(true,ccest,minest,micestny,micesty,micestiy,micestiyt,micestyc,micests,micestiya,micestiyat,missF,micestyrf,micestycart)
alldata<-as.data.table(alldata)
alldata[,mean:=exp(estimate)]
alldata[,LCI:=exp(`2.5 %`)]
alldata[,UCI:=exp(`97.5 %`)]
alldata[,Scenario:=as.factor(Scenario)]
alldata[,Scenario:=factor(Scenario,levels=c("Full","MCAR","MAR","MART","MART2","MNAR"))]  
alldata[,T_effect:=as.factor(Treatment)] 
alldata[,T_effect:=factor(T_effect,levels=c("Homogeneus","Heterogeneus"))]
alldata[,type:=as.factor(type)]
alldata[,type:=factor(type,levels=c("True","CC","MInd","MICE.ny","MICE.y","MICE.y.mind","MICE.sep","MICE.intyt","MICE.intyt.ac","MICE.yrf","MICE.ycart","MissForest"))]
```

##Complete case analysis (CC)
```{r}
theme_set(theme_bw())
pd = position_dodge(.4)
ggplot(alldata[type%in%c("CC")],aes(x= estimand, y= mean,color =Scenario)) +
  geom_point(shape = 1,
             size  = 1,
             position = pd) +
  geom_errorbar(aes(ymin  = LCI,
                    ymax  = UCI),
                width = 0.2,
                size  = 0.5,
                position = pd)+
  geom_hline(data = alldata[type=="True"], aes(yintercept =mean))+
  facet_grid(method ~ T_effect+estimand,scales="free")+ see::scale_color_flat()+theme_light()       
```

## Missing indicator (MI)
```{r}
ggplot(alldata[type%in%c("CC","MInd","MICE.ny","MICE.y","MICE.y.mind")&Scenario!="MART"],aes(x= Scenario, y= mean,color =type)) +
  geom_point(shape = 1,
             size  = 1,
             position = pd) +
  geom_errorbar(aes(ymin  = LCI,
                    ymax  = UCI),
                width = 0.2,
                size  = 0.5,
                position = pd)+
  geom_hline(data = alldata[type=="True"&Scenario!="MART"], aes(yintercept =mean))+
  facet_grid(T_effect+method ~ estimand,scales="free")+ see::scale_color_flat()+theme_light()       

```

## Multiple Imputation (MICE)
```{r}
ggplot(alldata[type%in%c("MICE.ny","MICE.y","MICE.intyt","MICE.intyt.ac","MICE.sep","MICE.y.mind")&Scenario!="MART"],aes(x= Scenario, y= mean,color =type)) +
  geom_point(shape = 1,
             size  = 1,
             position = pd) +
  geom_errorbar(aes(ymin  = LCI,
                    ymax  = UCI),
                width = 0.2,
                size  = 0.5,
                position = pd)+
  geom_hline(data = alldata[type=="True"], aes(yintercept =mean))+
  facet_grid(T_effect+method ~ estimand,scales="free")+ see::scale_color_flat()+theme_light()       
```

## Machine learning imputation methods
```{r}
ggplot(alldata[type%in%c("MICE.intyt","MICE.sep","MissForest","MICE.yrf","MICE.ycart")&Scenario!="MART"],aes(x= Scenario, y= mean,color =type)) +
  geom_point(shape = 1,
             size  = 1,
             position = pd) +
  geom_errorbar(aes(ymin  = LCI,
                    ymax  = UCI),
                width = 0.2,
                size  = 0.5,
                position = pd)+
  geom_hline(data = alldata[type=="True"], aes(yintercept =mean))+
  facet_grid(T_effect+method ~ estimand,scales="free")+ see::scale_color_flat()+theme_light()       
```


