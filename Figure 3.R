########RFM##################
install.packages("C:/Users/shirl/Documents/R/win-library/4.1/runjags_2.2.2-4.zip",repos = NULL,type = "win.binary")
install.packages("C:/Users/shirl/Documents/R/win-library/4.1/swfscMisc_1.6.5.tar.gz",repos = NULL,type = "source")
install.packages("C:/Users/shirl/Documents/R/win-library/4.1/rfPermute_2.5.2.tar.gz",repos = NULL,type = "source")
install.packages("C:/Users/shirl/Documents/R/win-library/4.1/rfUtilities_2.1-5.tar.gz",repos = NULL,type = "source")

library(rfPermute)
library(randomForest)
library(rfUtilities)

setwd("C:/graduate/WH赤霉病/扩增子")
rm(list=ls())


fun_env<-read.csv("all.RF3.csv",head=T,row.names=1)

#new-all
fun.rf <- randomForest(fun_env$Symptom~AVD+stability+shannon+dispersion+AVD_fun+stability_fun+shannon_fun+dispersion_fun+RA_Fusarium,
                       data=fun_env,num.cores = 4, nrep = 1000)##ntree=500 mtry = 6 default
rf.significance(fun.rf,fun_env[,c(4:11)])##change the combined table with different columns,?????????? the same result P=0.001,R2=0.399
fun.rf##the details of model

fun.rfP1 <- rfPermute(fun_env$Symptom~AVD+stability+shannon+dispersion+AVD_fun+stability_fun+shannon_fun+dispersion_fun+RA_Fusarium,
                      data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP1)
plotImportance(fun.rfP1,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)


##############SEM#############
library(vegan)
library(nlme)
library(labdsv)
require(ecodist)
library(devtools)
library(piecewiseSEM)
library(multcompView)
library(lavaan)
library(semPlot)

dat <- read.csv("结构方程all.csv")

######bacteria#######
sem_1 <- psem(
  lm(State ~ Bacterial_robustness+Bacterial_dispersion+RA_Fusarium, data = dat),
  lm(Bacterial_robustness ~ log(Bacterial_diversity)+Bacterial_dispersion+Bacterial_AVD+RA_Fusarium, data = dat),
  lm(Bacterial_AVD ~ log(Bacterial_diversity)+Bacterial_dispersion+RA_Fusarium, data = dat)
)


# Get the basis set
basisSet(sem_1)
# Conduct d-sep tests
dTable_sem_1 <- dSep(sem_1)
# Get Fisher’s C
fisherC(dTable_sem_1)
# Get coefficients
coefs(sem_1)
# Get R-squared
rsquared(sem_1)
# summary
summary(sem_1)


#######fungi############
sem_1 <- psem(
  lm(State ~ Fungal_robustness+Fungal_dispersion+RA_Fusarium, data = dat),
  lm(Fungal_robustness ~Fungal_dispersion+Fungal_diversity+RA_Fusarium, data = dat)
)


# Get the basis set
basisSet(sem_1)
# Conduct d-sep tests
dTable_sem_1 <- dSep(sem_1)
# Get Fisher’s C
fisherC(dTable_sem_1)
# Get coefficients
coefs(sem_1)
# Get R-squared
rsquared(sem_1)
# summary
summary(sem_1)