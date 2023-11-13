require(OpenMx)   #Loads OpenMx
require(psych)   #Loads Psych package

library(devtools)
library(umx)
library (OpenMx)
library (haven)
library (dplyr)
library(car)
library(psych)
library(lme4)
library(nlme)
library(lmerTest)
library(ggplot2)
library(tidyverse)
library(lme4)
library(car)
library(multilevelTools)
source("GenEpiHelperFunctions.R")
source("polychoricMeansMatrix3.R")
library(haven)
#install.packages("JWileyMisc")
library(JWileymisc)
library(effects)
#install.packages("sjPlot")
library(sjPlot)
#install.packages("sjmisc")
library(sjmisc)
library(sjPlot)
library(lme4)
library(ggplot2)
library("bestNormalize")

mxOption(NULL, "Default optimizer", "SLSQP")  

#READ IN DATA; 
setwd("Z:/Tyler/Tyler Master VETSA File/")

VETSA <- read_sav("Tyler Bell_VETSA_Master_Data_05232023.sav")

# PULL FROM MASTER DATA. 
names(VETSA)
names(VETSA)<- toupper(names(VETSA))

twins<-VETSA
#twins<-subset(SMdata, SMdata$RMCI_CONS_V1 != 1 | SMdata$RMCI_CONS_V1 !=2 | SMdata$RMCI_CONS_V1  != 3 | SMdata$RMCI_CONS_V1  != 4)
#load("C:/Users/tyler/Downloads/Tyler08092022.RData")
#names(newtwins
#View(twins)
names(twins)
names(twins)<- toupper(names(twins))
View(twins)


## Creating MZ and DZ data sets ##
twinA <- twins[twins$TWIN=="1",]
twinB <- twins[twins$TWIN=="2",]

df2 <- merge(twinA, twinB, by=c("CASE","ZYG2019"),all.x=TRUE, all.y=TRUE,suffixes=c("_T1","_T2"))
names(df2)

#save.image(file = "SM GENETIC PAPER MAIN MODELS.RData")
#load("/Volumes/ngillespie/Documents/work/projects/2016/2016. 8. VETSA/analyses/1. Tyler/subjective memory/SM GENETIC PAPER MAIN MODELS.RData")

df2<-newtwins


############################################################
#POLYCHOR CORRELATIONS ACCOUNTING FOR TWINNESS AND ZYGOSITY
############################################################

#POLYCHOR WITH CONTINUOUS VARIABLES

selVarsMZ	= c("ZLOG_RESIDUALS_AB42_MEAN_T1","CNR_ROSTRAL2_T1",
              "ZLOG_RESIDUALS_AB42_MEAN_T2")
selVarsDZ	= c("ZLOG_RESIDUALS_AB42_MEAN_T1","CNR_ROSTRAL2_T1",
              "ZLOG_RESIDUALS_AB42_MEAN_T2","CNR_ROSTRAL2_T2")


mzdata 	= subset(df2, ZYG2019==1, selVarsMZ)
dzdata 	= subset(df2, ZYG2019==2, selVarsDZ)

describe(mzdata)
describe(dzdata)

nv     	= length(selVarsMZ)/2     
ntv    	= length(selVarsMZ)


labelsMZ =c(
  "a", 	
  "c", "b", 	
  "d", "c", "a")

labelsDZ =c(
  "a", 	
  "c", "b", 	
  "d", "f", "a", 
  "f", "e", "c", "b")


freeMZ =
  c(
    T,
    T,T,	
    T,T,T)

freeDZ =
  c(
    T,
    T,T,	
    T,T,T,	
    T,T,T,T)

valsMZ =
  c(
    5,
    0.5,5,	
    0.5,0.5,5)

valsDZ =
  c(
    5,
    0.5,5,	
    0.5,0.5,5,	
    0.5,0.5,0.5,5)

temp = mxModel("corr",
               mxModel("top",
                       mxMatrix( name="MeansMZ", 	type="Full", nrow=1, ncol=3, free=T, labels=c("u1","u2","u1")),
                       mxMatrix( name="MeansDZ", 	type="Full", nrow=1, ncol=4, free=T, labels=c("u1","u2","u1","u2")),  
                       mxMatrix( name="CovMZ",		type="Symm", nrow=3, ncol=3, free=freeMZ, labels=labelsMZ, values=valsMZ,byrow=T, lbound=-20,ubound=50),
                       mxMatrix( name="CovDZ",		type="Symm", nrow=4, ncol=4, free=freeDZ, labels=labelsDZ, values=valsDZ,byrow=T, lbound=-20,ubound=50),
                       mxAlgebra(name="corrMZ",	expression=cov2cor(CovMZ)),
                       mxAlgebra(name="corrDZ",	expression=cov2cor(CovDZ)),
                       mxCI(c("corrDZ[4,3]")) ),    
               mxModel("MZ", mxData(		mzdata, type ="raw"),mxExpectationNormal("top.CovMZ", means="top.MeansMZ", dimnames=selVarsMZ), mxFitFunctionML() ),
               mxModel("DZ", mxData(		dzdata, type ="raw"),mxExpectationNormal("top.CovDZ", means="top.MeansDZ", dimnames=selVarsDZ), mxFitFunctionML() ),
               mxFitFunctionMultigroup(c("MZ","DZ")) )
omxGetParameters(temp)

install.OpenMx("NPSOL")


model_fit = mxTryHard(temp)   
model_fit = mxRun(model_fit,intervals=T)     
summary( model_fit)

model_fit$top$CovMZ
model_fit$top$CovDZ  
model_fit$top$corrMZ
model_fit$top$corrDZ

model_fit = paste0(
  format(round(print(summary(SM_38_fit)$CI)[1,2],2),nsmall=2)," (",
  format(round(print(summary(SM_38_fit)$CI)[1,1],2),nsmall=2),", ",
  format(round(print(summary(SM_38_fit)$CI)[1,3],2),nsmall=2),")")
model_fit
