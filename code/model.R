# Author: Jingwei Zhang
# Date: July 28，2025
# Purpose: Statistical Model for the Development and 
#          Validation of a Stroke Heat Risk Grading Prediction Model

# NOTE: this code is a guide for transparency and 
#       reproducibility and is not able to be run

# Key packages used
library(dlnm)       # For distributed lag nonlinear models (DLNM)
library(splines)    # For spline functions
library(mvmeta)     # For multivariate meta-analysis
library(survival)   # For survival analysis (used in case-crossover)
library(tidyverse)  # For data transformation and structuring

# 1.Construction of the Stroke Heat Risk Grading Prediction Model  
# Stage 1. Time-series analysis using DLNM across counties  
# NOTE: time-series analyses were conducted using quasi-Poisson regression combined with DLNM, 
#       separately within age- and gender-specific subgroups, based on stroke mortality, 
#       meteorological, and air pollution data during the summer months from 2013 to 2018.
#      
# Construct a list of data frames for each selected county
counties <- as.character(unique(
  #mcode is the unique identifier for each county
  #data1 is county-level time-series data on stroke mortality from 2013 to 2018
  data1$mcode))
datalist <- lapply(counties, function(x) data1[data1$mcode==counties,])
names(datalist) <- counties

#Calculate temperature ranges
ranges <- t(sapply(datalist,
                   #Tmeanc is the daily mean temperature
                   function(x) range(x$Tmeanc,na.rm=T)))
bound <- colMeans(ranges)

# Model temperature using a natural cubic spline with three internal knots
perc <- c(10,75,90)
varknots <- rowMeans(sapply(datalist,function(x) {
  quantile(x$Tmeanc,perc/100,na.rm=T)
}))
argvar <- list(fun="ns", degree=3, knots=quantile(data1$Tmeanc, perc/100,na.rm=T))

# Define lag structure on log scale
arglag <- list(knots=logknots(lag=3,lagnk=2))

# Create crossbasis function for temperature
cb <- crossbasis(data1$Tmeanc, lag=3, argvar=argvar, arglag=arglag)

# Fit the quasi-Poisson model with adjustment variables
# NOTE: Fit the quasi-Poisson model with adjustment variables 
#       for one selected city as an example
model <- glm(
  # indicating the county-level daily stroke mortality 
  str ~ 
    cb + 
    
    # adjust for daily mean relative humidity with 3 df
    ns(Rhu_mean, df=3) + 
    
    # adjust for daily mean wind speed with 3 df
    ns(Wind_mean, df=3) +   
    
    # adjust for daily average concentrations of PM2.5 and O3
    PM2.5 + M8H_O3 +   
    
    # adjust for day of week
    as.factor(dow) +     
    
    # adjust for long-term trend (7 df per year for 6 years)(2013-2018)
    ns(time, df=7*6),      
  #data1 is county-level time-series data on stroke mortality from 2013 to 2018
  data=data1,
  family = quasipoisson(link = "log"),
  control = glm.control(epsilon = 10E-8, maxit = 5000))

# Stage 2: Random-effects meta-analysis based on all the study counties
# NOTE: meta-analyses were conducted separately within 
#       age- and gender-specific subgroups.

# Perform random-effects meta-analysis
mv <- mvmeta(
  #ymat and Slist represent the estimated 
  #coefficients and variance–covariance matrices, respectively
  ymat, Slist, 
  # use "Maximum likelihood estimation"
  method="ml"
)
summary(mv)
#Stage3: Derive temperature–mortality associations and corresponding risk estimates
# Create a sequence of temperature
Tmeanc <- seq(bound[1],bound[2],length=30)

#Construct the basis matrix for the temperature variable
bTmeanc <- onebasis(Tmeanc,fun="ns",degree=3,knots=varknots,Bound=bound)


#Generate temperature percentiles
percTmeanc<-rowMeans(sapply(datalist,function(x){quantile(x$Tmeanc,c(0:1/2,1:100)/100,na.rm=T)}))
perTmeancdata<-data.frame(percent=names(percTmeanc),percTmeanc=round(percTmeanc,1))

#Prediction of temperature–mortality associations
cprel<-crosspred(bTmeanc,coef=coef(mv),vcov=vcov(mv),at=c(percTmeanc,varknots),model.link="log")

#Find MMT (minimum mortality temperature)
mmt=round(as.numeric(names(which.min(cprel$allRRfit))),1)

#Reprediction of temperature–mortality associations using the MMT as the reference point
cen <- perTmeancdata[perTmeancdata$percTmeanc==mmt,]
cprel<-crosspred(bTmeanc,coef=coef(mv),vcov=vcov(mv),cen=mmt,at=c(percTmeanc,varknots),model.link="log")

#Obtain the relative risk (RR) value corresponding to each temperature
RR<-data.frame(Tem=list(cprel$predvar),B=cprel$allfit,SE=cprel$allse,RR=cprel$allRRfit,low=cprel$allRRlow,high=cprel$allRRhigh)##获取各参数
colnames(RR)<-c("Tem","b","se","RR","low","High")

#Obtain the relative risk (RR) value corresponding to each temperature percentile
pt<-data.frame(percent=names(percTmeanc),Tem=percTmeanc)
RR_perc<-left_join(pt,RR,by=c("Tem"))

# 2.Validation of the Stroke Heat Risk Grading Prediction Model
# NOTE: As a validation analysis, case-crossover analyses were conducted based on stroke mortality, 
#       meteorological, and air pollution data during the summer months from 2019 to 2022.

# Fit conditional logistic regression
model <- clogit(
  case ~ 
    #an ordinal variable with “1” for identified "moderate risk" level, 
    #“2” for identified "high risk" level,
    #“3” for identified "extremely high risk" level,
    #and “0” for identified "low risk" level
    risk_level + 
    
    # adjust for daily mean relative humidity with 3 df
    ns(Rhu_mean, df = 3) + 
    
    # adjust for daily mean wind speed with 3 df
    ns(Wind_mean, df = 3) + 
    
    # adjust for daily average concentrations of PM2.5 and O3
    PM2.5 + M8H_O3 + 
    
    # Match on case-control sets
    strata(ID),  
  #data2 is individual-level stroke mortality data from 2019 to 2022
  data = data2,
  #breslow method is a maximum likelihood approach for matched case-control studies 
  method = "breslow"
)

# Calculate OR for each risk level
or_results <- summary(model)$coefficients %>%
  as.data.frame() %>%
  mutate(
    OR = exp(Estimate),
    OR_low = exp(Estimate - 1.96 * summary(model)$coefficients[,"Std. Error"]),
    OR_high = exp(Estimate + 1.96 * summary(model)$coefficients[,"Std. Error"])
  )
