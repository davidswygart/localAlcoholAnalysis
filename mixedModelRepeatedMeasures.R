#---------------------------------------------------------------------------------------------------------
#R Code for Chapter 14 of:
#
#Field, A. P., Miles, J. N. V., & Field, Z. C. (2012). Discovering Statistics Using R: and Sex and Drugs and Rock 'N' Roll. #London Sage
#
#(c) 2011 Andy P. Field, Jeremy N. V. Miles & Zoe C. Field
#-----------------------------------------------------------------------------------------------------------

#Set the working directory (you will need to edit this to be the directory where you have stored the data files for this Chapter)

setwd("C:/Users/david/Desktop/Ch14_mixed")
imageDirectory<-"C:/Users/david/Desktop/Ch14_mixed"

setwd("C:/Users/david/Desktop/Ch14_mixed")
imageDirectory<-"C:/Users/david/Desktop/Ch14_mixed"

######Install packages

install.packages("compute.es")
install.packages("ez")
install.packages("ggplot2")
install.packages("multcomp")
install.packages("nlme")
install.packages("pastecs")
install.packages("reshape")
install.packages("WRS", repos="http://R-Forge.R-project.org")


#Initiate packages
library(compute.es)
library(ez)
library(ggplot2)
library(multcomp)
library(nlme)
library(pastecs)
library(reshape)
library(WRS)
source("http://www-rcf.usc.edu/~rwilcox/Rallfun-v14")


#--------Speed Dating Example ----------


data<-read.delim("alcoholInject.csv", header = TRUE, sep=",")

#using ezAnova
#SomevsNone<-c(1, 1, -2)
#HivsAv<-c(1, -1, 0)

#AttractivevsUgly<-c(1, 1, -2)
#AttractvsAv<-c(1, -1, 0)

#contrasts(speedData$personality)<-cbind(SomevsNone, HivsAv)
#contrasts(speedData$looks)<-cbind(AttractivevsUgly, AttractvsAv)

options(digits = 3)
model<-ezANOVA(data = data, dv = .(spk_z), wid = .(clusterID),  between = .(treatment), within = .(bin), type = 3, detailed = TRUE)
model
options(digits = 7)

#Using lme


HighvsAv<-c(1, 0, 0)
DullvsAv<-c(0, 0, 1)

AttractivevsAv<-c(1, 0, 0)
UglyvsAv<-c(0, 0, 1)

contrasts(speedData$personality)<-cbind(HighvsAv, DullvsAv)
contrasts(speedData$looks)<-cbind(AttractivevsAv, UglyvsAv)
contrasts(speedData$gender)<-c(1, 0)

#Building the model
baseline<-lme(dateRating ~ 1, random = ~1|participant/looks/personality, data = speedData, method = "ML")
looksM<-update(baseline, .~. + looks)
personalityM<-update(looksM, .~. + personality)
genderM<-update(personalityM, .~. + gender)
looks_gender<-update(genderM, .~. + looks:gender)
personality_gender<-update(looks_gender, .~. + personality:gender)
looks_personality<-update(personality_gender, .~. + looks:personality)
speedDateModel<-update(looks_personality, .~. + looks:personality:gender)


anova(baseline, looksM, personalityM, genderM, looks_gender, personality_gender, looks_personality, speedDateModel)
summary(speedDateModel)



#-------Effect Sizes

rcontrast<-function(t, df)
{r<-sqrt(t^2/(t^2 + df))
print(paste("r = ", r))
}



t<-summary(speedDateModel)$tTable[,4]
df<-summary(speedDateModel)$tTable[,3]

rcontrast(-1.20802, 108)
rcontrast(3.85315, 108)
rcontrast(-7.53968, 108)
rcontrast(-0.97891, 108)