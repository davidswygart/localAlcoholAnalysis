library(ez)

setwd("C:/Users/david/OneDrive - Indiana University/localAlcohol/ephysData")
data<-read.delim("injectFiring.csv", header = TRUE, sep=",")

#using ezAnova
#SomevsNone<-c(1, 1, -2)
#HivsAv<-c(1, -1, 0)

#contrasts(speedData$personality)<-cbind(SomevsNone, HivsAv)
#contrasts(speedData$looks)<-cbind(AttractivevsUgly, AttractvsAv)

options(digits = 3)
model<-ezANOVA(
  data = data,
  wid = .(clusterID), # within group ID
  within = .(bin), # within group predictor variable
  between = .(treatment), # between group predictor variable
  dv = .(spk), # dependent variable
  type = 3,
  detailed = FALSE
  )
model


#ezDesign(
#  data = data,
#  x = .(treatment),
#  y = .(recordingID),
#)
