library(ez)

setwd("C:/Users/david/OneDrive - Indiana University/localAlcohol/ephysData")
data<-read.delim("effectOfInject.csv", header = TRUE, sep=",")

#using ezAnova
#SomevsNone<-c(1, 1, -2)
#HivsAv<-c(1, -1, 0)

#contrasts(speedData$personality)<-cbind(SomevsNone, HivsAv)
#contrasts(speedData$looks)<-cbind(AttractivevsUgly, AttractvsAv)

options(digits = 3)
model<-ezANOVA(
  data = data,
  dv = .(spk),
  wid = .(clusterID),
  between = .(treatment),
  within = .(bin),
  type = 3,
  detailed = TRUE
  )
model
