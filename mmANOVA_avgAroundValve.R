library(ez)

setwd("C:/Users/david/OneDrive - Indiana University/localAlcohol/ephysData")
data<-read.delim("averageAroundValve.csv", header = TRUE, sep=",")


#controlVsAlcohol<-c(-2, 1, 1)
#contrasts(data$treatment)<-cbind(controlVsAlcohol)

options(digits = 3)
model<-ezANOVA(
  data = data,
  dv = .(spk),
  wid = .(clusterID),
  between = .(treatment),
  within = .(bin),
  type = 3,
  detailed = FALSE
  )
model

