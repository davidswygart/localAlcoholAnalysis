library(ez)

setwd("C:/Users/david/OneDrive - Indiana University/localAlcohol/ephysData")
#controlVsAlcohol<-c(-2, 1, 1)
#contrasts(data$treatment)<-cbind(controlVsAlcohol)


data<-read.delim("ExcPeak.csv", header = TRUE, sep=",")
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

data<-read.delim("InhPeak.csv", header = TRUE, sep=",")
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

