################################################################################
#
#
#       Calibration curves for IWPVP
#       Written by J.F. Devlin 
#       Modified by Ellis Fangmann
#       Feb 22, 2026
#
#
################################################################################
inputpvp<-read.csv("FileName.csv", stringsAsFactors = F)
x<-(inputpvp$SandFlux)
y<-(inputpvp$ProbeFlux)


regcurrent<-lm(formula=y~x)
m<-round(regcurrent$coefficients[2],1)
#Vector for calculating confidence intervals
xconfint<-seq(0,70,by=0.1)

conf_interval <- predict(regcurrent, newdata=(data.frame(x=xconfint)), interval="confidence",
                         level = 0.95)

par(mfrow=c(1,2))
plot(x,y, xlab="outside well", ylab="in probe",
     pch=16,main=paste("calibration factor (slope) = ",m))
abline(lm(y~x),lty=3,lwd=2)
lines(xconfint, conf_interval[,2], col="black", lty=2)
lines(xconfint, conf_interval[,3], col="black", lty=2)
rsq <- summary(regcurrent)$r.squared
rsq
text(10, 1050, pos = 4, paste("R squared =", round(rsq, digits = 3) ) )

ci<-as.data.frame(conf_interval)
ci$x<-xconfint
ci<-ci[c(4,1,2,3)]

#Determine MDL
b_up<-min(ci$upr)
#find MDL
i=1
while(ci[i,3]<b_up){
  yy<-ci[i+1,3]
  mdl<<-round(ci[i+1,1],1)
  i=i+1
}
mdl

#Find limit of quantification (threshold value)
i=1
while(ci[i,2]<b_up){
  yy<-ci[i+1,2]
  tc<<-round(ci[i+1,1],1)
  i=i+1
}
tc
b_up<-round(min(ci$upr),1)

plot(x,y,
     #pch=16,xlim=c(0,2*mdl),ylim=c(0,2*b_up),
     pch=16,xlim=c(0,0.005),ylim=c(0,2*b_up),
     main=paste("LoQ = ",tc,
                " MDL = ",mdl,"\n"," probe mdl  =",b_up))
lines(xconfint, ci[,2], col="black", lty=2,
      xlim=c(0,2*mdl))
lines(xconfint, ci[,3], col="black", lty=2,
      xlim=c(0,2*mdl))
lines(xconfint, ci[,4], col="black", lty=2,
      xlim=c(0,2*mdl),ylim=c(0,2*b_up))
abline(h=b_up,col="lightgray")
abline(v=tc,col="pink")
abline(v=mdl,col="lightblue")

summary(regcurrent)