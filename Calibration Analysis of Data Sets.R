################################################################################
#
#       Calibration curves for IWPVP
#       Written by J.F. Devlin 
#       Modified by Ellis Fangmann
#       April 16, 2026
#
################################################################################

# Load and prepare data
inputpvp<-read.csv("4ChannelData2ExpDAY.csv", stringsAsFactors = F)
x<-(inputpvp$SandFlux)
y<-(inputpvp$ProbeFlux)

regcurrent<-lm(formula=y~x)

# Extract slope and intercept
b <- round(regcurrent$coefficients[1], 4)  # intercept
m <- round(regcurrent$coefficients[2], 4)  # slope

# Vector for calculating confidence intervals
xconfint<-seq(0,70,by=0.0001)
conf_interval <- predict(regcurrent, newdata=(data.frame(x=xconfint)), 
                         interval="confidence", level = 0.95)
par(mfrow=c(1,2))

# Plot Calibration line
plot(x,y,xlab="",ylab="", xaxt="n", yaxt="n",
     pch=16,cex.main=2,main=("Calibration Line"))
title(xlab="Outside Well", line=2.5, cex.lab=2, font.lab=1)   # move up/down
title(ylab="In Probe", line=2.3, cex.lab=2, font.lab=1)

# Axis major ticks
axis(1, at = pretty(x, n = 8))
axis(2, at = pretty(y, n = 8))
# Optional axis minor ticks
#axis(1, at = pretty(x, n = 20), tcl = -0.25, labels = FALSE)
#axis(2, at = pretty(y, n = 20), tcl = -0.25, labels = FALSE)

abline(lm(y~x), col="black", lty=3, lwd=2)
lines(xconfint, conf_interval[,2], col="#00CD00", lty=2, lwd=4)
lines(xconfint, conf_interval[,3], col="black", lty=2, lwd=2)

# R-squared
rsq <- summary(regcurrent)$r.squared
text(min(x), max(y)*0.95, pos=4, font=4, paste("R² =",(rsq)))

# Equation label (handles + / - automatically)
eqn <- paste0("y = ", m, "x ", ifelse(b >= 0, "+ ", "- "), abs(b))
text(min(x), max(y)*0.98, pos=4, font=4, eqn)

ci<-as.data.frame(conf_interval)
ci$x<-xconfint
ci<-ci[c(4,1,2,3)]

# Determine MDL
b_up<-min(ci$upr)

# Find MDL
i=11
while(ci[i,3]<b_up){
  yy<-ci[i+1,3]
  mdl<<-round(ci[i+1,1],4)
  i=i+1
}
mdl

# Find limit of quantification (threshold value)
i=1
while(ci[i,2]<b_up){
  yy<-ci[i+1,2]
  tc<<-round(ci[i+1,1],4)
  i=i+1
}
tc

# Find probe MDL
b_up<-round(min(ci$upr),4)

# Plot calibration thresholds
plot(x,y,
     xlab="", ylab="", xaxt="n", yaxt="n",
     pch=16,xlim=c(0,2*mdl),ylim=c(0,2*b_up))
title(xlab="Outside Well", line=2.5, cex.lab=2, font.lab=1)   # move up/down
title(ylab="In Probe", line=2.3, cex.lab=2, font.lab=1)

# Axis major ticks
axis(1, at = pretty(c(0,2*mdl), n = 6))
axis(2, at = pretty(c(0,2*b_up), n = 6))
# Optional axis minor ticks
#axis(1, at = pretty(c(0,2*mdl), n = 20), tcl = -0.25, labels = FALSE)
#axis(2, at = pretty(c(0,2*b_up), n = 20), tcl = -0.25, labels = FALSE)

# Build multi-colored title
mtext("Calibration Thresholds", side=3, line=1, col="black", font=2, cex=2)
mtext(paste("LoQ =", tc), side=3, line=0, col="red", font=4, cex=1.3, adj=0)
mtext(paste("MDL =", mdl), side=3, line=0, col="blue", font=4, cex=1.3 ,adj=0.5)
mtext(paste("ProbeMDL =", b_up), side=3, line=0, col="black", font=4, cex=1.3, adj=1)

#Plot upper and lower bounds
lines(xconfint, ci[,2], col="black", lty=2, lwd=2)
lines(xconfint, ci[,3], col="#00CD00", lty=2, lwd=2)
lines(xconfint, ci[,4], col="#00CD00", lty=2, lwd=2)

# Plot probe mdl, loq, mdl
abline(h=b_up, col="black", lty=1, lwd=2)
abline(v=tc, col="red", lty=1, lwd=2)
abline(v=mdl, col="blue", lty=1, lwd=2)

# Print coefficients
summary(regcurrent)