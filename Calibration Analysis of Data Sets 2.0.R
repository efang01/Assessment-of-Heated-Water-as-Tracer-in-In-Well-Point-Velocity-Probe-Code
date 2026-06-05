################################################################################
#
#       Calibration curves for IWPVP
#       Written by J.F. Devlin
#       Modified by Ellis Fangmann
#       April 16, 2026
#
#       Fix: threshold lines (LoQ/MDL/ProbeMDL) now intersect the horizontal
#       ProbeMDL line exactly on the confidence curves. Crossings are solved
#       directly with approx() at the unrounded b_up height; rounding is applied
#       only to the printed labels, never to the plotted geometry.
#
################################################################################

# Load and prepare data
# fileEncoding handles a possible UTF-8 BOM on the header (otherwise the first
# column name comes in mangled and lm() fails)
inputpvp<-read.csv("4ChannelData3ExpDAY03.csv", stringsAsFactors = F,
                   fileEncoding = "UTF-8-BOM")
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
title(xlab = expression("Outside Well (cm/day)"), line=2.7, cex.lab=2, font.lab=1)
title(ylab = expression("In Probe (cm/day)"), line=1.7, cex.lab=2, font.lab=1)

# Axis major ticks
axis(1, at = pretty(x, n = 6))
axis(2, at = pretty(y, n = 6))
# Optional axis minor ticks
#axis(1, at = pretty(x, n = 20), tcl = -0.25, labels = FALSE)
#axis(2, at = pretty(y, n = 20), tcl = -0.25, labels = FALSE)

abline(lm(y~x), col="black", lty=3, lwd=2)
lines(xconfint, conf_interval[,2], col="#00CD00", lty=2, lwd=2)
lines(xconfint, conf_interval[,3], col="#00CD00", lty=2, lwd=2)

# R-squared
rsq <- (summary(regcurrent)$r.squared)

# Safe plotting positions
xpos <- min(x, na.rm=TRUE) + 0.05*diff(range(x, na.rm=TRUE))
ypos <- max(y, na.rm=TRUE) - 0.05*diff(range(y, na.rm=TRUE))

# Equation label (handles + / - automatically)
eqn <- paste0("y = ", m, "x ", ifelse(b >= 0, "+ ", "- "), abs(b))
text(xpos, ypos, labels = eqn, pos=4, font=4, cex=1.5)

# R²
text(xpos, ypos - 0.05*diff(range(y, na.rm=TRUE)),
     labels = paste("R² =", format(rsq, digits=12, scientific = FALSE)),
     pos=4, font=4, cex=1.5)

ci<-as.data.frame(conf_interval)
ci$x<-xconfint
ci<-ci[c(4,1,2,3)]          # columns now ordered: x, fit, lwr, upr

################################################################################

b_up <- min(ci$upr)                            # unrounded -> used for geometry
tc  <- approx(x = ci$fit, y = ci$x, xout = b_up)$y   # LoQ
mdl <- approx(x = ci$lwr, y = ci$x, xout = b_up)$y   # MDL

# Rounded values for labels ONLY
tc_lab   <- signif(tc, 1)
mdl_lab  <- round(mdl, 6)
b_up_lab <- round(b_up, 6)

tc
mdl
b_up

# Plot calibration thresholds
plot(x,y,
     xlab="", ylab="", xaxt="n", yaxt="n",
     pch=16,xlim=c(0,2*mdl),ylim=c(0,2*b_up))
title(xlab = expression("Outside Well (cm/day)"), line=2.7, cex.lab=2, font.lab=1)
title(ylab = expression("In Probe (cm/day)"), line=1.8, cex.lab=2, font.lab=1)

# Axis major ticks
axis(1, at = pretty(c(0,2*mdl), n = 6))
axis(2, at = pretty(c(0,2*b_up), n = 6))
# Optional axis minor ticks
#axis(1, at = pretty(c(0,2*mdl), n = 20), tcl = -0.25, labels = FALSE)
#axis(2, at = pretty(c(0,2*b_up), n = 20), tcl = -0.25, labels = FALSE)

# Build multi-colored title (labels use rounded values)
mtext("Calibration Thresholds", side=3, line=1, col="black", font=2, cex=2)
mtext(paste("LoQ =", tc_lab), side=3, line=0, col="red", font=4, cex=1.3, adj=0)
mtext(paste("MDL =", mdl_lab), side=3, line=0, col="blue", font=4, cex=1.3 ,adj=0.5)
mtext(paste("ProbeMDL =", b_up_lab), side=3, line=0, col="lightblue", font=4, cex=1.3, adj=1)

#Plot upper and lower bounds
lines(xconfint, ci[,2], col="black", lty=2, lwd=2)     # fit
lines(xconfint, ci[,3], col="#00CD00", lty=2, lwd=2)   # lwr
lines(xconfint, ci[,4], col="#00CD00", lty=2, lwd=2)   # upr

# Plot probe mdl, loq, mdl -- drawn at the SAME unrounded values used above
abline(h=b_up, col="lightblue", lty=3, lwd=2)
abline(v=tc, col="red", lty=1, lwd=2)
abline(v=mdl, col="blue", lty=1, lwd=2)

# Print coefficients
summary(regcurrent)