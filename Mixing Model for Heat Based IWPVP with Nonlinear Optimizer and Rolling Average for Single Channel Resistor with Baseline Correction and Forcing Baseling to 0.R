################################################################################
# 
#
#       Mixing Model for Heat Based IWPVP w/ Nonlinear Fit + Rolling Average + Baseline Correction for Single Channel Resistor + Forcing to 0
#       Written by Ellis Fangmann
#       Feb 12, 2026
#
#
################################################################################

# DONT FORGET TO SELECT CHANNEL!

# Load required libraries
library(dplyr)
library(lubridate)
library(stringr)
library(zoo)

# Load and prepare data
ExpData <- read.table("CR1000_TempSensAug25Test6.dat", header = TRUE, sep = ",", skip = 1)

ExpData <- ExpData %>%
  filter(str_detect(TIMESTAMP, "^\\d{4}-\\d{2}-\\d{2}")) %>%
  mutate(
    TIMESTAMP = ymd_hms(TIMESTAMP),
    mins = as.numeric(difftime(TIMESTAMP, first(TIMESTAMP), units = "mins")),
    hrs = mins / 60,
    Tmeas = as.numeric(Temperature_C)
  )

# Apply rolling mean after time processing
ExpData <- ExpData %>%
  mutate(Tmeas_roll = rollmean(Tmeas, k = 5, fill = NA, align = "right")) %>%
  filter(!is.na(Tmeas_roll))  # Remove NAs caused by rolling mean

################################################################################

#Independent model parameters
V <- 3             # Volume (mL)
Kh <- 1            # Heat transfer coefficient
t_stop <- 0.076    # Time when heating stops (hrs)

#Dependent model parameters which you can change on a per channel basis
Qin <- 250         # Flow rate (mL/hrs)
Rin <- 50          # Heat input rate (mass/hrs)
Tb <- 0            # Background temperature (°C)
T0 <- 0.03          # Plateau temperature (°C)
#T0 <- V * Rin / Qin + Tb; Other option to double check

################################################################################

#Baseline Correction
#Edit lines 48 & 49 to select the range of data you want to examine
BaselineCorrect <- as.data.frame(ExpData$hrs[1:1000])
BaselineCorrect$data <- as.numeric(ExpData$Tmeas_roll[1:1000])
colnames(BaselineCorrect)[1] <- "hrs"

Tb_BTC <- as.numeric(min(BaselineCorrect$data)) #Background temp for BaselineCorrection (°C)
BaselineCorrect$data <- BaselineCorrect$data - Tb_BTC
#Edit line 54 to select first and last data point of the range selected
BaselineCorrect2 <- BaselineCorrect[-(2:999),]
colnames(BaselineCorrect2)[1] <- "hrs"

linecof <- lm(data ~ hrs, BaselineCorrect2)
b <- linecof$coefficients[1]
m <- linecof$coefficients[2]
BaselineCorrect$line <- BaselineCorrect$hrs * m + b
BaselineCorrect$FinalData <- BaselineCorrect$data - BaselineCorrect$line
BaselineCorrect$FinalData <- pmax(BaselineCorrect$FinalData, 0)

#plot(BaselineCorrect$hrs,BaselineCorrect$FinalData, main = "Smoothed and Corrected Data", xlab = "Time (hrs)", ylab = "Temperature (°C)")

################################################################################

# Precompute model constants
TermOne <- (V * Rin) / Qin
TermTwo <- Qin / ((1 + Kh) * V)
Temp_heating_end <- Tb + TermOne * (1 - exp(-TermTwo * t_stop))

# Recalculate mask vectors after filtering
t <- BaselineCorrect$hrs
tmaskHeating <- ifelse(t < t_stop, 1, 0)
tmaskPostHeating <- ifelse(t >= t_stop, 1, 0)

# Model functions
MixModelHeating <- function(BaselineCorrect, t_stop, Tb, Qin, Kh, V, Rin, tmaskHeating, Temp_heating_end) {
  t <- BaselineCorrect$hrs
  TcalcHeat <- Tb + (V * Rin) / Qin * (1 - exp((-1) * Qin / ((1 + Kh) * V) * t))
  TcalcHeat * tmaskHeating
}

MixModelPostHeating <- function(BaselineCorrect, t_stop, Tb, Qin, Kh, V, Rin, tmaskPostHeating, Temp_heating_end) {
  t <- BaselineCorrect$hrs
  TermOne <- (V * Rin) / Qin
  TermTwo <- Qin / ((1 + Kh) * V)
  Temp_heating_end <- Tb + TermOne * (1 - exp(-TermTwo * t_stop))
  delta_t <- t - t_stop
  TcalcPostHeat <- Tb + (Temp_heating_end - Tb) * exp((-1) * Qin / ((1 + Kh) * V) * delta_t)
  TcalcPostHeat * tmaskPostHeating
}

MixModelTotal <- function(BaselineCorrect, t_stop, Tb, Qin, Kh, V, Rin, tmaskHeating, tmaskPostHeating, Temp_heating_end) {
  Tcalc <- MixModelHeating(BaselineCorrect, t_stop, Tb, Qin, Kh, V, Rin, tmaskHeating, Temp_heating_end) +
    MixModelPostHeating(BaselineCorrect, t_stop, Tb, Qin, Kh, V, Rin, tmaskPostHeating, Temp_heating_end)
  return(Tcalc)
}

################################################################################

# Run the initial model
BaselineCorrect$Tcalc <- MixModelTotal(BaselineCorrect, t_stop, Tb, Qin, Kh, V, Rin, 
                                       tmaskHeating, tmaskPostHeating, Temp_heating_end)

# Plot initial model
plot(BaselineCorrect$hrs, BaselineCorrect$Tcalc, type = "l", col = "blue",
     main = "Breakthrough Curve - Initial Model",
     xlab = "Time (hrs)", ylab = "Temperature (°C)",
     ylim = c(min(BaselineCorrect$FinalData), max(BaselineCorrect$FinalData)))
lines(BaselineCorrect$hrs, BaselineCorrect$FinalData, col = "grey", lty = 1, lwd = 1)
abline(v = t_stop, col = "red", lty = 2)
text(0.5, 0.673, paste0('Qin = ', round(Qin, 2), '\nRin = ', round(Rin, 5)))

################################################################################

#Weights to force line on the falling side
BaselineCorrect$Weight <- 1
dfrows <- nrow(BaselineCorrect)
i <- 1 
while(BaselineCorrect$hrs[i] < t_stop){
  BaselineCorrect$Weight[i] <- 0
  i = i + 1
}

# Nonlinear fitting using smoothed data
Weight <- BaselineCorrect$Weight
nlmod <- nls(BaselineCorrect$FinalData ~ MixModelTotal(BaselineCorrect, t_stop, Tb, Qin, Kh, V, Rin, 
                                                       tmaskHeating, tmaskPostHeating, Temp_heating_end),
             data = BaselineCorrect,
             start = list(Qin = Qin, Rin = Rin),
             control = nls.control(maxiter = 100),
             weights = Weight)

# Update parameters with fitted values
Rin <- coef(nlmod)["Rin"]
Qin <- coef(nlmod)["Qin"]

# Recompute model terms
TermOne <- (V * Rin) / Qin
TermTwo <- Qin / ((1 + Kh) * V)
Temp_heating_end <- Tb + TermOne * (1 - exp(-TermTwo * t_stop))

# Recalculate model
BaselineCorrect$Tcalc <- MixModelTotal(BaselineCorrect, t_stop, Tb, Qin, Kh, V, Rin, 
                                       tmaskHeating, tmaskPostHeating, Temp_heating_end)

# Plot fitted model
plot(BaselineCorrect$hrs, BaselineCorrect$Tcalc, type = "l", col = "blue", lwd = 2,
     main = "Breakthrough Curve - NLS",
     xlab = "Time (hrs)", ylab = "Temperature (°C)",
     ylim = c(min(BaselineCorrect$FinalData), max(BaselineCorrect$FinalData)))
lines(BaselineCorrect$hrs, BaselineCorrect$FinalData, col = "grey", lty = 1, lwd = 1)
abline(v = t_stop, col = "red", lty = 2)
text(0.5, 0.673, paste0('Qin = ', round(Qin, 2), '\nRin = ', round(Rin, 5)))

# Print coefficients
coef(nlmod)