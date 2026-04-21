################################################################################
# 
#
#       Determining Flow Direction and Velocity in Heat IWPVP
#       Written by Ellis Fangmann
#       Feb 12, 2026
#
#
################################################################################
# DONT FORGET TO SELECT CHANNELS!

# Load and prepare data set 1
ExpData1 <- read.table("CR1000_TempSensApril9Test1.dat", header = TRUE, sep = ",", skip = 1)

ExpData1 <- ExpData1 %>%
  filter(str_detect(TIMESTAMP, "^\\d{4}-\\d{2}-\\d{2}")) %>%
  mutate(
    TIMESTAMP = ymd_hms(TIMESTAMP),
    mins = as.numeric(difftime(TIMESTAMP, first(TIMESTAMP), units = "mins")),
    hrs = mins / 60,
    Tmeas = as.numeric(Temperature_C.2.) #Change the # in Temperature_C.#. to pick which channel you want analyzed
  )

# Load and prepare data set 2
ExpData2 <- read.table("CR1000_TempSensApril9Test1.dat", header = TRUE, sep = ",", skip = 1)

ExpData2 <- ExpData2 %>%
  filter(str_detect(TIMESTAMP, "^\\d{4}-\\d{2}-\\d{2}")) %>%
  mutate(
    TIMESTAMP = ymd_hms(TIMESTAMP),
    mins = as.numeric(difftime(TIMESTAMP, first(TIMESTAMP), units = "mins")),
    hrs = mins / 60,
    Tmeas = as.numeric(Temperature_C.3.) #Change the # in Temperature_C.#. to pick which channel you want analyzed
)

# Apply rolling mean after time processing
ExpData1 <- ExpData1 %>%
  mutate(Tmeas_roll = rollmean(Tmeas, k = 5, fill = NA, align = "right")) %>%
  filter(!is.na(Tmeas_roll))  # Remove NAs caused by rolling mean

# Apply rolling mean after time processing
ExpData2 <- ExpData2 %>%
  mutate(Tmeas_roll = rollmean(Tmeas, k = 5, fill = NA, align = "right")) %>%
  filter(!is.na(Tmeas_roll))  # Remove NAs caused by rolling mean

################################################################################

# Calculate direction of water and velocity 

# Calculate area under curve for each channel
#Background temp for BTC (°C), don't forget to match the channels selected above
Tb_BTC1 <- as.numeric(min(ExpData1$Temperature_C.2.)) 
Tb_BTC2 <- as.numeric(min(ExpData2$Temperature_C.3.))

#Temps for BTC (°C), don't forget to match the channels selected above
T_BTC1 <- as.numeric(ExpData1$Temperature_C.2.) 
T_BTC2 <- as.numeric(ExpData2$Temperature_C.3.) 

################################################################################

#Correct the baseline for W1 + Qin1
#Edit lines 63 & 64 to select the range of data you want to examine
BaselineCorrectW1 <- as.data.frame(ExpData1$hrs[1:1000])
BaselineCorrectW1$data <- as.numeric(ExpData1$Tmeas_roll[1:1000])
colnames(BaselineCorrectW1)[1] <- "hrs"

BaselineCorrectW1$data <- BaselineCorrectW1$data - Tb_BTC1
#Edit line 69 to select first and last data point of the range selected
BaselineCorrectW12 <- BaselineCorrectW1[-(2:999),]

colnames(BaselineCorrectW12)[1] <- "hrs"

linecof1 <- lm(data ~ hrs, BaselineCorrectW12)
b1 <- linecof1$coefficients[1]
m1 <- linecof1$coefficients[2]
BaselineCorrectW1$line <- BaselineCorrectW1$hrs * m1 + b1
BaselineCorrectW1$FinalData <- BaselineCorrectW1$data - BaselineCorrectW1$line
BaselineCorrectW1$FinalData <- pmax(BaselineCorrectW1$FinalData, 0)

#plot(BaselineCorrectW1$hrs,BaselineCorrectW1$FinalData)
#plot(BaselineCorrectW1$hrs,BaselineCorrectW1$data)

#Correct the baseline for W2 + Qin2
#Edit lines 85 & 86 to select the range of data you want to examine
BaselineCorrectW2 <- as.data.frame(ExpData2$hrs[1:1000])
BaselineCorrectW2$data <- as.numeric(ExpData2$Tmeas_roll[1:1000])
colnames(BaselineCorrectW2)[1] <- "hrs"

BaselineCorrectW2$data <- BaselineCorrectW2$data - Tb_BTC2
#Edit line 91 to select first and last data point of the range selected
BaselineCorrectW22 <- BaselineCorrectW2[-(2:999),]

colnames(BaselineCorrectW22)[1] <- "hrs"

linecof2 <- lm(data ~ hrs, BaselineCorrectW22)
b2 <- linecof2$coefficients[1]
m2 <- linecof2$coefficients[2]
BaselineCorrectW2$line <- BaselineCorrectW2$hrs * m2 + b2
BaselineCorrectW2$FinalData <- BaselineCorrectW2$data - BaselineCorrectW2$line
BaselineCorrectW2$FinalData <- pmax(BaselineCorrectW2$FinalData, 0)

#plot(BaselineCorrectW2$hrs,BaselineCorrectW2$FinalData)
#plot(BaselineCorrectW2$hrs,BaselineCorrectW2$data)

################################################################################

Δt1 <- diff(BaselineCorrectW1$hrs) #Time step (hrs)
Δt2 <- diff(BaselineCorrectW2$hrs)

#Area of total curve
Atot1 <- sum((BaselineCorrectW1$FinalData[-1]) * Δt1)
Atot2 <- sum((BaselineCorrectW2$FinalData[-1]) * Δt2)

#Calculate weights for each channel
W1 <- Atot1/(Atot1+Atot2)
W2 <- Atot2/(Atot1+Atot2) 

#Calculate apparent theta angle
#Edit lines 120 & 121 to input the calculated Qins for BOTH channels
Qin1 <-  1355.95 #Velocity measured (mL/hr)
Qin2 <-  1229.28 #Velocity measured (mL/hr)

#The angle is from whichever weight and velocity are on the num/Qin1/W1 of equation!!!
num <- Qin1 * W1
den <- Qin2 * W2
Θapp_rad <- atan((num)/(den))
Θapp_deg <- Θapp_rad * 180 / pi

#Calculate velocity of IWPVP
Viwpvp <- sqrt((Qin1 * sin(Θapp_rad))^2 + (Qin2 * cos(Θapp_rad))^2)

# Print values
print(Viwpvp)
print(Θapp_deg)