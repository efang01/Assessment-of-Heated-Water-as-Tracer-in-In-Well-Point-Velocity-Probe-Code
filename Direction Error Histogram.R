################################################################################
#
#   Histogram of Directional Error (apparent angle - expected 45 deg)
#   4-Channel IWPVP
#
################################################################################

library(ggplot2)

# Load data
df <- read.csv("4ChannelData2Rotate.csv")

# Signed deviation from the expected 45 deg angle
df$signed <- df$FlowRotate - 45

# Histogram
ggplot(df, aes(x = signed)) +
  # +/-15 deg accuracy band (Osorno et al., 2018)
  annotate("rect", xmin = -15, xmax = 15, ymin = 0, ymax = Inf,
           fill = "darkgrey", alpha = 0.50) +
  geom_histogram(breaks = seq(-40, 15, by = 5),
                 fill = "#4C72B0", colour = "black", alpha = 0.85) +
  # Expected angle line at 0 deviation
  geom_vline(xintercept = 0, color = "red",
             linetype = "dashed", linewidth = 0.8) +
  scale_x_continuous(breaks = seq(-40, 15, by = 5)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1)) +
  labs(
    x = "Deviation from Expected Angle",
    y = "Number of Tests",
    title = "Distribution of Directional Error in Multi-Channel IWPVP"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))