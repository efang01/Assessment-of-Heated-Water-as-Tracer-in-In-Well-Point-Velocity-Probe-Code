# ============================================================
# Plot: Probe Velocity vs Sand Tank Velocity
# Colored by probe rotation angle
# Data sets: 0°, 45°, 90°
# ============================================================

library(ggplot2)

# --- Load data -----------------------------------------------
ds1 <- read.csv("4ChannelData1ExpDAY027.csv")
ds2 <- read.csv("4ChannelData2ExpDAY027.csv")
ds3 <- read.csv("4ChannelData3ExpDAY027.csv")

# Keep only the two relevant columns (ds3 has trailing empty cols)
ds1 <- ds1[, c("ProbeFlux", "SandFlux")]
ds2 <- ds2[, c("ProbeFlux", "SandFlux")]
ds3 <- ds3[, c("ProbeFlux", "SandFlux")]

# --- Tag each dataset with its rotation label ----------------
ds1$Rotation <- "0°"
ds2$Rotation <- "45°"
ds3$Rotation <- "90°"

# --- Combine into one data frame -----------------------------
all_data <- rbind(ds1, ds2, ds3)

# Set factor order so legend appears 0° → 45° → 90°
all_data$Rotation <- factor(all_data$Rotation, levels = c("0°", "45°", "90°"))

# --- Plot ----------------------------------------------------
ggplot(all_data, aes(x = SandFlux, y = ProbeFlux, color = Rotation, shape = Rotation)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(
    name  = "Probe Rotation",
    values = c("0°"  = "#1b7837",   # green
               "45°" = "#762a83",   # purple
               "90°" = "#d6604d")   # red-orange
  ) +
  scale_shape_manual(
    name  = "Probe Rotation",
    values = c("0°" = 16, "45°" = 17, "90°" = 15)  # circle, triangle, square
  ) +
  labs(
    title = "Probe Velocity vs. Sand Tank Velocity by Rotation Angle",
    x     = "Sand Tank Velocity (cm/day)",
    y     = "Probe Flux (cm³/day)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position   = c(0.15, 0.80),
    legend.background = element_rect(fill = "white", color = "grey70"),
    legend.key        = element_rect(fill = "white"),
    plot.title        = element_text(hjust = 0.5, face = "bold")
  )