################################################################################
# 
#
#       Rose Diagram Generator 
#       Written by Ellis Fangmann
#       April 20, 2026
#
#
################################################################################
library(ggplot2)
library(dplyr)
library(patchwork)

df <- read.csv("4ChannelData2Rotate.csv")

make_polar <- function(flow_rot, test_id, pump_rate, wedge_width = 5) {
  angle_offset <- 45
  rose_xy <- function(angle_deg, r = 1) {
    data.frame(
      x = r * sin((angle_deg + angle_offset) * pi / 180),
      y = r * cos((angle_deg + angle_offset) * pi / 180)
    )
  }
  
  rings <- lapply(seq(0.2, 1, by = 0.1), function(r) {
    a <- seq(0, 360, length.out = 300)
    data.frame(
      x = r * sin(a * pi / 180),
      y = r * cos(a * pi / 180),
      r = r
    )
  }) |> bind_rows()
  
  spokes <- data.frame(angle = seq(0, 350, by = 10)) |>
    rowwise() |>
    do(data.frame(
      x = c(0, sin(.$angle * pi / 180)),
      y = c(0, cos(.$angle * pi / 180)),
      g = .$angle
    ))
  
  wedge_angles <- seq(
    flow_rot - wedge_width / 2,
    flow_rot + wedge_width / 2,
    length.out = 40
  )
  
  wedge <- rose_xy(wedge_angles, r = 1)
  wedge <- rbind(data.frame(x = 0, y = 0), wedge)
  
  ggplot() +
    geom_path(data = rings, aes(x, y, group = r), color = "black", linewidth = 0.25) +
    geom_path(data = spokes, aes(x, y, group = g), color = "black", linewidth = 0.25) +
    geom_segment(aes(x = -1.05, xend = 1.05, y = 0, yend = 0), linewidth = 0.75) +
    geom_segment(aes(x = 0, xend = 0, y = -1.05, yend = 1.05), linewidth = 0.75) +
    geom_polygon(data = wedge, aes(x, y), fill = "red", alpha = 0.75) +
    coord_fixed(xlim = c(-1.08, 1.08), ylim = c(-1.08, 1.08)) +
    labs(
      title = paste("Test", test_id),
      caption = paste0(pump_rate, " mL/min")   #pump rate added here
    ) +
    theme_void() +
    theme(
      plot.title = element_text(
        size = 9,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 2)
      ),
      plot.caption = element_text(    #style the caption
        size = 9,
        face = "bold",
        hjust = 0.5,
        margin = margin(t = 3)
      )
    )
}

df_plot <- df[1:36, ]

all_tests <- lapply(1:nrow(df_plot), function(i) {
  make_polar(
    df_plot$FlowRotate[i],
    i,
    df_plot$PumpRate[i]   # <-- pass pump rate
  )
})

final_plot <- wrap_plots(all_tests, ncol = 9)
final_plot