# Set seed for reproducibility
set.seed(1)

require(ggplot2)
require(hexSticker)

# Generate sample data for a smooth functional curve with credible intervals
q <- 100
time <- seq(1, q, length.out = q)
mean_curve <- 1+sin(time / 17)+1.2*cos(time/12.5)
lower_band <- mean_curve - 0.9
upper_band <- mean_curve + 0.8

# Generate fewer points for added visual interest
nb_points <- 6
points_indices <- sample(1:q, nb_points)
points_data <- data.frame(
  time = time[points_indices],
  value = mean_curve[points_indices] + rnorm(nb_points, mean = 0, sd = 0.5)
)

# Ensure most points fall within the credible band
points_data$value[points_data$value < lower_band[points_indices]] <- lower_band[points_indices]
points_data$value[points_data$value > upper_band[points_indices]] <- upper_band[points_indices]


curve_data <- data.frame(
  time = time,
  mean = mean_curve,
  lower = lower_band,
  upper = upper_band
)

traj_data <- data.frame(
  time = time,
  mean = 1+sin(time / 17)+1.6*cos(time/12)
)

traj_data_2 <- data.frame(
  time = time,
  mean = 0.3+sin(time / 20)+1.2*cos(time/14)
)

traj_data_3 <- data.frame(
  time = time,
  mean = 0.9+sin(time / 19)+1.6*cos(time/13)
)

traj_data_4 <- data.frame(
  time = time,
  mean = 1.8+sin(time / 18.5)+1.3*cos(time/13.5)
)

traj_data_5 <- data.frame(
  time = time,
  mean = 0.7+sin(time / 18)+1.2*cos(time/12.5)
)

traj_data_6 <- data.frame(
  time = time,
  mean = 0.5+sin(time / 19)+2*cos(time/13)
)

traj_data_7 <- data.frame(
  time = time,
  mean = 0.2+sin(time / 13)+1.5*cos(time/12)
)

traj_data_8 <- data.frame(
  time = time,
  mean = 0.6+sin(time / 15)+1.6*cos(time/11)
)

traj_data_9 <- data.frame(
  time = time,
  mean = 1.25+sin(time / 17)+2.2*cos(time/14)
)

col_traj <- "grey15"
col_factor <- "#88BDE6"

# Plot the functional curves with credible intervals and points
plot_curves <- ggplot() +
  geom_line(data = traj_data, aes(x = time, y = mean), color = col_traj, size = 0.4, alpha = 1) + # darkblue
  geom_line(data = traj_data_2, aes(x = time, y = mean), color = col_traj, size = 0.4, alpha = 1) + # darkblue
  geom_line(data = traj_data_3, aes(x = time, y = mean), color = col_traj, size = 0.4, alpha = 1) + # darkblue
  geom_line(data = traj_data_4, aes(x = time, y = mean), color = col_traj, size = 0.4, alpha = 1) + # darkblue
  geom_line(data = traj_data_5, aes(x = time, y = mean), color = col_traj, size = 0.4, alpha = 1) + # darkblue
  geom_line(data = traj_data_6, aes(x = time, y = mean), color = col_traj, size = 0.4, alpha = 1) + # darkblue
  geom_line(data = traj_data_7, aes(x = time, y = mean), color = col_traj, size = 0.4, alpha = 1) + # darkblue
  geom_line(data = traj_data_8, aes(x = time, y = mean), color = col_traj, size = 0.4, alpha = 1) + # darkblue
  geom_line(data = traj_data_9, aes(x = time, y = mean), color = col_traj, size = 0.4, alpha = 1) + # darkblue
  geom_ribbon(data = curve_data, aes(x = time, ymin = lower, ymax = upper), fill = "gray30", alpha = 0.4) +
  geom_line(data = curve_data, aes(x = time, y = mean), color = col_factor, size = 1.5) + # darkblue
  geom_point(data = points_data, aes(x = time, y = value), color = "white", size = 1.5) +
  theme_void() +
  theme(legend.position = "none")

# Create the hex sticker
sticker(
  subplot = plot_curves,
  package = "bayesSYNC",
  p_size = 64,  # Larger text size
  s_x = 1,
  s_y = 0.89,
  s_width = 1.6,
  s_height = 1.5,
  p_x = 0.85,
  p_y = 1.44,
  u_color = "white",
  u_size = 1,
  h_fill = "gray70",  # lighter grey background
  h_color = "gray70",
  filename = "man/figures/bayesSYNC_logo.png",
  dpi = 1200
)
