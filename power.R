# Power ########################################################################
# Packages
library(tidyverse)
library(gridExtra)

# Data
load("Data/delta_000.RData")
delta_000 <- result_df
load("Data/delta_010.RData")
delta_010 <- result_df
load("../Stat Praktikum/delta_013.RData")
delta_013 <- result_df
load("../Stat Praktikum/delta_015.RData")
delta_015 <- result_df
load("../Stat Praktikum/delta_020.RData")
delta_020 <- result_df
load("../Stat Praktikum/delta_025.RData")
delta_025 <- result_df
load("../Stat Praktikum/delta_100.RData")
delta_100 <- result_df
rm(result_df)
load("Data/N_100.RData")
load("Data/N_200.RData")
load("Data/N_400.RData")
load("Data/N_600.RData")
load("Data/delta_01.RData")
# rm(list = ls(pattern = "n200"))

# Function for power calculation ###############################################
calculate_power <- function(data = Null) {
  
  # A list with all Data Frames of the global environment
  global <- ls(envir = .GlobalEnv)
  data_list <- lapply(
    X = global[sapply(
      X = global,
      FUN = function(object) {is.data.frame(get(object))})],
    FUN = get
  )
  
  # Only treatment p values
  p_treatment_list <- lapply(
    X = data_list, 
    FUN = select, 
    neg_bin_p_value_treatment, 
    # neg_bin_log_p_value_treatment, 
    logrank_p_value, 
    cox_p_value_treatment, 
    # cox_p_value_log_treatment, 
    # cox3_p_value_treatment, 
    logit_p_value_treatment, 
    # logit_p_value_log_treatment, 
    chi_square_p_value
  )
  
  # Data frame with power.
  # Each column one Dataset.
  # Each row one test
  as_tibble(t(bind_cols(lapply(X = p_treatment_list, FUN = function(data) {
    apply(X = data, MARGIN = 2, FUN = function(column) {
      nrow(filter(.data = data, column <= 0.05)) / nrow(data)
    })
  }))), .name_repair = ~ c(
    "Negative.Binomial", 
    # "Negative.Binomial.Log", 
    "Log.Rank", 
    "Cox", 
    # "Cox.Log", 
    # "Cox.only.treatment", 
    "Logit",
    # "Logit.Log",
    "Chi.Square"
  ))
}

# Function to plot power curves ################################################
plot_power <- function(
  data = power,     # Name of data frame with power values
  x = x,            # Scale of x axis
  group = "Method", # "Method": Each line gets individual color
  # "Group": Colors are grouped in duration, count and binary
  smooth = FALSE,   # TRUE the lines are smoothed 
  title,            # Title of plot
  x_label           # Name of x-axis
) {
  
  # Create data frame, which can be grouped by method or group
  power_plot <- stack(power) %>% 
    rename(Power = values, Method = ind) %>% 
    mutate(
      Group = rep("Duration", times = length(Power)), 
      Group = if_else(
        condition = str_detect(string = Method, pattern = 'Logit|Chi'), 
        true = "Binary", 
        false = Group
      ), 
      Group = if_else(
        condition = str_detect(string = Method, pattern = 'Negative.Binomial'), 
        true = "Count", 
        false = Group
      ), 
      x = rep(x, times = ncol(power))
    )
  
  # Plot
  if (smooth == TRUE) {
    plot <- ggplot(data = power_plot) + 
      stat_smooth(
        mapping = aes_string(x = "x", y = "Power", group = "Method", 
                             color = group),
        method = "loess", formula = y ~ x, level = 0) + 
      geom_line(mapping = aes(x = x, y = 0.8)) + 
      scale_color_viridis_d() + 
      xlab(label = x_label) + 
      ggtitle(label = title)
  } else {
    plot <- ggplot(data = power_plot) + 
      geom_line(mapping = aes_string(x = "x", y = "Power", group = "Method", 
                                     color = group)) + 
      geom_line(mapping = aes(x = x, y = 0.8)) + 
      scale_color_viridis_d() + 
      xlab(label = x_label) +
      ggtitle(label = title)
  }
  x_plot <<- ggplot_build(plot = plot)$data[[1]]$x
  y_plot <<- ggplot_build(plot = plot)$data[[1]]$y
  plot
}

# Function to calculate values of x-axis for a specific power ##################
calculate_x_values <- function(power = 0.8, x = x, data = power, 
                               smooth = FALSE) {
  
  if (smooth == TRUE) {
    
    intersect <- which(diff(y_plot > 0.8) > 0)
    
    y1 <- y_plot[intersect]
    y2 <- y_plot[intersect + 1]
    
    slope <- (y2 - y1) / (x_plot[intersect + 1] - x_plot[intersect])
    
    as_tibble(
      data.frame(x_values = x_plot[intersect] + ((0.8 - y1) / (slope - 0)))
    )
    
  } else {
    index <- 1:ncol(data)
    
    # Vectors are split into values above and below a power
    intersect <- unlist(apply(X = diff(data > power) != 0, 
                              MARGIN = 2, FUN = which))
    
    # Point below the split.
    y1 <- unlist(sapply(X = index, FUN = function(i) {data[intersect[i], i]}))
    y1 <- y1[!is.na(y1)]
    # Point above the split.
    y2 <- unlist(sapply(X = index, FUN = function(i) {data[intersect[i] + 1, i]}))
    y2 <- y2[!is.na(y2)]
    
    # Slope of intersection
    slope <- (y2 - y1) / (x[intersect + 1] - x[intersect])
    
    # Intersection
    as_tibble(rownames_to_column(
      data.frame(x_values = x[intersect] + ((0.8 - y1) / (slope - 0)))
    )) %>% rename(Method = rowname)
  }
}

# Execution ####################################################################
# Calculate power
remove(epilepsy)
power <- calculate_power(data = NULL)
# save(list = "power_100", file = "Data/power_100.RData", envir = .GlobalEnv)
power

x <- seq(from = 100, to = 1000, by = 50)

# Plot
plot_power(data = power, x = x, group = "Method", title = "Delta = 0.1", 
           x_label = "n", smooth = TRUE)
ggsave(filename = "delta_01_smooth.pdf", path = "plots")

# X-values for specific power
calculate_x_values(power = 0.8, x = x, data = power, smooth = TRUE)

# Other Plots ##################################################################
epi_seizures <- ggplot(data = epilepsy) + 
  geom_density(mapping = aes(x = seizures_treatment, fill = 1)) + 
  guides(fill = FALSE) + 
  theme_classic()
dat_seizures <- ggplot(data = dataset[[11]]) + 
  geom_density(mapping = aes(x = seizures_treatment, fill = 1)) + 
  guides(fill = FALSE) +
  theme_classic()
epi_time <- ggplot(data = epilepsy) + 
  geom_density(mapping = aes(x = time_baseline, fill = 1)) + 
  guides(fill = FALSE) + 
  theme_classic()
dat_time <- ggplot(data = dataset[[11]]) + 
  geom_density(mapping = aes(x = time_baseline, fill = 1)) + 
  guides(fill = FALSE) + 
  theme_classic()
epi_point <- ggplot(data = epilepsy) + 
  geom_point(mapping = aes(x = time_study, y = seizures_treatment, color = 1)) + 
  guides(color = FALSE) + 
  theme_classic()
dat_point <- ggplot(data = dataset[[11]]) + 
  geom_point(mapping = aes(x = time_study, y = seizures_treatment, color = 1)) + 
  guides(color = FALSE) + 
  theme_classic()

# jpeg(filename = "plots/plot.jpg")
grid.arrange(epi_seizures, dat_seizures, epi_time, dat_time, epi_point, 
             dat_point)
# dev.off()

ggplot(data = power) + 
  geom_line(mapping = aes(x = x, y = neg_bin_p_value_treatment, 
                          color = "Negative Binomial")) + 
  geom_line(mapping = aes(x = x, y = neg_bin_log_p_value_treatment, 
                          color = "Negative Binomial Log")) +
  geom_line(mapping = aes(x = x, y = logrank_p_value, color = "Log Rank")) +
  geom_line(mapping = aes(x = x, y = cox_p_value_treatment, 
                          color = "Cox (full)")) +
  geom_line(mapping = aes(x = x, y = cox_p_value_log_treatment, 
                          color = "Cox Log (full)")) +
  geom_line(mapping = aes(x = x, y = cox3_p_value_treatment, 
                          color = "Cox (only treatment)")) +
  geom_line(mapping = aes(x = x, y = logit_p_value_treatment, 
                          color = "Logit")) +
  geom_line(mapping = aes(x = x, y = logit_p_value_log_treatment, 
                          color = "Logit Log")) +
  geom_line(mapping = aes(x = x, y = chi_square_p_value, 
                          color = "Chi Square")) +
  geom_line(mapping = aes(x = x, y = 0.05, color = "0.05")) +
  geom_line(mapping = aes(x = x, y = 0.8, color = "0.8")) + 
  scale_color_viridis_d(name = "Method") +
  xlab(label = "delta") +
  ylab(label = "Power") +
  ggtitle(label = "N = 400")