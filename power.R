# Power ########################################################################
# Packages
library(tidyverse)
library(gridExtra)

# Data
load(file = "data/delta_25.RData")

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
  x_label,          # Name of x-axis
  power = 0.8       # One specific power which should stand out
) {
  
  # Create data frame, which can be grouped by method or group
  power_plot <- stack(data) %>% 
    rename(Power = values, Method = ind) %>% 
    mutate(
      Group = rep("Duration", times = nrow(data)*ncol(data)), 
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
      x = rep(x, times = ncol(data))
    )
  
  # Plot
  if (smooth == TRUE) {
    plot <- ggplot(data = power_plot) + 
      stat_smooth(
        mapping = aes_string(x = "x", y = "Power", group = "Method", 
                             color = group),
        method = "loess", formula = y ~ x, level = 0) + 
      geom_line(mapping = aes(x = x, y = power)) + 
      scale_color_viridis_d() + 
      xlab(label = x_label) + 
      ggtitle(label = title)
  } else {
    plot <- ggplot(data = power_plot) + 
      geom_line(mapping = aes_string(x = "x", y = "Power", group = "Method", 
                                     color = group)) + 
      geom_line(mapping = aes(x = x, y = power)) + 
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
power_delta_30 <- calculate_power(data = NULL)


# Plot
n <- c(seq(from = 50, to = 600, by = 25), 700, 800, 900, 1000)
delta <- seq(from = 0, to = 0.3, by = 0.05)
plot_n_200 <- plot_power(data = power_n_200, x = delta, group = "Method", 
           title = "n = 200", x_label = expression(delta), 
           smooth = FALSE, power = 0.8)
# expression(delta ~ "= 0.30") If plots with fixed delta for title

# X-values for specific power
calculate_x_values(power = 0.8, x = x, data = power_delta_30, 
                              smooth = FALSE)
x_delta <- bind_cols(x_delta_10, x_delta_15, x_delta_20, x_delta_25, 
                     x_delta_30) %>% 
  rename(delta.10 = x_values, delta.15 = x_values1, delta.20 = x_values2, 
         delta.25 = x_values3, delta_30 = x_values4) %>% 
  select(-Method1, -Method2, -Method3, -Method4)
bind_cols(x_n, x_delta) %>% select(-Method1)

# Power = 0.8
x_n_plot <- x_n %>% 
  t %>% 
  as_tibble(.name_repair = ~ x_n$Method) %>% 
  slice(2:ncol(x_n)) %>% 
  stack %>% 
  as_tibble %>%
  rename(delta = values, Method = ind) %>%
  mutate(n = rep(c(200, 340, 750), times = nrow(x_n)), 
         delta = as.numeric(delta))
x_delta_plot <- x_delta %>% 
  t %>% 
  as_tibble(.name_repair = ~ x_delta$Method) %>% 
  slice(2:ncol(x_delta)) %>% 
  stack %>% 
  as_tibble %>%
  rename(n = values, Method = ind) %>%
  mutate(delta = rep(c(0.1, 0.15, 0.2, 0.25, 0.3), times = nrow(x_n)), 
         n = as.numeric(n))
power_80_n <- ggplot(x_n_plot) + 
  geom_line(mapping = aes(x = n, y = delta, group = Method, color = Method)) + 
  scale_color_viridis_d() +
  ylab(label = expression(delta)) + 
  ggtitle(label = "Power = 0.8")

# Saving
n <- grid.arrange(plot_n_200, plot_n_340, plot_n_750)
plot_delta <- grid.arrange(plot_delta_00, plot_delta_05, plot_delta_10, 
                           plot_delta_15, plot_delta_20, plot_delta_25, 
                           plot_delta_30)
ggsave(filename = "n.svg", path = "plots", plot = n)
save(list = "n", file = "data/n.RData")
