# Test #########################################################################
n <- 200
delta <- 0.13

# seizures_baseline with individual lambda and more than three seizures. #######
output <- tibble(
  lambda_baseline = rgamma(n = n * 1.5, shape = 5, scale = 2),
  seizures_baseline = sapply(X = lambda_baseline, FUN = rpois, n = 1)
) %>% 
  filter(seizures_baseline > 3) %>% 
  slice(1:n)
lambda_baseline <- output$lambda_baseline
seizures_baseline <- output$seizures_baseline
ggplot() + 
  geom_density(mapping = aes(x = seizures_baseline, colour = 1, fill = 2)) + 
  guides(colour = FALSE, fill = FALSE)

# treatment ####################################################################
treatment = as.factor(rbinom(n = n, size = 1, prob = 0.5))
ggplot() + 
  geom_bar(mapping = aes(x = treatment, colour = 1, fill = 2)) + 
  guides(colour = FALSE, fill = FALSE)

# time_study ###################################################################
time_study = round(pmin(rexp(n = length(seizures_baseline),
                             rate = -log(0.8) / 56), 56))
time_study = if_else(condition = time_study != 0, true = time_study, 
                     false = 1)
ggplot() + 
  geom_density(mapping = aes(x = time_study, colour = 1, fill = 2)) + 
  guides(colour = FALSE, fill = FALSE)

# seizures_baseline and time_baseline with individual lambda ###################
lambda_treatment <- lambda_baseline / 
  (exp((as.numeric(treatment) - 1) * delta + 0.2) * 28)
seizures_treatment <- c()
time_baseline <- time_study
for (i in 1:n) {
  count <- 0
  time <- 0
  while (time < time_study[i]) {
    if (count == seizures_baseline[i]) {
      time_baseline[i] <- time
    }
    time <- time + rexp(n = 1, rate = lambda_treatment[i])
    count <- count + 1
  }
  seizures_treatment[i] <- count - 1
}

duration_times <- lapply(
  X = lapply(X = lambda_treatment, FUN = rexp, n = 60), 
  FUN = cumsum
)
seizures_treatment <- mapply(
  FUN = function(x, y) {length(x[x <= y])}, 
  x = duration_times, 
  y = time_study
) 
time_baseline <- round(unlist(lapply(
  X = mapply(FUN = function(x, y) {x[1:y]}, x = duration_times, 
             y = seizures_baseline), 
  FUN = function(x) {x[length(x)]}))
)
time_baseline <- if_else(condition = time_baseline <= time_study, 
                         true = time_baseline, false = time_study)

summary(object = seizures_treatment)
summary(object = epilepsy$seizures_treatment)
summary(time_baseline)
summary(object = epilepsy$time_baseline)
ggplot() + 
  geom_density(mapping = aes(x = time_baseline, colour = 1, fill = 2)) + 
  guides(colour = FALSE, fill = FALSE)
ggplot() + 
  geom_density(mapping = aes(x = seizures_treatment, colour = 1, fill = 2)) + 
  guides(colour = FALSE, fill = FALSE)
ggplot() + 
  geom_point(mapping = aes(x = time_study, y = seizures_treatment, 
                           colour = 1)) + 
  guides(colour = FALSE)
seizures_treatment[26]
time_baseline[26]
time_baseline[[400]][length(time_baseline[[400]])]
seizures_treatment <- mapply(
  FUN = function(x, y) {length(x[x < y])}, 
  x = time, 
  y = time_study
)
time[[400]][1:seizures_baseline[400]]
time_baseline <- round(unlist(mapply(
  FUN = function(x, y) {length(x[x == y])},
  x = time, 
  y = seizures_baseline
)))
time_baseline <- if_else(condition = time_baseline <= time_study, 
                         true = time_baseline, false = time_study)

seizures_treatment <- mapply(
  FUN = function(x, y) {length(x[x < y])}, 
  x = lapply(X = mapply(FUN = rexp, n = time_study, rate = lambda_treatment), 
             FUN = cumsum), 
  y = time_study
)


time_baseline <- round(unlist(lapply(
  X = mapply(FUN = rexp, n = seizures_baseline, rate = lambda_treatment),
  FUN = sum
)))
time_baseline <- if_else(condition = time_baseline <= time_study, 
                         true = time_baseline, false = time_study)

seizures_treatment <- mapply(
  FUN = function(x, y) {length(x[x < y])}, 
  x = lapply(X = mapply(FUN = rexp, n = time_study, rate = lambda_treatment), 
             FUN = cumsum), 
  y = time_study
)

any(time_baseline >= time_study)

summary(seizures_treatment)
summary(round(time_baseline))

any(time_baseline > time_study)
summary(epilepsy$seizures_treatment)
summary(epilepsy$time_baseline)
cor(seizures_treatment, time_study)

if (is.na(time_baseline[i])) {
  time_baseline[i] <- time_study[i]
}
if (time_baseline[i] > time_study[i]) {
  time_baseline[i] <- time_study[i]
}
time_baseline <- round(if_else(condition = time_baseline < time_study,
                               true = time_baseline, false = time_study))

ggplot() + 
  geom_density(mapping = aes(x = seizures_treatment, colour = 1, fill = 2)) + 
  guides(colour = FALSE, fill = FALSE)

# seizures_treatment
lambda_treatment = lambda_baseline / 
  (exp((as.numeric(treatment) - 1) * delta + 0.2) * 28)
seizures_treatment = sapply(
  X = lambda_treatment, 
  FUN = rpois, 
  n = 1
)
ggplot() + 
  geom_density(mapping = aes(x = seizures_treatment, colour = 1, fill = 2)) + 
  guides(colour = FALSE, fill = FALSE)

# time_baseline with individual lambda.
lambda_time = if_else(condition = seizures_treatment != 0, 
                      true = seizures_treatment / time_study, 
                      false = 1e-10)
time_baseline = round(unlist(lapply(X = mapply(FUN = rexp, 
                                               n = seizures_baseline, 
                                               rate = lambda_time), 
                                    FUN = sum)))
time_baseline = if_else(condition = time_baseline <= time_study, 
                        true = time_baseline, false = time_study)
time_baseline = if_else(condition = time_baseline > 0, 
                        true = time_baseline, false = 1)
ggplot() + 
  geom_density(mapping = aes(x = time_baseline, colour = 1, fill = 2)) + 
  guides(colour = FALSE, fill = FALSE)

# Power = 0.8 ##################################################################
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
calculate_x_values(power = 0.8, x = x, data = power, smooth = TRUE)
index <- 1:ncol(power)
intersect <- unlist(apply(X = diff(power > 0.8) != 0, MARGIN = 2, FUN = which))
y1 <- unlist(sapply(X = index, FUN = function(i) {power[intersect[i], i]}))
y1 <- y1[!is.na(y1)]
y2 <- unlist(sapply(X = index, FUN = function(i) {power[intersect[i] + 1, i]}))
y2 <- y2[!is.na(y2)]
slope <- (y2 - y1) / (x[intersect + 1] - x[intersect])
p <- plot_power(data = power, x = x, group = "Method", title = "Delta = 0.1", 
           x_label = "n", smooth = TRUE)
x_plot <- ggplot_build(plot = p)$data[[1]]$x
y_plot <- ggplot_build(plot = p)$data[[1]]$y
x_plot[intersect != max(x_plot)]
intersect <- which(diff(y_plot > 0.8) > 0)
y1 <- y_plot[intersect]
y2 <- y_plot[intersect + 1]
slope <- (y2 - y1) / (x_plot[intersect + 1] - x_plot[intersect])
as_tibble(rownames_to_column(
  data.frame(x_values = x_plot[intersect] + ((0.8 - y1) / (slope - 0)))
)) %>% rename(Method = rowname)
# Different Deltas and Ns ######################################################
# function to simulate for different values of delta
simulate_data <- function(number_datasets = 10000, n = 200, delta = 0.13, 
                          name = "n_200_delta_") {
  if (length(n) > 1) {
    results <- lapply(X = n, FUN = function(n) {
      bind_rows(parLapply(cl = cl, X = 1:number_datasets, fun = sampling, 
                          n = n, delta = delta))
    })
    # save each data.frame with an individual name
    lapply(X = 1:length(x), FUN = function(i) {
      nam <- results[[i]]
      if (nam$n[1] < 1000 & nam$n[1] >= 100) {
        assign(paste0(name, "0", nam$n[1], sep = ""), value = nam,
               envir = globalenv())
      } else if (nam$n[1] < 100) {
        assign(paste0(name, "00", nam$n[1], sep = ""), value = nam,
               envir = globalenv())
      } else {
        assign(paste0(name, nam$n[1], sep = ""), value = nam,
               envir = globalenv())
      }
    })
  } 
  if (length(delta) > 1) {
    results <- lapply(X = delta, FUN = function(delta) {
      bind_rows(parLapply(cl = cl, X = 1:number_datasets, fun = sampling, 
                          n = n, delta = delta))
    })
    # save each data.frame with an individual name
    lapply(X = 1:length(x), FUN = function(i) {
      nam <- results[[i]]
      assign(paste0(name, nam$delta[1], sep = ""), value = nam, 
             envir = globalenv())
    })
  }
}

results <- lapply(X = x, FUN = function(n) {
  bind_rows(parLapply(cl = cl, X = 1:10, fun = sampling, 
                      n = n, delta = 0.13))
})
lapply(X = 1:length(x), FUN = function(i) {
  nam <- results[[i]]
  assign(paste0("name", nam$n[1], sep = ""), value = nam, 
         envir = globalenv())
})

# Power ########################################################################
global <- ls(envir = .GlobalEnv)
data_list <- lapply(
  X = global[sapply(
    X = global,
    FUN = function(object) {is.data.frame(get(object))})],
  FUN = get
)
p_treatment_list <- lapply(
  X = data_list, 
  FUN = select, 
  neg_bin_p_value_treatment, 
  neg_bin_log_p_value_treatment, 
  logrank_p_value, 
  cox_p_value_treatment, 
  cox_p_value_log_treatment, 
  cox3_p_value_treatment, 
  logit_p_value_treatment, 
  logit_p_value_log_treatment, 
  chi_square_p_value
)
as_tibble(t(bind_cols(lapply(X = p_treatment_list, FUN = function(data) {
  apply(X = data, MARGIN = 2, FUN = function(column) {
    nrow(filter(.data = data, column <= 0.05)) / nrow(data)
  })
}))), .name_repair = ~ c(
  "neg_bin_p_value_treatment", 
  "neg_bin_log_p_value_treatment", 
  "logrank_p_value", 
  "cox_p_value_treatment", 
  "cox_p_value_log_treatment", 
  "cox3_p_value_treatment", 
  "logit_p_value_treatment",
  "logit_p_value_log_treatment",
  "chi_square_p_value"
))

power <- bind_cols(lapply(X = p_treatment_list, FUN = function(data) {
  apply(X = data, MARGIN = 2, FUN = function(column) {
    nrow(filter(.data = data, column <= 0.05)) / nrow(data)
  })
})) %>% 
  mutate(Method = c(
    "neg_bin_p_value_treatment", 
    "neg_bin_log_p_value_treatment", 
    "logrank_p_value", 
    "cox_p_value_treatment", 
    "cox_p_value_log_treatment", 
    "cox3_p_value_treatment", 
    "logit_p_value_treatment",
    "logit_p_value_log_treatment",
    "chi_square_p_value"
  ))
power_plot <- as_tibble(stack(power)) %>% 
  rename(Power = values, Method = ind) %>% 
  mutate(
    Group = if_else(condition = Method == "neg_bin_p_value_treatment" ||
                      "neg_bin_log_p_value_treatment", 
                    true = "count", false = "duration"), 
    Group = if_else(condition = Method == "logit_p_value_treatment" || 
                      "logit_p_value_log_treatment", 
                    true = "binary", false = Group), 
    Group = if_else(condition = Method == "chi_square_p_value", 
                    true = "binary", false = Group)
  )

if_else(condition = power_plot$Method == "neg_bin_p_value_treatment" |
          "neg_bin_log_p_value_treatment", 
        true = "count", false = "duration")

for (i in 1:nrow(power_plot)) {
  if (power_plot$Method == "logit_p_value_treatment" || 
      "logit_p_value_log_treatment") {
    group[i] <- "count"
  }
}          

x <- seq(from = 0, to = 0.75, by = 0.05)
power_plot <- stack(power) %>% 
  rename(Power = values, Method = ind) %>% 
  mutate(
    Group = rep("duration", times = length(Power)), 
    Group = if_else(
      condition = str_detect(string = Method, pattern = 'Logit|Chi'), 
      true = "binary", 
      false = Group
    ), 
    Group = if_else(
      condition = str_detect(string = Method, pattern = 'Negative.Binomial'), 
      true = "count", 
      false = Group
    ), x = rep(x, times = ncol(power))
  )
power_plot <- power_plot %>%
  filter(str_detect(string = Method, pattern = 'Logit|Chi')) %>%
  mutate(group = "binary") 
power_plot <- power_plot %>%
  filter(grepl(pattern = 'Negative.Binomial', x = power_plot$Method)) %>%
  mutate(group = "count")
str_detect(string = power_plot$Method, pattern = 'Logit|Chi')
power_plot$group[grepl(pattern = "Logit", x = power_plot$Method)] <- "binary"
power_plot$group[grep(pattern = "chi", x = power_plot$Method)] <- "binary"
power_plot$group[grep(pattern = "neg_bin", x = power_plot$Method)] <- "count"
ggplot(data = power_plot) + 
  geom_line(mapping = aes(x = x, y = Power, group = Method, color = Group)) + 
  geom_line(mapping = aes(x = x, y = 0.8)) + 
  scale_color_viridis_d() + 
  xlab(label = "Delta") +
  ggtitle(label = "N = 400")
plot_power <- function(data = power, x = x, group = "Method", title, x_label) {
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
  ggplot(data = power_plot) + 
    geom_line(mapping = aes_string(x = "x", y = "Power", group = "Method", 
                            color = group)) + 
    geom_line(mapping = aes(x = x, y = 0.8)) + 
    scale_color_viridis_d() + 
    xlab(label = x_label) +
    ggtitle(label = title)
}
plot_power(data = power, x = x, group = "Method", title = "N = 400", 
           x_label = "Delta")
ggsave(filename = "n_400.svg", path = "plots")
x <- seq(from = 100, to = 1000, by = 50)
ggplot(data = power_plot) + 
  stat_smooth(mapping = aes(x = x, y = Power, group = Method, color = Method), 
              method = "loess", formula = y ~ x, level = 0, fullrange = TRUE) + 
  geom_line(mapping = aes(x = x, y = 0.8)) + 
  scale_color_viridis_d() + 
  xlab(label = "n") + 
  ggtitle(label = "Delta = 0.1")

ggplot_build(p)$data[[1]]
  approx(x = x, y = Power)
predict(object = loess(formula = Power ~ x, data = power_plot), 
        newdata = data.frame(Power = rep(0.8, times = 19)))
ncol(power_plot)
loess(formula = Power ~ x, data = power_plot)



# plot_power ###################################################################
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
    plot
  } else {
    plot <- ggplot(data = power_plot) + 
      geom_line(mapping = aes_string(x = "x", y = "Power", group = "Method", 
                                     color = group)) + 
      geom_line(mapping = aes(x = x, y = 0.8)) + 
      scale_color_viridis_d() + 
      xlab(label = x_label) +
      ggtitle(label = title)
    plot
  }
  x_plot <<- ggplot_build(plot = plot)$data[[1]]$x
  y_plot <<- ggplot_build(plot = plot)$data[[1]]$y
}
plot_power(data = power, x = x, group = "Method", title = "Delta = 0.1", 
           x_label = "n", smooth = TRUE)
power_plot <- stack(power) %>% 
  rename(Power = values, Method = ind) %>% 
  mutate(
    Group = rep("duration", times = length(Power)), 
    Group = if_else(
      condition = str_detect(string = Method, pattern = 'Logit|Chi'), 
      true = "binary", 
      false = Group
    ), 
    Group = if_else(
      condition = str_detect(string = Method, pattern = 'Negative.Binomial'), 
      true = "count", 
      false = Group
    ), x = rep(x, times = ncol(power))
  )
ggplot(data = power_plot) + 
  stat_smooth(mapping = aes(x = x, y = Power, group = Method, color = Method), 
              method = "loess", formula = y ~ x, level = 0) + 
  geom_line(mapping = aes(x = x, y = 0.8)) + 
  scale_color_viridis_d() + 
  xlab(label = "n") + 
  ggtitle(label = "Delta = 0.1")
# Indifference Plot ############################################################
powerdat <- data.frame(
  Negative.Binomial = c(1700, 454, 206, 115, 55), 
  Log.Rank = c(2454, 758, 343, 197, 91),
  Cox = c(2437, 755, 340, 194, 90), 
  Logit = c(4239, 1291, 576, 314, 138),
  Chi2 = c(4408, 1395, 614, 346, 156)
)

as_tibble(stack(powerdat)) %>% 
  rename(n = values, Method = ind) %>%
  mutate(delta = rep(c(0.05, 0.1, 0.15, 0.2, 0.3), times = ncol(powerdat))) %>% 
  ggplot() + 
  geom_line(mapping = aes(x = delta, y = n, group = Method, color = Method)) + 
  scale_color_viridis_d() +
  ggtitle(label = "Power = 0.8")
ggsave(filename = "power_08.png", device = "png", dpi = "retina", width = 20, 
       height = 12, units = "cm", path = "plots")

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
