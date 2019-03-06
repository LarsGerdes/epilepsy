# Test #########################################################################
n <- 400
delta <- 0.13

# seizures_baseline with individual lambda and more than three seizures. #######
tib <- tibble(
  lambda_baseline = rgamma(n = n * 1.5, shape = 5, scale = 2),
  seizures_baseline = sapply(X = lambda_baseline, FUN = rpois, n = 1)
) %>% 
  filter(seizures_baseline > 3) %>% 
  slice(1:n)
lambda_baseline <- tib$lambda_baseline
seizures_baseline <- tib$seizures_baseline
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
ggplot() + 
  geom_density(mapping = aes(x = time_baseline, colour = 1, fill = 2)) + 
  guides(colour = FALSE, fill = FALSE)
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

ggplot() + 
  geom_point(mapping = aes(x = time_study, y = seizures_treatment, 
                           colour = 1)) + 
  guides(colour = FALSE)

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
x <- seq(from = 0, to = 0.75, by = 0.05)
power_0.8 <- function(power = 0.8, x = x, data = power) {
  
  index <- 1:ncol(data)
  
  # Vectors are split into values above and below 0.8
  intersect <- apply(X = diff(data > power) != 0, MARGIN = 2, FUN = which)
  
  # Point below the split.
  y1 <- unlist(sapply(X = index, FUN = function(i) {data[intersect[i], i]}))
  # Point above the split.
  y2 <- unlist(sapply(X = index, FUN = function(i) {data[intersect[i] + 1, i]}))
  
  # Slope of intersection
  slope <- (y2 - y1) / (x[intersect + 1] - x[intersect])
  
  # Intersection
  data.frame(x_values = x[intersect] + ((0.8 - y1) / (slope - 0)))
}
power_0.8(power = 0.8, x = x, data = power)
cat("Power =", power)
index <- 1:ncol(power)
intersect <- apply(X = diff(power > 0.8) != 0, MARGIN = 2, FUN = which)
y1 <- unlist(sapply(X = index, FUN = function(i) {power[intersect[i], i]}))
n <- 1:ncol(power)
y  <- rep(0.8, length(x))
y <- 0.8
1:n
dim(power)
dim(tibble(1:ncol(power)))
intersect <- apply(X = diff(power > 0.8) != 0, MARGIN = 2, FUN = which)
y1 <- sapply(X = 1:length(intersect), 
            FUN = function(x) {power[intersect[x], x]})
power[intersect[1], 1]
y2 <- unlist(sapply(X = 1:n, FUN = function(x) {power[intersect[x] + 1, x]}))
slope_1 <- (y2 - y1) / (x[intersect + 1] - x[intersect])
slope_2 <- rep(0, times = ncol(power))
slope_2 <- 0
point_x <- x[intersect] + 
  ((0.8 - y1) / (slope - 0))

y <- power$neg_bin_p_value_treatment
y1 <- power$neg_bin_p_value_treatment
y2  <- rep(0.8, length(x))

above <- y1 > y2
intersect <- which(diff(above) != 0)

slope_1 <- (y1[intersect + 1] - y1[intersect]) / 
  (x[intersect + 1] - x[intersect])
slope_2 <- (y2[intersect + 1] - y2[intersect]) /
  (x[intersect + 1] - x[intersect])

point_x <- x[intersect] + 
  ((y2[intersect] - y1[intersect]) / (slope_1 - slope_2))
point_y <- y1[intersect] + (slope_1 * (point_x - x[intersect]))