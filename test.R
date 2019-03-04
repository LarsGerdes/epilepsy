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
time_baseline <- c()
for (i in 1:n) {
  count <- 0
  time <- 0
  while (time < 56) {
    count <- count + 1
    time <- time + rexp(n = 1, rate = lambda_treatment[i])
    if (count == seizures_baseline[i]) {
      time_baseline[i] <- time
    }
  }
  seizures_treatment[i] <- count - 1
  if (is.na(time_baseline[i])) {
    time_baseline[i] <- time_study[i]
  }
}
time_baseline <- round(if_else(condition = time_baseline < time_study,
                               true = time_baseline, false = time_study))
summary(seizures_treatment)
summary(round(time_baseline))

if (is.na(time_baseline[i])) {
  time_baseline[i] <- time_study[i]
}
if (time_baseline[i] > time_study[i]) {
  time_baseline[i] <- time_study[i]
}


ggplot() + 
  geom_density(mapping = aes(x = seizures_treatment, colour = 1, fill = 2)) + 
  guides(colour = FALSE, fill = FALSE)
ggplot() + 
  geom_density(mapping = aes(x = time_baseline, colour = 1, fill = 2)) + 
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

# Trials of vecorization #######################################################
seizures <- mapp
lapply(X = lapply(X = lambda_treatment, FUN = rexp, n = 56), FUN = sum)
cumsum(seizures)
summary(seizures)

sum(sapply(
  X = rexp(
    n = 56, 
    lambda_treatment[2]
  ), 
  FUN = function(S, t) {
    sum(S < t)
  }, 
  t = 56
))
sum(
  rexp(
    n = 56, 
    lambda_treatment[1]
  ) < 56
)



cumsum(sapply(X = lambda_treatment, FUN = rexp, n = 56) < 56)
n_func <- function(t, S) sapply(t, function(t) sum(S <= t))
n_func(t = 56, S = lapply(X = lambda_treatment, FUN = rexp, n = 1))
vap

sum






round(unlist(lapply(X = mapply(FUN = rexp, 
                               n = seizures_baseline, 
                               rate = lambda_treatment), 
                    FUN = sum)))
