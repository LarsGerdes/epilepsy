#### Replication ####
# Packages
library(tidyverse)
library(GGally)

# Replication
replicate_data <- function(n = 200, delta = 0.13){
  output <- tibble(
    
    # seizures_baseline with individual lambda and more than three seizures
    lambda_baseline = rgamma(n = n * 1.5, shape = 5, scale = 2),
    seizures_baseline = sapply(X = lambda_baseline, FUN = rpois, n = 1)
  ) %>% 
    filter(seizures_baseline > 3) %>% 
    slice(1:n) %>% 
    mutate(
      
      # treatment
      treatment = as.factor(rbinom(n = n, size = 1, prob = 0.5)),  
     
      # time_study
      time_study = round(pmin(rexp(n = n,
                                   rate = -log(0.8) / 56), 56)),
      time_study = if_else(condition = time_study != 0, true = time_study, 
                           false = 1), 
      
      # drop_out
      drop_out = as.factor(if_else(condition = time_study == 56, true = 0, 
                                   false = 1))
    )
  # seizures_baseline and time_baseline with individual lambda
  lambda_treatment <- output$lambda_baseline / 
    (exp((as.numeric(output$treatment) - 1) * delta + 0.2) * 28)
  seizures_treatment <- c()
  time_baseline <- c()
  for (i in 1:n) {
    count <- 0
    time <- 0
    while (time < 56) {
      count <- count + 1
      time <- time + rexp(n = 1, rate = lambda_treatment[i])
      if (count == output$seizures_baseline[i]) {
        time_baseline[i] <- time
      }
    }
    seizures_treatment[i] <- count - 1
    if (is.na(time_baseline[i])) {
      time_baseline[i] <- output$time_study[i]
    }
  }
  mutate(
    .data = output,
    seizures_treatment = seizures_treatment,
    time_baseline = round(if_else(condition = time_baseline < time_study,
                                  true = time_baseline, 
                                  false = time_study)),
    # censor
    censor = as.factor(if_else(
      condition = seizures_treatment < seizures_baseline,
      true = 0,
      false = 1
    )),
    
    # response
    response = as.factor(if_else(
      condition = seizures_treatment <= seizures_baseline & drop_out == 0,
      true = 1,
      false = 0
    )),
    
    # Log-transformations
    seizures_baseline_log = log(seizures_baseline),
    time_study_log = log(time_study)
  ) %>%
    
    # Deleting columns of lambdas
    select(-lambda_baseline)
}

set.seed(seed = 42)
dataset <- replicate_data(n = 200, delta = 0.13)
dataset
summary(object = dataset)

# Visualization
# pdf(file = "plots/replicated_data.pdf", width = 13)
ggpairs(
  data = dataset, 
  columns = 1:8,
  mapping = aes(colour = 1), 
  title = "Replicated data",
  upper = list(continuous = "points", combo = "box", discrete = "facetbar"),
  lower = list(continuous = "points", combo = "box", discrete = "facetbar")
)
# dev.off()