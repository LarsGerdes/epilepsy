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
      treatment = rbinom(n = n, size = 1, prob = 0.5),  
     
      # time_study
      time_study = round(pmin(rexp(n = n, rate = -log(0.8) / 56), 56)),
      time_study = if_else(condition = time_study != 0, true = time_study, 
                           false = 1), 
      
      # drop_out
      drop_out = if_else(condition = time_study == 56, true = 0, false = 1),
      
      # seizures_baseline and time_baseline with individual lambda
      # lambda
      lambda_treatment = lambda_baseline / (exp(treatment * delta + 0.2) * 28)
    )
  # List with duration_times between seizures of each patient 
  duration_times <- lapply(
    X = lapply(X = output$lambda_treatment, FUN = rexp, n = 60), 
    FUN = cumsum
  )
  mutate(
    .data = output,
    # seizures_treatment
    seizures_treatment = mapply(
      FUN = function(x, y) {length(x[x <= y])},
      x = duration_times,
      y = time_study
    ),
    # time_baseline
    time_baseline = round(unlist(lapply(
      X = mapply(FUN = function(x, y) {x[1:y]}, x = duration_times, 
                 y = seizures_baseline),
      FUN = function(x) {x[length(x)]}
    ))),
    time_baseline = if_else(condition = time_baseline <= time_study, 
                            true = time_baseline, false = time_study),
    
    # censor
    censor = if_else(condition = seizures_treatment < seizures_baseline,
                     true = 0, false = 1),
    
    # response
    response = if_else(
      condition = seizures_treatment <= seizures_baseline & drop_out == 0,
      true = 1,
      false = 0
    ),
    
    # Log-transformations
    seizures_baseline_log = log(seizures_baseline),
    time_study_log = log(time_study),
    
    # Factors for visualization
    treatment = as.factor(treatment),
    drop_out = as.factor(drop_out),
    censor = as.factor(censor),
    response = as.factor(response)
  ) %>%
    
    # Deleting columns of lambdas
    select(-lambda_baseline, -lambda_treatment)
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