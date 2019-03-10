# Simulation ###################################################################
# Packages
library(MASS)
library(tidyverse)
library(survival)
library(parallel)
library(GGally)

# Function to replicate data ###################################################
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
    time_study_log = log(time_study)
  )
}

# Function for simulation and regression of multiple datasets ##################
regress <- function(x, n, delta) {
  
  # Copy of "replicate_data", since it has to work with parallel
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
      time_study_log = log(time_study)
    )
  }
  
  regress <- function(data) {
    # negativ binomial model
    neg_bin_summary <- summary(object = neg_bin <- glm.nb(
      formula = data$seizures_treatment ~ 
        data$treatment + data$seizures_baseline + offset(data$time_study_log)
    ))
    # neg bin seizures_baseline_log
    neg_bin_summary2 <- summary(object = neg_bin2 <- glm.nb(
      formula = data$seizures_treatment ~ data$treatment + 
        data$seizures_baseline_log + offset(data$time_study_log)
    ))
    
    # survival
    survtime <- Surv(time = data$time_baseline, event = data$censor)
    # log rank
    surv_diff <- survdiff(formula = survtime ~ data$treatment)
    # cox
    cox_summary <- summary(object = cox <- coxph(
      formula =  survtime ~ data$treatment + data$seizures_baseline
    ))
    # cox seizures_baseline_log
    cox_summary2 <- summary(object = cox2 <- coxph(
      formula = survtime ~ data$treatment + data$seizures_baseline_log
    ))
    # cox only treatment
    cox_summary3 <- summary(object = cox3 <- coxph(
      formula = survtime ~ data$treatment
    ))
    
    # binary
    # logit
    logit_summary <- summary(object = logit <- glm(
      formula = data$response ~ data$treatment + data$seizures_baseline,
      family = binomial(link = "logit")
    ))
    # logit seizures_baseline_log
    logit_summary2 <- summary(object = logit2 <- glm(
      formula = data$response ~ data$treatment + data$seizures_baseline_log,
      family = binomial(link = "logit")
    ))
    # chi square
    chisq <- chisq.test(x = data$treatment, y = data$response)
    
    tibble(
      # negativ binomial model
      neg_bin_coefficient_treatment = exp(neg_bin$coefficients[2]),
      neg_bin_coefficient_seizures_baseline = exp(neg_bin$coefficients[3]),
      neg_bin_p_value_treatment = neg_bin_summary$coefficients[2, 4],
      neg_bin_p_value_seizures_baseline = neg_bin_summary$coefficients[3, 4],
      neg_bin_conf_IRR_low = NA,
      neg_bin_conf_IRR_up = NA,
      neg_bin_se_treatment = neg_bin_summary$coefficients[2,2],
      # neg bin seizures_baseline_log
      neg_bin_log_coefficient_treatment = exp(neg_bin2$coefficients[2]),
      neg_bin_log_coefficient_seizures_baseline = exp(neg_bin2$coefficients[3]),
      neg_bin_log_p_value_treatment = neg_bin_summary2$coefficients[2, 4],
      neg_bin_log_p_value_seizures_baseline = 
        neg_bin_summary2$coefficients[3, 4],
      neg_bin_log_conf_IRR_low = NA,
      neg_bin_log_conf_IRR_up = NA,
      neg_bin_log_se_treatment = neg_bin_summary2$coefficients[2,2],
      
      # survival
      # log rank
      logrank_p_value = 1 - pchisq(surv_diff$chisq, df = 1),
      # cox
      cox_p_value_treatment = cox_summary$coefficients[1, 5], 
      cox_p_value_seizures_baseline = cox_summary$coefficients[2,5], 
      cox_HR_treatment = cox_summary$coefficients[1,2],  #hazard ratio
      cox_HR_seizures_baseline = cox_summary$coefficients[2, 2],
      cox_conf_treatment_low = NA,
      cox_conf_treatment_up = NA,
      cox_se_treatment = cox_summary$coefficients[1,3],
      cox_wald = cox$wald.test,
      # cox seizures_baseline_log
      cox_p_value_log_treatment = cox_summary2$coefficients[1, 5], 
      cox_p_value_log_seizures_baseline = cox_summary2$coefficients[2,5], 
      cox_HR_log_treatment = cox_summary2$coefficients[1,2],       
      cox_HR_log_seizures_baseline = cox_summary2$coefficients[2, 2],
      cox_log_conf_treatment_low = NA,
      cox_log_conf_treatment_up = NA,
      cox_log_se_treatment = cox_summary2$coefficients[1,3],
      cox_log_wald = cox2$wald.test,
      # cox only treatment
      cox3_p_value_treatment = cox_summary3$coefficients[1, 5],
      cox3_HR_treatment = cox_summary3$coefficients[1, 2],
      cox3_conf_treatment_low = NA,
      cox3_conf_treatment_up = NA,
      cox3_se_treatment = cox_summary3$coefficients[1,3],
      cox3_wald = cox3$wald.test,
      
      # binary
      # logit
      logit_p_value_treatment = logit_summary$coefficients[2, 4],
      logit_p_value_seizures_baseline = logit_summary$coefficients[3, 4],
      oddsratio_treatment = exp(logit$coefficients[2]),
      oddsratio_seizures_baseline = exp(logit$coefficients[3]),
      conf_treatment_low = NA,
      conf_treatment_up = NA,
      logit_se_treatment = logit_summary$coefficients[2, 2],
      # logit seizures_baseline_log
      logit_p_value_log_treatment = logit_summary2$coefficients[2, 4],
      logit_p_value_log_seizures_baseline = logit_summary2$coefficients[3, 4],
      oddsratio_log_treatment = exp(logit2$coefficients[2]),
      oddsratio_log_seizures_baseline = exp(logit2$coefficients[3]),
      conf_log_treatment_low = NA,
      conf_log_treatment_up = NA,
      logit_log_se_treatment = logit_summary2$coefficients[2, 2],
      # chi square
      chi_square_p_value = chisq$p.value
    )
  }
  data <- replicate_data(n = n, delta = delta)
  bind_cols(tibble(delta = delta, n = n), regress(data = data))
}

# Function to simulate different values of delta and n #########################
simulate_data <- function(
  number_datasets = 10000, # Number of Dataframes 
  n = 200,                 # Number of observations. If this should be variied, 
                           # it has to be a vector
  delta = 0.13,            # If observations are variied, this has to be a 
                           # constant
  name                     # Name of data frames (e.g. "n_200_delta_")
) {
  if (length(n) > 1) {
    results <- lapply(X = n, FUN = function(n) {
      bind_rows(parLapply(cl = cl, X = 1:number_datasets, fun = regress, 
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
      bind_rows(parLapply(cl = cl, X = 1:number_datasets, fun = regress, 
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

# Without parralel computation for testing #####################################
set.seed(42)
result <- bind_rows(lapply(X = 1:10, FUN = regress, n = 200, delta = 0.13))

# parallel computing for faster computation ####################################
cl <- makeCluster(spec = detectCores())
clusterEvalQ(cl = cl, expr = lapply(X = c("MASS", "tidyverse", "survival"),
                                    FUN = require, character.only = TRUE))
set.seed(seed = 42)
system.time(results <- bind_rows(
  parLapply(cl = cl, X = 1:10, fun = regress, n = 200, delta = 0.13)
))
stopCluster(cl = cl)
# Visualisation of one dataset #################################################
set.seed(seed = 18)
dataset <- replicate_data(n = 200, delta = 0.13)

# Factors for visualisation
dataset <- dataset %>% 
  mutate(
    treatment = as_factor(treatment), 
    drop_out = as_factor(drop_out), 
    censor = as_factor(censor), 
    response = as_factor(response)
  ) %>%
  select(
    -lambda_baseline, 
    -lambda_treatment, 
    -seizures_baseline_log, 
    -time_study_log
  )

dataset
summary(object = dataset)

ggpairs(
  data = dataset,
  mapping = aes(colour = 1), 
  title = "Replicated dataset",
  upper = list(continuous = "points", combo = "box", discrete = "facetbar"),
  lower = list(continuous = "points", combo = "box", discrete = "facetbar")
)
# ggsave(filename = "replicated_dataset.svg", path = "plots", scale = 2)

# Different Deltas and Ns ######################################################
x <- seq(from = 0, to = 0.35, by = 0.05)
x <- seq(from = 100, to = 450, by = 50)

cl <- makeCluster(spec = detectCores())
clusterEvalQ(cl = cl, expr = lapply(X = c("MASS", "tidyverse", "survival"),
                                    FUN = require, character.only = TRUE))


system.time(simulate_data(number_datasets = 100, delta = 0.1, n = x, 
                          name = "delta_0.1_n_"))

stopCluster(cl = cl)

# save dataframes
# save(list = ls(all.names = TRUE), file = "TestN_400.RData", 
#      envir = .GlobalEnv) 

rm(list = ls(pattern = "name"))
