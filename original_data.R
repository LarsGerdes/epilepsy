# Original Dataset #############################################################
# Packages
library(tidyverse)
library(GGally)

# Data
epilepsy <- as_tibble(read.csv(file = "epilepsy.csv", sep = ";")) %>% 
  rename(
    seizures_baseline = number.of.seizures.at.baseline, 
    seizures_treatment = number.of.seizures.under.treatment, 
    time_study = time.in.study..days.,
    drop_out = drop.out,
    time_baseline = time.to.baseline.number
  ) %>% 
  mutate(
    treatment = as.factor(treatment), 
    drop_out = as.factor(drop_out), 
    censor = as.factor(if_else(condition = censor == 0, true = 1, false = 0)), 
    response = as.factor(if_else(
      condition = seizures_treatment <= seizures_baseline & drop_out == 0, 
      true = 1, 
      false = 0
    )), 
    seizures_baseline_log = log(seizures_baseline), 
    time_study_log = log(time_study)
  )
epilepsy
summary(object = epilepsy)

# Visualization
ggpairs(
  data = epilepsy, 
  columns = 1:9,
  mapping = aes(colour = 1), 
  title = "Original dataset",
  upper = list(continuous = "points", combo = "box", discrete = "facetbar"),
  lower = list(continuous = "points", combo = "box", discrete = "facetbar")
)
# ggsave(filename = "original_dataset.svg", path = "plots", scale = 2)
