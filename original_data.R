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
# ggsave(filename = "original_dataset_test.svg", path = "plots", scale = 2)
# Regression ###################################################################
# Logit
summary(object = logit <- glm(
  formula = response ~ treatment + seizures_baseline, 
  family = binomial(link = "logit"), 
  data = epilepsy
))
summary(object = logit_log <- glm(
  formula = response ~ treatment + seizures_baseline_log, 
  family = binomial(link = "logit"), 
  data = epilepsy
))
# AIC of logit_log is smaller at 0.01
# -> We choose the model without transformation.
# -> Only intercept is significant.
# -> Treatment has no significant effect.

# Chi^2
table(epilepsy$response, epilepsy$treatment, dnn = c("response", "treatment"))
chisq.test(x = epilepsy$treatment, y = epilepsy$response)
# -> Do not reject H0.
# -> Independence.

# Cox
cox_summary <- summary(object = cox <- coxph(
  formula =  Surv(time = time_baseline, event = as.numeric(censor)) ~ 
    treatment + seizures_baseline, 
  data = epilepsy
))
cox_summary <- summary(object = cox <- coxph(
  formula =  Surv(time = time_baseline, event = as.numeric(censor)) ~ 
    treatment + seizures_baseline_log, 
  data = epilepsy
))
exp(confint(cox))

# Log-Rank
survdiff(
  formula = Surv(time = time_baseline, event = as.numeric(censor)) ~ 
    treatment, 
  data = epilepsy 
)

# Negative Binomial
neg_bin_summary <- summary(object = neg_bin <- glm.nb(
  formula = seizures_treatment ~ 
    treatment + seizures_baseline + offset(time_study_log), 
  data = epilepsy
))
neg_bin_summary <- summary(object = neg_bin <- glm.nb(
  formula = seizures_treatment ~ 
    treatment + seizures_baseline_log + offset(time_study_log), 
  data = epilepsy
))
# AIC of are identical
# -> We choose the model without transformation.
# -> All parameters are significant at 0.05.