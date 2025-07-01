#Load tidyverse and ggplot2
library(tidyverse) 
library(ggplot2)

#clear environment
rm(list = ls())

#Read raw data
raw_ES <- read_csv("Input/ES_MeanDifferences.CSV")

#filter for wild studies only 
raw_ES <- raw_ES %>%
  filter(location!='lab')


# Creating functions to calculate effect sizes
SD_to_Var <- function(SD) {
  Var = SD^2
  return(Var)
}

SE_to_SD <- function(SE, N) {
  SD = (SE)*(sqrt(N))
  return(SD)
}

CI_to_SD <- function(upper, lower, N) {
  SD = ((upper-lower)*(sqrt(N)))/3.92
  return(SD)
}

SD_pooled <- function(n1, s1, n2, s2) {
  
  s_pooled = sqrt((((n1-1)*s1^2) +  ((n2-1)*s2^2)/((n1-1)+(n2-1))))
  return(s_pooled)
}

SMD_between <- function(x1, x2, s_pooled) {
  
  SMD = (x1-x2)/s_pooled
  return(SMD)
}
  
SMD_between_SE <- function(n1, n2, SMD_between) {
  
  SMD_SE = sqrt((((n1+n2)/(n1*n2))+((SMD_between^2)/(2*(n1+n2)))))
  return(SMD_SE)
  
}
               


# Turning CI's to SD's
CI_rows <- raw_ES %>%
  dplyr::filter(surv_CI_upper!="NA") %>%
  mutate(orig_sample_size = as.numeric(orig_sample_size)) %>%
  mutate(surv_sample_size = as.numeric(surv_sample_size)) %>%
  mutate(SD_orig = CI_to_SD(upper=orig_CI_upper, lower=orig_CI_lower, N=orig_sample_size)) %>%
  mutate(SD_surv = CI_to_SD(upper=surv_CI_upper, lower=surv_CI_lower, N=surv_sample_size))
  
# Turning SE's to SD's
SE_rows <- raw_ES %>%
  dplyr::filter(surv_SE!="NA") %>%
  mutate(orig_sample_size = as.numeric(orig_sample_size)) %>%
  mutate(surv_sample_size = as.numeric(surv_sample_size)) %>%
  mutate(SD_orig = SE_to_SD(orig_SE, orig_sample_size)) %>%
  mutate(SD_surv = SE_to_SD(surv_SE, surv_sample_size))
  
# Creating column to match for values that are already SD's
SD_rows <- raw_ES %>%
  dplyr::filter(orig_SD!="NA") %>%
  mutate(orig_sample_size = as.numeric(orig_sample_size)) %>%
  mutate(surv_sample_size = as.numeric(surv_sample_size)) %>%
  mutate(SD_orig = orig_SD) %>%
  mutate(SD_surv = surv_SD)

# Calculate SMD, and SMD SE  
  dat_es <- SD_rows %>%
    rbind.data.frame(., SE_rows, CI_rows) %>%
    mutate(pooled_SD = SD_pooled(orig_sample_size, SD_orig, surv_sample_size, SD_surv)) %>%
    mutate(SMD = SMD_between(orig_mean, surv_mean, pooled_SD), 
           SMD_se = SMD_between_SE(orig_sample_size, surv_sample_size, SMD)) %>%
    mutate(pooled_sample_size = surv_sample_size + orig_sample_size)
    
    saveRDS(dat_es, 'Input/effect_sizes.rds')
  
# Plot ES's
dat_es %>%
  ggplot(aes(factor(ES_id), SMD)) +
  geom_point(aes()) +
  geom_errorbar(aes(ymin=SMD-SMD_se, ymax=SMD+SMD_se), width=.8) +
  theme_classic() +
  geom_hline(yintercept = 0) +
  ylim(c(-1, 1)) +
  labs(x="Effect size ID", y= "Effect Size")
    
    
