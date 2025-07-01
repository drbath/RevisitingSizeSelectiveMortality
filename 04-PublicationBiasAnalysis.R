# Load packages 
require(ggplot2)
require(tidybayes)
library(brms)
library(tidyverse)
library(posterior) 
library(ggsci)
library(glue)
library(pals)
library(metafor)


#clear environment




#Get relevant column info
PubBiasDat <- readRDS('Input/effect_sizes.rds') %>%
  filter(location!='lab') %>%
  select(`Publication Year`, study_id,study_name, SMD,SMD_se, pooled_sample_size)
  


#Funnel plot Standard Error
png(file = 'Figures/FunnelPlotStandardError.png', width = 16, height = 16, units = 'cm', res = 300)
funnel(PubBiasDat$SMD, sei = PubBiasDat$SMD_se, xlab = 'Standardized mean difference', ylab = 'Standard error')
dev.off()


#Linear Regression for year 
PubBiasDat %>%
  ggplot(aes(y=SMD, x= `Publication Year`)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(
    y='Standardized Mean Difference', 
    x='Publication year'
  ) +
  theme_minimal(base_size = 18, base_family = 'Roboto')
  
ggsave('Figures/PublicationBiasYear.png', width = 16, height = 22, units = 'cm', dpi=300)







