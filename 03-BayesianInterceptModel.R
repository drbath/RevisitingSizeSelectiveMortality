# Load packages 
require(ggplot2)
require(tidybayes)
require(bayesplot)
library(brms)
library(tidyverse)
library(posterior) 
library(ggsci)
library(glue)
library(pals)
library(ggthemes)
library(paletteer)

options(brms.backend = "cmdstanr")

#clear environment
rm(list = ls())

# Load data
dat_es <- readRDS('Input/effect_sizes.RDS') %>%
  filter(study_id!=254) # Filtering out Jensen et al. 2018

#Change larval to larvae
dat_es$LHS <- replace(dat_es$LHS, dat_es$LHS=='larval', 'larvae')


#Specify priors
priors <- c(prior(normal(0,1), class = Intercept),
            prior(cauchy(0,0.5), class = sd))

priors_info <- get_prior(SMD | se(SMD_se) ~ 1 + (1 | study_id/ES_id), data = dat_es, prior = priors)
#Fit the model 
set.seed(99)
BMod <- brm(SMD|se(SMD_se) ~1 + (1|study_id/ES_id), 
            data = dat_es, 
            prior = priors, 
            iter=10000, 
            warmup = 5000,
            chains = 4, 
            cores = 4,
            control = list(adapt_delta = 0.99, max_treedepth = 20))

#saveRDS(BMod, 'Input/FullBayesMod.RDS')
BMod <- readRDS('Input/FullBayesMod.RDS')

#### Assessing Convergence and Model Validity ####

#Checking traceplots 
post_samp <- as_draws_df(BMod)
names(post_samp)[1:3] <- c('Pooled Effect Size', 'Between study heterogeneity', 'Within study heterogeneity')

mcmc_trace(post_samp, pars = c('Pooled Effect Size', 'Between study heterogeneity', 'Within study heterogeneity')) 

ggsave('Figures/TracePlot.png', height = 4, width = 20, units = 'cm')


#Using potential scale reduction factor (Rhat) less than 1.01 to signify the convergence
summary(BMod)

#Poster predictive check figure
y_rep <- posterior_predict(BMod)

y_obs <- dat_es$SMD

ppc_dens_overlay(y_obs, y_rep[1:100, ]) +
  xlab("Cohen's d") +
  scale_color_manual(values = c('black', 'grey'), guide = 'none') +
  theme_default() 

ggsave('Figures/PosterPredict.png', height = 20, width = 16, units = 'cm')


#Group level parameter estimates
ranef(BMod)

#Grabbing population level parameter estimates
post.draws <- as_draws_df(BMod) %>%
  dplyr::select(b_Intercept, sd_study_id__Intercept, `sd_study_id:ES_id__Intercept`)

pop_effects <- as.data.frame(post.draws)

names(pop_effects) <- c("smd", "tau", "epsilon")

#Figure for population level parameter estimates
# SMD mu
ggplot(aes(x = smd), data = pop_effects) +
  geom_density(fill = "lightblue",                
               color = "lightblue", alpha = 0.7) +  
  geom_point(y = 0,                              
             x = mean(pop_effects$smd)) +
  labs(x = 'Pooled effect size',
       y = element_blank()) +
  theme_default()
ggsave('Figures/PooledEff_ProbDist.png', height = 20, width = 16, units = 'cm')


# Interpreting probablity of paramter estimates being smaller according to ECDF
smd.ecdf <- ecdf(pop_effects$smd)
smd.ecdf(-0.3) # 0% chance of being smaller then -0.3
smd.ecdf(0) # 90% chance of being smaller than 0
smd.ecdf(-0.1) # ~15% chance of being smaller than -0.1



# SMD tau
ggplot(aes(x = tau), data = pop_effects) +
  geom_density(fill = "lightgreen",               
               color = "lightgreen", alpha = 0.7) +  
  geom_point(y = 0, 
             x = mean(pop_effects$tau)) +        
  labs(x = 'Between study heterogeneity',
       y = element_blank()) +
  theme_default()

ggsave('Figures/BetweenHet_ProbDist.png', height = 20, width = 16, units = 'cm')

# SMD epsilon 
ggplot(aes(x = epsilon), data = pop_effects) +
  geom_density(fill = "violet",               
               color = "violet", alpha = 0.7) +  
  geom_point(y = 0, 
             x = mean(pop_effects$epsilon)) +        
  labs(x = "Within study heterogeneity",
       y = element_blank()) +
  theme_default()
ggsave('Figures/WithinHet_ProbDist.png', height = 20, width = 16, units = 'cm')



#### Creating forest plot 

BMod <- readRDS('Input/FullBayesMod.RDS')

#calculate each ES of each study by adding pooled ES to estimate dev. of each study
study.draws <- spread_draws(BMod, r_study_id[study_id,], b_Intercept) %>% 
  mutate(b_Intercept = r_study_id + b_Intercept) %>%
  mutate(study_id = as.character(study_id))

#calculating pooled effect of ES's
pooled.effect.draws <- spread_draws(BMod, b_Intercept) %>%
  mutate(study_id = "Pooled Effect")

#Binding study draws and pooled effects draws
forest.data <- bind_rows(study.draws, 
                         pooled.effect.draws) %>%
  ungroup() %>%
  mutate(study_id = reorder(study_id, b_Intercept))

#Create study names dataframe - author names
StudyNames <- distinct(readRDS('Input/effect_sizes.rds')[c('study_id', 'study_name')])
StudyNames$study_id <- factor(StudyNames$study_id) #Need to change to factor to combine with forest data

#Add credible intervals
forest.data.summary <- group_by(forest.data, study_id) %>%
  mean_qi(b_Intercept) %>%
  left_join(StudyNames) %>%
  mutate(study_name = replace_na(study_name, 'Pooled Effect'))

#Save forest plot data
saveRDS(forest.data.summary, 'Input/StudyEffectSizes.RDS')

color_palettte <- c(rep('#FC7D0BFF', 3), '#1170AAFF', rep('#FC7D0BFF', 7))

#Create forest plot 
forest.data.summary %>%
ggplot(aes(b_Intercept, fct_relevel(study_name, "Pooled Effect", 
                                after = Inf)), color='black') +
  geom_vline(xintercept = fixef(BMod)[1, 1], 
             color = "grey", linewidth = 1) +
  geom_vline(xintercept = fixef(BMod)[1, 3:4], 
             color = "grey", linetype = 2) +
  geom_vline(xintercept = 0, color = "black", 
             linewidth = 1) +
  geom_pointrange(aes(xmin=.lower, xmax=.upper), fatten = 8, linewidth = 1.25, color=color_palettte) +
  geom_point(color='white') +
  labs(x="Standardized Mean Difference", y=element_blank()) +
  xlim(c(-0.35,0.45)) +
  theme(
    panel.background = element_rect(fill="white", color="black"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(family = "Times New Roman", size=15),
    axis.title = element_text(family = "Times New Roman", size=18), 
    legend.position = "none"
  ) +
  geom_text(data = mutate_if(forest.data.summary, 
                             is.numeric, round, 2),
            aes(label = glue("{b_Intercept} [{.lower}, {.upper}] "), 
                x = Inf), hjust = "inward") 
  scale_colour_paletteer_d("ggthemes::Color_Blind")

ggsave('Figures/ForestPlot.png', width = 20, height = 16, units = 'cm', dpi = 300)

###### Sensitivity analysis #######
#Remvoving study 296 which has 36/76 effect sizes 

SensAnal <- dat_es %>%
  filter(study_id != 254) %>%
  filter(study_id != 296)




Sens_BMod <- brm(SMD|se(SMD_se) ~1 + (1|study_id/ES_id), 
           data = SensAnal, 
           prior = priors, 
           iter=10000, 
           warmup = 5000, 
           chains = 4,          
           cores = 4,           
           threads = threading(2),
           control = list(adapt_delta = 0.99, max_treedepth = 20))

#saveRDS(Sens_BMod, 'Input/Sens_BMod.RDS')


Sens_BMod <- readRDS('Input/Sens_BMod.RDS')

#### Assessing Convergence and Model Validity ####

#Checking traceplots 
plot(Sens_BMod)

#Using potential scale reduction factor (Rhat) less than 1.01 to signify the convergence
summary(Sens_BMod)

#Poster predictive check figure
pp_check(Sens_BMod)

#Group level parameter estimates
ranef(Sens_BMod)

#Grabbing population level parameter estimates
sens_post.draws <- as_draws_df(Sens_BMod) %>%
  dplyr::select(b_Intercept, sd_study_id__Intercept, `sd_study_id:ES_id__Intercept`)

sens_pop_effects <- as.data.frame(sens_post.draws)

names(sens_pop_effects) <- c("smd", "tau", "epsilon")

#Figure for population level parameter estimates
# SMD mu
ggplot(aes(x = smd), data = sens_pop_effects) +
  geom_density(fill = "darkblue",                
               color = "darkblue", alpha = 0.7) +  
  geom_point(y = 0,                              
             x = mean(sens_pop_effects$smd)) +
  labs(x = 'Pooled effect size',
       y = element_blank()) +
  theme_default()

ggsave('Figures/PooledEff_ProbDist_SensAnalys.png', height = 20, width = 16, units = 'cm')


Sens_smd.ecdf <- ecdf(sens_pop_effects$smd)
Sens_smd.ecdf(-0.3) # 0% chance of being smaller then -0.3
Sens_smd.ecdf(0) # 78% chance of being smaller than 0
Sens_smd.ecdf(-0.1) # 1.1% chance of being smaller than -0.1


# SMD tau
ggplot(aes(x = tau), data = sens_pop_effects) +
  geom_density(fill = "darkgreen",               
               color = "darkgreen", alpha = 0.7) +  
  geom_point(y = 0, 
             x = 0.04) +        
  labs(x = 'Between study heterogeneity',
       y = element_blank()) +
  theme_default()

ggsave('Figures/BetweenHet_ProbDist_SensAnalys.png', height = 20, width = 16, units = 'cm')

# SMD epsilon 
ggplot(aes(x = epsilon), data = sens_pop_effects) +
  geom_density(fill = "purple",               
               color = "purple", alpha = 0.7) +  
  geom_point(y = 0, 
             x = mean(pop_effects$epsilon)) +        
  labs(x = "Within study heterogeneity",
       y = element_blank()) +
  theme_default()

ggsave('Figures/WithinHet_ProbDist_SensAnalys.png', height = 20, width = 16, units = 'cm')


#### Creating forest plot 

#calculate each ES of each study by adding pooled ES to estimate dev. of each study
sens_study.draws <- spread_draws(Sens_BMod, r_study_id[study_id,], b_Intercept) %>% 
  mutate(b_Intercept = r_study_id + b_Intercept) %>%
  mutate(study_id = as.character(study_id))

#calculating pooled effect of ES's
sens_pooled.effect.draws <- spread_draws(Sens_BMod, b_Intercept) %>%
  mutate(study_id = "Pooled Effect")

#Binding study draws and pooled effects draws
sens_forest.data <- bind_rows(sens_study.draws, 
                              sens_pooled.effect.draws) %>%
  ungroup() %>%
  mutate(study_id = reorder(study_id, b_Intercept))

#Add credible intervals
sens_forest.data.summary <- group_by(sens_forest.data, study_id) %>%
  mean_qi(b_Intercept) %>%
  left_join(StudyNames) %>%
  mutate(study_name = replace_na(study_name, 'Pooled Effect'))

color_palettte2 <- c(rep('#FC7D0BFF', 4), '#1170AAFF', rep('#FC7D0BFF', 5))

#Create forest plot 
sens_forest.data.summary %>%
  ggplot(aes(b_Intercept, fct_relevel(study_name, "Pooled Effect", 
                                  after = Inf)), color='black') +
  geom_vline(xintercept = fixef(Sens_BMod)[1, 1], 
             color = "grey", linewidth = 1) +
  geom_vline(xintercept = fixef(Sens_BMod)[1, 3:4], 
             color = "grey", linetype = 2) +
  geom_vline(xintercept = 0, color = "black", 
             linewidth = 1) +
  geom_pointrange(aes(xmin=.lower, xmax=.upper), fatten = 8, linewidth = 1.25, color=color_palettte2) +
  labs(x="Standardized Mean Difference", y=element_blank()) +
  xlim(c(-0.2,0.3)) +
  theme(
    panel.background = element_rect(fill="white", color="black"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(family = "Times New Roman", size=15),
    axis.title = element_text(family = "Times New Roman", size=18), 
    legend.position = "none"
  ) +
  geom_point(color='white') +
  geom_text(data = mutate_if(sens_forest.data.summary, 
                             is.numeric, round, 2),
            aes(label = glue("{b_Intercept} [{.lower}, {.upper }] "), 
                x = Inf), hjust = "inward") 
  

ggsave('Figures/ForestPlotSensitivityAnalysis.png', width = 20, height = 16, units = 'cm')






















