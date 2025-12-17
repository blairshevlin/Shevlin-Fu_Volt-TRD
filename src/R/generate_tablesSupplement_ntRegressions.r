# Copyright (C) 2025 Blair Shevlin <blairshevlin@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Record of Revisions
#
# Date            Programmers                         Descriptions of Change
# ====         ================                       ======================
# 2025/11/14    Blair Shevlin                          finalized code for rebuttal

rm(list = ls())

library(tidyverse)
library(fs)
library(here)
library(patchwork)
library(ggpubr)
library(RColorBrewer)
library(lmerTest)
library(GGally)
library(ggeffects)
library(readxl)
library(purrr)
library(broom)
library(caTools)
library("ggnewscale")
library(rstatix)
library(scales)
library(ggrepel)
library(effectsize)
library(emmeans)
library(broom.mixed)
library(knitr)
library(kableExtra)
library(flextable)
library(officer)
library(cowplot)


# Paths
dir = path(here())
nt_dir = dir / "data" / "nt" / "processed"
beh_dir = dir / "data" / "behavior"
clin_dir = dir / "data" / "clinical"
res_dir = dir / "results" # Updated

# Specify colors
NT_colors = data.frame(id = c("DA","5-HT"),
                       color = c("#cb181d","#2171b5"))

# Subject-level NT Data
load(nt_dir / "UG_RL_NT-Continuous_9-23-25.RData")

tasks.EST.alt = tasks.EST.alt %>%
    mutate( stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
            task = factor(task, levels = c("RL","UG")))

# IDs
ids_final = unique(ug.EST.alt$idx)

# Load clinical data
cl.df = read.csv(clin_dir / "clinical-data_deid_07-10-25.csv" ) 
cl.df.baseline = cl.df %>% filter(session == "pre stim") %>% select(idx,HDRS)

# Load behavioral data
rl.beh = read.csv(file = beh_dir / "rl-data_deid_07-10-25.csv")

ug.beh = read.csv(file = beh_dir / "ug-data_deid_07-10-25.csv")


rl.EST.Reward = rl.EST %>%
  filter(event == "Reward") %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))

ug.EST.Offer = ug.EST %>%
  filter(event == "Offer") %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim"))) 

# Regression formats
rl.EST.pr = rl.EST.Reward %>%
  mutate(block_type = factor(cond, levels = c("Mixed","Negative","Positive")),
         prev_rew_raw_f = factor(prev_rew_raw, levels = c(0,1)),
         prev_rew = ifelse(prev_rew_raw == 1,
                           ifelse(block_type=="Positive" | block_type == "Mixed",10,0),
                           ifelse(block_type=="Negative" | block_type == "Mixed",-10,0)),
         prev_rew_f = factor(prev_rew, levels = c(0,-10,10)),
         outcome_f = factor(outcome, levels = c(0,1),labels = c("Loss","Win")),
         rew = ifelse(outcome==1,
                      ifelse(block_type=="Negative",0,10),
                      ifelse(block_type=="Positive",0,-10)),
         rew_f = factor(rew, levels=c(0,-10,10)),          
         trial_z = scale(trial)[,1],
         w_trial_z = scale(trial_within_block)[,1]
  ) 

ug.EST.pr = ug.EST.Offer %>%
  mutate(offer_bin_f = factor(offer_bin, levels = c("Middle","Low","High")),
         offer_z = scale(offer)[,1],
         offer_change = offer - prev_offer_raw,
         offer_change_z = scale(offer_change)[,1],
         offer_change_f = factor(ifelse(offer_change > 0, "improve","worse"),
                                 levels=c("worse","improve")),
          trial_z = scale(trial)[,1]
  )

###########################
# Trial-level regressions #
###########################

rl_trial_data <- merge(rl.EST.pr, cl.df.baseline)  %>%
  group_by(idx) %>%
  mutate(
    # RT bins for analysis
    rt_slow = ifelse(rt > median(rt, na.rm = TRUE), 1, 0),
  ) %>%
  ungroup() %>%
  mutate(
    # Create meaningful predictors
    block_type_factor = block_type,
    log_rt = log(rt),
    rt_z = scale(rt)[,1],
    trial_z = scale(trial)[,1],

    # Create reversal indicators
    reversal_trial = factor(reversal),
    post_reversal = ifelse(reversal == 1 | 
                          lag(reversal, 1, default = 0) == 1 | 
                          lag(reversal, 2, default = 0) == 1 | 
                          lag(reversal, 3, default = 0) == 1 | 
                          lag(reversal, 4, default = 0) == 1, 1, 0), # reversal trial + 4 trials after
    pre_reversal = ifelse(lead(reversal, 1, default = 0) == 1 | 
                         lead(reversal, 2, default = 0) == 1 | 
                         lead(reversal, 3, default = 0) == 1 | 
                         lead(reversal, 4, default = 0) == 1 | 
                         lead(reversal, 5, default = 0) == 1, 1, 0), # 5 trials before reversal
    
    # Error trials (assuming optimal choice exists)
    error_trial = ifelse(opt == 0, 1, 0), # Did not choose optimal option
    
    # Depression severity
    HDRS_z = scale(HDRS)[,1] ,

    # Prediction error: current reward - previous reward (simple version)
    pred_error = rew - prev_rew,
    # Categorize prediction errors
    pred_error_type = case_when(
      pred_error > 0 ~ "Positive_PE",   # Better than expected
      pred_error < 0 ~ "Negative_PE",   # Worse than expected  
      pred_error == 0 ~ "No_PE",        # As expected
      TRUE ~ "Missing"
    ),
    pred_error_type = factor(pred_error_type, levels = c("No_PE", "Negative_PE", "Positive_PE"))
  )

contrasts(rl_trial_data$stim)
# Model 1: Basic trial-level effects
rl.DA.trial.lm1 <- lmer(Oz ~ outcome_f * stim + 
                       block_type_factor * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = rl_trial_data %>% filter(nt == "DA"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "nmkbw"))
rl.SE.trial.lm1 <- lmer(Oz ~ outcome_f * stim + 
                       block_type_factor * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = rl_trial_data %>% filter(nt == "SE"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "nmkbw"))
rl.NE.trial.lm1 <- lmer(Oz ~ outcome_f * stim + 
                       block_type_factor * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = rl_trial_data %>% filter(nt == "NE"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "nmkbw"))

#Model 2: Error and reversal effects
rl.DA.trial.lm2 <- lmer(Oz ~ error_trial * stim + 
               #        reversal_trial * stim + 
                       post_reversal * stim +
                       rt_slow * stim +
                       trial_z * stim +
                       (1 + stim | idx),
                     data = rl_trial_data %>% filter(nt == "DA"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "nmkbw"))
rl.SE.trial.lm2 <- lmer(Oz ~ error_trial * stim + 
       #                reversal_trial * stim + 
                       post_reversal * stim +
                       rt_slow * stim +
                       trial_z * stim +
                       (1 + stim | idx),
                     data = rl_trial_data %>% filter(nt == "SE"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "nmkbw"))
rl.NE.trial.lm2 <- lmer(Oz ~ error_trial * stim + 
               #        reversal_trial * stim + 
                       post_reversal * stim +
                       rt_slow * stim +
                       trial_z * stim +
                       (1 + stim | idx),
                     data = rl_trial_data %>% filter(nt == "NE"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))

#Model 3: Current and previous reward effects
rl.DA.trial.lm3 <- lmer(Oz ~ rew_f * stim + 
                    prev_rew_f * stim + 
                    rew_f * prev_rew_f * stim +  # interaction between current and previous reward
                    trial_z * stim +
                    (1 + stim | idx),
                  data = rl_trial_data %>% filter(nt == "DA", !is.na(prev_rew)),
                  REML = FALSE,
                  control = lmerControl(optimizer = "nmkbw"))
rl.SE.trial.lm3 <- lmer(Oz ~ rew_f * stim + 
                    prev_rew_f * stim + 
                    rew_f * prev_rew * stim +  # interaction between current and previous reward
                    trial_z * stim +
                    (1 + stim | idx),
                  data = rl_trial_data %>% filter(nt == "SE", !is.na(prev_rew)),
                  REML = FALSE,
                  control = lmerControl(optimizer = "nmkbw"))
rl.NE.trial.lm3 <- lmer(Oz ~ rew_f * stim + 
                    prev_rew_f * stim + 
                    rew_f * prev_rew_f * stim +  # interaction between current and previous reward
                    trial_z * stim +
                    (1 + stim | idx),
                  data = rl_trial_data %>% filter(nt == "NE", !is.na(prev_rew)),
                  REML = FALSE,
                  control = lmerControl(optimizer = "nmkbw"))

# Model 4: Prediction error effects (explicitly modeled)
rl.DA.trial.lm4 <- lmer(Oz ~ pred_error_type * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = rl_trial_data %>% filter(nt == "DA", !is.na(prev_rew)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "nmkbw"))
rl.SE.trial.lm4 <- lmer(Oz ~ pred_error_type * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = rl_trial_data %>% filter(nt == "SE", !is.na(prev_rew)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "nmkbw"))
rl.NE.trial.lm4 <- lmer(Oz ~ pred_error_type * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = rl_trial_data %>% filter(nt == "NE", !is.na(prev_rew)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "nmkbw"))

summary(rl.DA.trial.lm1)
# Only Stim, Trial x Stim
summary(rl.SE.trial.lm1)
# Block Type, Stim x Block Type
summary(rl.NE.trial.lm1)
# RT, Trial, Stim x Trial

summary(rl.DA.trial.lm2)
# Stim, Trial, Stim x Trial
summary(rl.SE.trial.lm2)
# Trial, Error x Stim
summary(rl.NE.trial.lm2)
# Trial, Stim x Trial

summary(rl.DA.trial.lm3)
# Stim, Trial, rew x prev_rew, stim x trial, rew x stim x prev_rew
summary(rl.SE.trial.lm3)
# Trial, rew, prev_rew, Stim x prev_reward
summary(rl.NE.trial.lm3)
# Trial, Trial x Stim, rew x stim x prev_rew

summary(rl.DA.trial.lm4)
# Stim, Trial, trial, Stim x Trial
summary(rl.SE.trial.lm4)
# Trial
summary(rl.NE.trial.lm4)
# Trial, Stim x Trial


# Test specific contrasts

# Win vs Loss effects by stimulation condition
rl_da_outcome_effects <- emmeans(rl.DA.trial.lm1, ~ outcome_f * stim)
rl_se_outcome_effects <- emmeans(rl.SE.trial.lm1, ~ outcome_f * stim)
rl_ne_outcome_effects <- emmeans(rl.NE.trial.lm1, ~ outcome_f * stim)

pairs(rl_da_outcome_effects, by = "stim") # NS
pairs(rl_se_outcome_effects, by = "stim") # NS
pairs(rl_ne_outcome_effects, by = "stim") # NS

# Block type effects
rl_da_block_effects <- emmeans(rl.DA.trial.lm1, ~ block_type_factor * stim)
rl_se_block_effects <- emmeans(rl.SE.trial.lm1, ~ block_type_factor * stim)
rl_ne_block_effects <- emmeans(rl.NE.trial.lm1, ~ block_type_factor * stim)

pairs(rl_da_block_effects, by = "stim") # Pre-Stim: Negative - Positive
pairs(rl_se_block_effects, by = "stim") # Pre-Stim: Mixed - Positive
pairs(rl_ne_block_effects, by = "stim") # Pre-Stim: Negaitve - Positive

# Slow trial effects
rl_da_rt_effects <- emmeans(rl.DA.trial.lm2, ~ rt_slow * stim)
rl_ne_rt_effects <- emmeans(rl.SE.trial.lm2, ~ rt_slow * stim)
rl_se_rt_effects <- emmeans(rl.NE.trial.lm2, ~ rt_slow * stim)

pairs(rl_da_rt_effects, by = "stim") # NS
pairs(rl_se_rt_effects, by = "stim") # NS
pairs(rl_ne_rt_effects, by = "stim") # NS

# Error trial effects
rl_da_error_effects <- emmeans(rl.DA.trial.lm2, ~ error_trial * stim)
rl_se_error_effects <- emmeans(rl.SE.trial.lm2, ~ error_trial * stim)
rl_ne_error_effects <- emmeans(rl.NE.trial.lm2, ~ error_trial * stim)

pairs(rl_da_error_effects, by = "stim") # NS
pairs(rl_se_error_effects, by = "stim") # NS
pairs(rl_ne_error_effects, by = "stim") # NS

# Reversal trial effects
rl_da_reversal_effects <- emmeans(rl.DA.trial.lm2, ~ post_reversal * stim)
rl_se_reversal_effects <- emmeans(rl.SE.trial.lm2, ~ post_reversal * stim)
rl_ne_reversal_effects <- emmeans(rl.NE.trial.lm2, ~ post_reversal * stim)

pairs(rl_da_reversal_effects, by = "stim") # NS
pairs(rl_se_reversal_effects, by = "stim") # NS
pairs(rl_ne_reversal_effects, by = "stim") # NS

# Reward effects
rl_da_rew_effects <- emmeans(rl.DA.trial.lm3, ~ prev_rew_f * rew_f * stim)
rl_se_rew_effects <- emmeans(rl.SE.trial.lm3, ~ prev_rew_f * rew_f * stim)
rl_ne_rew_effects <- emmeans(rl.NE.trial.lm3, ~ prev_rew_f * rew_f * stim)

pairs(rl_da_rew_effects, by = c("prev_rew_f","stim"))
# Pre-Stim, prev_rew = 0
# 0 vs -10, 0 vs 10
pairs(rl_se_rew_effects, by = c("prev_rew_f","stim")) # NS
pairs(rl_ne_rew_effects, by = c("prev_rew_f","stim")) # NS
# Post-Stim, prev_rew = 10
# -10 vs 10

# Prediction error effects
rl_da_pe_effects <- emmeans(rl.DA.trial.lm4, ~ pred_error_type * stim)
rl_se_pe_effects <- emmeans(rl.SE.trial.lm4, ~ pred_error_type * stim)
rl_ne_pe_effects <- emmeans(rl.NE.trial.lm4, ~ pred_error_type * stim)

pairs(rl_da_pe_effects, by = "stim") # NS
pairs(rl_se_pe_effects, by = "stim") # NS
pairs(rl_ne_pe_effects, by = "stim") # NS

# RT correlations with neurotransmitters
rl_rt_correlation_data <- rl_trial_data %>%
  group_by(idx, nt, stim) %>%
  summarise(
    rt_nt_cor = cor(rt, Oz, use = "complete.obs"),
    .groups = "drop"
  )
rl_rt_cor_test <- rl_rt_correlation_data %>%
  pivot_wider(names_from = stim, values_from = rt_nt_cor) %>%
  mutate(cor_diff = `Post-Stim` - `Pre-Stim`)

# RT-Neurotransmitter Correlation Changes
rl_rt_cor_test %>% 
  group_by(nt) %>% 
  summarise(
    mean_cor_change = mean(cor_diff, na.rm = TRUE),
    se_cor_change = sd(cor_diff, na.rm = TRUE) / sqrt(n()),
    t_stat = mean_cor_change / se_cor_change,
    p_value = 2 * (1 - pt(abs(t_stat), df = n() - 1))
  ) # NS


## UG ##
########
#ug.EST.pr %>% head()

ug_trial_data <- merge(ug.EST.pr, cl.df.baseline)  %>%
  group_by(idx) %>%
  mutate(
    # RT bins for analysis
    rt_slow = ifelse(rt > median(rt, na.rm = TRUE), 1, 0),
  ) %>%
  mutate(
    # Create meaningful predictors
    outcome_binary = 1*(rej==0),  # Accept = 1, Reject = 0
    offer_bin_f = factor(offer_bin_f, levels = c("Middle","Low","High")),
    prev_offer_bin_f = ifelse(prev_offer_raw < 4, "Low", ifelse(prev_offer_raw > 7, "High", "Medium")),
    prev_offer_bin_f = factor(prev_offer_bin_f, levels = c("Medium","Low","High")),
    log_rt = log(rt),
    rt_z = scale(rt)[,1],
    trial_z = scale(trial)[,1],
    stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
    
    # Accept trials
    accept_trial = ifelse(outcome_binary == 1, 1, 0), # Accepted

    # Previous offer and outcome
    offer_change = offer - prev_offer_raw,

    # Categorize prediction errors
    offer_change_type = case_when(
      offer_change > 0 ~ "Offer_Better",   # Better than expected
      offer_change < 0 ~ "Offer_Worse",   # Worse than expected  
      offer_change == 0 ~ "Same_Offer",        # As expected
      TRUE ~ "Missing"
    ),
    offer_change_type = factor(offer_change_type, levels = c("Same_Offer", "Offer_Worse", "Offer_Better")),

    # Depression severity
    HDRS_z = scale(HDRS)[,1] 
  )

# Model 1: Basic trial-level effects
ug.DA.trial.lm1 <- lmer(Oz ~ offer_bin_f * stim + 
                       #log_rt * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "DA"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
ug.SE.trial.lm1 <- lmer(Oz ~ offer_bin_f * stim + 
                       #log_rt * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "SE"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
ug.NE.trial.lm1 <- lmer(Oz ~ offer_bin_f * stim + 
                       #log_rt * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "NE"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))

#Model 2: Error and reversal effects
ug.DA.trial.lm2 <- lmer(Oz ~ accept_trial * stim + 
                       rt_slow * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "DA"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
ug.SE.trial.lm2 <- lmer(Oz ~ accept_trial * stim + 
                       rt_slow * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "SE"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
ug.NE.trial.lm2 <- lmer(Oz ~ accept_trial * stim + 
                       rt_slow * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "NE"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))

#Model 3: Current and previous reward effects
ug.DA.trial.lm3 <- lmer(Oz ~ offer_bin_f * stim + 
                    prev_offer_bin_f * stim + 
                    offer_bin_f * prev_offer_bin_f * stim +  # interaction between current and previous reward
                    trial_z * stim +
                    (1 + stim | idx),
                  data = ug_trial_data %>% filter(nt == "DA", !is.na(prev_offer_raw)),
                  REML = FALSE,
                  control = lmerControl(optimizer = "bobyqa"))
ug.SE.trial.lm3 <- lmer(Oz ~ offer_bin_f * stim + 
                    prev_offer_bin_f * stim + 
                    offer_bin_f * prev_offer_bin_f * stim +  # interaction between current and previous reward
                    trial_z * stim +
                    (1 + stim | idx),
                  data = ug_trial_data %>% filter(nt == "SE", !is.na(prev_offer_raw)),
                  REML = FALSE,
                  control = lmerControl(optimizer = "bobyqa"))
ug.NE.trial.lm3 <- lmer(Oz ~ offer_bin_f * stim + 
                    prev_offer_bin_f * stim + 
                    offer_bin_f * prev_offer_bin_f * stim +  # interaction between current and previous reward
                    trial_z * stim +
                    (1 + stim | idx),
                  data = ug_trial_data %>% filter(nt == "NE", !is.na(prev_offer_raw)),
                  REML = FALSE,
                  control = lmerControl(optimizer = "nmkbw"))

# Model 4: Prediction error/offer change effects (explicitly modeled)
ug.DA.trial.lm4 <- lmer(Oz ~ offer_change_type * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "DA", !is.na(prev_offer_raw)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
ug.SE.trial.lm4 <- lmer(Oz ~ offer_change_type * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "SE", !is.na(prev_offer_raw)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
ug.NE.trial.lm4 <- lmer(Oz ~ offer_change_type * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "NE", !is.na(prev_offer_raw)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "nmkbw"))

summary(ug.DA.trial.lm1)
# Trial, Trial x Stim
summary(ug.SE.trial.lm1)
# Stim, RT, stim x RT
summary(ug.NE.trial.lm1)
# Trial, Stim x Trial

summary(ug.DA.trial.lm2)
# Trial, Stim x Trial
summary(ug.SE.trial.lm2)
# Stim
summary(ug.NE.trial.lm2)
# Trial, Stim x Trial

summary(ug.DA.trial.lm3)
# Stim, Trial, Stim x Trial
summary(ug.SE.trial.lm3)
# Bunch
summary(ug.NE.trial.lm3)
# Trial, Trial x Stim

summary(ug.DA.trial.lm4)
# Trial, Stim x Trial
summary(ug.SE.trial.lm4)
# Stim
summary(ug.NE.trial.lm4)
# Trial, Stim x Trial


# Test specific contrasts

# Offer Bin by stimulation condition
ug_da_outcome_effects <- emmeans(ug.DA.trial.lm1, ~ offer_bin_f* stim)
ug_se_outcome_effects <- emmeans(ug.SE.trial.lm1, ~ offer_bin_f * stim)
ug_ne_outcome_effects <- emmeans(ug.NE.trial.lm1, ~ offer_bin_f * stim)

pairs(ug_da_outcome_effects, by = "stim") # Pre-Stim, Middle - Low
pairs(ug_se_outcome_effects, by = "stim") # NS
pairs(ug_ne_outcome_effects, by = "stim") # NS

# Slow trial effects
ug_da_rt_effects <- emmeans(ug.DA.trial.lm2, ~ rt_slow * stim)
ug_ne_rt_effects <- emmeans(ug.SE.trial.lm2, ~ rt_slow * stim)
ug_se_rt_effects <- emmeans(ug.NE.trial.lm2, ~ rt_slow * stim)

pairs(ug_da_rt_effects, by = "stim") # NS
pairs(ug_se_rt_effects, by = "stim") # NS
pairs(ug_ne_rt_effects, by = "stim") # NS


pairs(ug_da_rt_effects, by = "rt_slow") # NS
pairs(ug_se_rt_effects, by = "rt_slow") # NS
pairs(ug_ne_rt_effects, by = "rt_slow") # NS

# Accept trial effects
ug_da_accept_effects <- emmeans(ug.DA.trial.lm2, ~ accept_trial * stim)
ug_se_accept_effects <- emmeans(ug.SE.trial.lm2, ~ accept_trial * stim)
ug_ne_accept_effects <- emmeans(ug.NE.trial.lm2, ~ accept_trial * stim)

pairs(ug_da_accept_effects, by = "stim") # NS
pairs(ug_se_accept_effects, by = "stim") # NS
pairs(ug_ne_accept_effects, by = "stim") # NS

# Offer effects
ug_da_offer_effects <- emmeans(ug.DA.trial.lm3, ~ prev_offer_bin_f * offer_bin_f * stim)
ug_se_offer_effects <- emmeans(ug.SE.trial.lm3, ~ prev_offer_bin_f * offer_bin_f * stim)
ug_ne_offer_effects <- emmeans(ug.NE.trial.lm3, ~ prev_offer_bin_f * offer_bin_f * stim)

pairs(ug_da_offer_effects, by = c("prev_offer_bin_f","stim")) # NS
pairs(ug_se_offer_effects, by = c("prev_offer_bin_f","stim")) # NS
pairs(ug_ne_offer_effects, by = c("prev_offer_bin_f","stim"))
# Pre-Stim, prev_offer = High
# Low-High
pairs(ug_ne_offer_effects, by = c("offer_bin_f","stim"))


# Offer change effects
ug_da_oc_effects <- emmeans(ug.DA.trial.lm4, ~ offer_change_type * stim)
ug_se_oc_effects <- emmeans(ug.SE.trial.lm4, ~ offer_change_type * stim)
ug_ne_oc_effects <- emmeans(ug.NE.trial.lm4, ~ offer_change_type * stim)

pairs(ug_da_oc_effects, by = "stim") # NS
pairs(ug_se_oc_effects, by = "stim") # NS
pairs(ug_ne_oc_effects, by = "stim") # NS

ug_correlation_data <- ug_trial_data %>%
  group_by(idx, nt, stim) %>%
  summarise(
    rt_nt_cor = cor(rt, Oz, use = "complete.obs"),
    .groups = "drop"
  )
ug_rt_cor_test <- ug_correlation_data %>%
  pivot_wider(names_from = stim, values_from = rt_nt_cor) %>%
  mutate(cor_diff = `Post-Stim` - `Pre-Stim`)
ug_rt_cor_test %>% 
  group_by(nt) %>% 
  summarise(
    mean_cor_change = mean(cor_diff, na.rm = TRUE),
    se_cor_change = sd(cor_diff, na.rm = TRUE) / sqrt(n()),
    t_stat = mean_cor_change / se_cor_change,
    p_value = 2 * (1 - pt(abs(t_stat), df = n() - 1))
  ) # NS


# Copy Seth's specific models

# Model 1: OZ ~ condition x choice + (1 + condition × choice | dataset)
# Model 2: OZ ~ condition × choice × order + (1 + condition × choice × order | dataset)
## Can't do because we had fixed order
# Model 3: RZ ~ value + value difference + (1 + value + value difference | dataset)
# Model 4: RZ ~ value + value difference + absolute value difference + (1 + value + value difference + absolute value difference | dataset)
# Model 5: RZ ~  1 + condition × (value + value difference) + (1 + condition × (value + value difference) | dataset)

ug.DA.trial.s1 <- lmer(Oz ~ stim * accept_trial +
                       (1 + stim  | idx),
                     data = ug_trial_data %>% filter(nt == "DA"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "nmkbw"))
ug.SE.trial.s1 <- lmer(Oz ~ stim * accept_trial +
                       (1 + stim  | idx),
                     data = ug_trial_data %>% filter(nt == "SE"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "nmkbw"))
ug.NE.trial.s1 <- lmer(Oz ~ stim * accept_trial +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "NE"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "nmkbw"))

ug.DA.trial.s3 <- lmer(Rz ~ offer + offer_change +
                       (1 + offer + offer_change | idx),
                     data = ug_trial_data %>% filter(nt == "DA",!is.na(prev_offer_raw)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "nmkbw"))
ug.SE.trial.s3 <- lmer(Rz ~ offer + offer_change +
                       (1 + offer + offer_change | idx),
                     data = ug_trial_data %>% filter(nt == "SE",!is.na(prev_offer_raw)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
ug.NE.trial.s3 <- lmer(Rz ~ offer + offer_change +
                       (1 + offer + offer_change | idx),
                     data = ug_trial_data %>% filter(nt == "NE",!is.na(prev_offer_raw)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))

ug.DA.trial.s3.stim <- lmer(Rz ~ offer * stim + offer_change * stim +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "DA",!is.na(prev_offer_raw)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "nmkbw"))
ug.SE.trial.s3.stim <- lmer(Rz ~ offer * stim + offer_change * stim +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "SE",!is.na(prev_offer_raw)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
ug.NE.trial.s3.stim <- lmer(Rz ~ offer * stim + offer_change * stim +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "NE",!is.na(prev_offer_raw)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
# Use the residuals for N-M3*
ug_trial_data_n3_star <- ug_trial_data %>%
  filter(!is.na(prev_offer_raw)) %>%
  mutate(offer_change_resid = residuals(lm(offer_change ~ offer)))
ug.DA.trial.s3_star <- lmer(Rz ~ offer + offer_change_resid +
                       (1 + offer + offer_change_resid | idx),
                            data = ug_trial_data_n3_star %>% filter(nt == "DA"),
                            REML = FALSE,
                            control = lmerControl(optimizer = "nmkbw"))
ug.SE.trial.s3_star <- lmer(Rz ~ offer + offer_change_resid +
                       (1 + offer + offer_change_resid | idx),
                            data = ug_trial_data_n3_star %>% filter(nt == "SE"),
                            REML = FALSE,
                            control = lmerControl(optimizer = "nmkbw"))
ug.NE.trial.s3_star <- lmer(Rz ~ offer + offer_change_resid +
                       (1 + offer + offer_change_resid | idx),
                            data = ug_trial_data_n3_star %>% filter(nt == "NE"),
                            REML = FALSE,
                            control = lmerControl(optimizer = "bobyqa"))

ug.DA.trial.s3_star.stim <- lmer(Rz ~ offer * stim + offer_change_resid * stim +
                       (1 + stim | idx),
                            data = ug_trial_data_n3_star %>% filter(nt == "DA"),
                            REML = FALSE,
                            control = lmerControl(optimizer = "nmkbw"))
ug.SE.trial.s3_star.stim <- lmer(Rz ~ offer * stim + offer_change_resid * stim +
                       (1 + stim | idx),
                            data = ug_trial_data_n3_star %>% filter(nt == "SE"),
                            REML = FALSE,
                            control = lmerControl(optimizer = "nmkbw"))
ug.NE.trial.s3_star.stim <- lmer(Rz ~ offer * stim + offer_change_resid * stim +
                       (1 + stim | idx),
                            data = ug_trial_data_n3_star %>% filter(nt == "NE"),
                            REML = FALSE,
                            control = lmerControl(optimizer = "bobyqa"))

ug.DA.trial.s4 <- lmer(Rz ~ offer +  offer_change + abs(offer_change) +
                       (1 + offer +  offer_change + abs(offer_change)  | idx),
                     data = ug_trial_data %>% filter(nt == "DA",!is.na(prev_offer_raw)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "nmkbw"))
ug.SE.trial.s4<- lmer(Rz ~ offer +  offer_change + abs(offer_change) +
                       (1 + offer +  offer_change + abs(offer_change)  | idx),
                     data = ug_trial_data %>% filter(nt == "SE",!is.na(prev_offer_raw)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
ug.NE.trial.s4 <- lmer(Rz ~ offer +  offer_change + abs(offer_change) +
                       (1 + offer +  offer_change + abs(offer_change)  | idx),
                     data = ug_trial_data %>% filter(nt == "NE",!is.na(prev_offer_raw)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))

ug.DA.trial.s4.stim <- lmer(Rz ~ stim * offer + stim * offer_change + stim * abs(offer_change) +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "DA",!is.na(prev_offer_raw)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "nmkbw"))
ug.SE.trial.s4.stim <- lmer(Rz ~ stim * offer + stim * offer_change + stim * abs(offer_change) +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "SE",!is.na(prev_offer_raw)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
ug.NE.trial.s4.stim <- lmer(Rz ~ stim * offer + stim * offer_change + stim * abs(offer_change) +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "NE",!is.na(prev_offer_raw)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))

summary(ug.DA.trial.s1) # NS
summary(ug.SE.trial.s1) # Stim only
summary(ug.NE.trial.s1) # NS

ug.DA.trial.s1 %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value, method = "BH") ) # Benjamini-Hochberg FDR)
ug.SE.trial.s1 %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value, method = "BH") )  
ug.NE.trial.s1 %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value, method = "BH") )

summary(ug.DA.trial.s3) # ns
summary(ug.SE.trial.s3) # ns
summary(ug.NE.trial.s3) # ns

ug.DA.trial.s3 %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value, method = "BH") )
ug.SE.trial.s3 %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value, method = "BH") )  
ug.NE.trial.s3 %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value, method = "BH") )

summary(ug.DA.trial.s3.stim) # ns
summary(ug.SE.trial.s3.stim) # ns
summary(ug.NE.trial.s3.stim) # ns

ug.DA.trial.s3.stim %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value, method = "BH") )
ug.SE.trial.s3.stim %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value, method = "BH") )  
ug.NE.trial.s3.stim %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value, method = "BH") )

summary(ug.DA.trial.s3_star) # ns
summary(ug.SE.trial.s3_star) # ns
summary(ug.NE.trial.s3_star) # ns

ug.DA.trial.s3_star %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value, method = "BH") )
ug.SE.trial.s3_star %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value, method = "BH") )  
ug.NE.trial.s3_star %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value, method = "BH") )

summary(ug.DA.trial.s3_star.stim) # ns
summary(ug.SE.trial.s3_star.stim) # ns
summary(ug.NE.trial.s3_star.stim) # ns

ug.DA.trial.s3_star.stim %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value, method = "BH") )
ug.SE.trial.s3_star.stim %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value, method = "BH") )  
ug.NE.trial.s3_star.stim %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value, method = "BH") )

summary(ug.DA.trial.s4) # ns
summary(ug.SE.trial.s4) # ns
summary(ug.NE.trial.s4) # ns

ug.DA.trial.s4 %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value, method = "BH") )
ug.SE.trial.s4 %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value, method = "BH") )  
ug.NE.trial.s4 %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value, method = "BH") )

summary(ug.DA.trial.s4.stim) # abs(offer_change)
summary(ug.SE.trial.s4.stim) # ns 
summary(ug.NE.trial.s4.stim) # ns

ug.DA.trial.s4.stim %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value, method = "BH") )
ug.SE.trial.s4.stim %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value, method = "BH") )  
ug.NE.trial.s4.stim %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value, method = "BH") )

# Function to extract and clean fixed effects - flexible for any model structure
extract_fixed_effects <- function(model, nt_name, task_type, model_type = "basic") {
  tidy(model, effects = "fixed") %>%
    mutate(
      neurotransmitter = nt_name,
      model_type = model_type,
      task =  task_type,
      # More comprehensive term cleaning to handle various model structures
      term_clean = case_when(
        # Basic terms
        term == "(Intercept)" ~ "Intercept",
        term == "stimPost-Stim" ~ "Stimulation: Post vs Pre",
        term == "trial_z" ~ "Trial Number (z-scored)",
        term == "log_rt" ~ "Log Response Time",
        term == "stimPost-Stim:trial_z" ~ "Stim × Trial Number",
        term == "stimPost-Stim:log_rt" ~ "Stim × Log RT",
        
        # Reversal Learning specific
        term == "outcome_fWin" ~ "Outcome: Win vs Loss",
        term == "outcome_fWin:stimPost-Stim" ~ "Outcome × Stimulation",
        term == "block_type_factorNegative" ~ "Block: Negative vs Mixed",
        term == "block_type_factorPositive" ~ "Block: Positive vs Mixed", 
        term == "stimPost-Stim:block_type_factorNegative" ~ "Stim × Block: Negative",
        term == "stimPost-Stim:block_type_factorPositive" ~ "Stim × Block: Positive",
        term == "error_trial" ~ "Error Trial",
        term == "error_trial:stimPost-Stim" ~ "Error Trial × Stimulation",
        term == "reversal_trial" ~ "Reversal Trial",
        term == "stimPost-Stim:reversal_trial" ~ "Reversal Trial × Stimulation",
        term == "post_reversal" ~ "Post-Reversal Period",
        term == "stimPost-Stim:post_reversal" ~ "Post-Reversal × Stimulation",
        term == "rt_slow" ~ "Slow Response Time",
        term == "stimPost-Stim:rt_slow" ~ "Slow RT × Stimulation",
        
        # Prediction Error models
        term == "pred_error_typeNegative_PE" ~ "Negative Prediction Error",
        term == "pred_error_typePositive_PE" ~ "Positive Prediction Error",
        term == "pred_error_typeNegative_PE:stimPost-Stim" ~ "Negative PE × Stimulation",
        term == "pred_error_typePositive_PE:stimPost-Stim" ~ "Positive PE × Stimulation",
        
        # Reward models - more specific
        term == "rew_f-10" ~ "Current Reward: -10 vs 0",
        term == "rew_f10" ~ "Current Reward: 10 vs 0",
        term == "prev_rew_f-10" ~ "Previous Reward: -10 vs 0", 
        term == "prev_rew_f10" ~ "Previous Reward: 10 vs 0",
        term == "rew_f-10:stimPost-Stim" ~ "Current Reward -10 × Stimulation",
        term == "rew_f10:stimPost-Stim" ~ "Current Reward 10 × Stimulation",
        term == "stimPost-Stim:prev_rew_f-10" ~ "Previous Reward -10 × Stimulation",
        term == "stimPost-Stim:prev_rew_f10" ~ "Previous Reward 10 × Stimulation",
        term == "rew_f-10:prev_rew_f-10" ~ "Both Rewards -10",
        term == "rew_f10:prev_rew_f-10" ~ "Current +10, Previous -10",
        term == "rew_f-10:prev_rew_f10" ~ "Current -10, Previous +10",
        term == "rew_f10:prev_rew_f10" ~ "Both Rewards +10",
        term == "rew_f-10:stimPost-Stim:prev_rew_f-10" ~ "3-way: Both -10 × Stimulation",
        term == "rew_f10:stimPost-Stim:prev_rew_f-10" ~ "3-way: Curr+10/Prev-10 × Stim",
        term == "rew_f-10:stimPost-Stim:prev_rew_f10" ~ "3-way: Curr-10/Prev+10 × Stim", 
        term == "rew_f10:stimPost-Stim:prev_rew_f10" ~ "3-way: Both +10 × Stimulation",
        
        # Ultimatum Game offers - more specific
        term == "offer_bin_fLow" ~ "Offer: Low vs Middle",
        term == "offer_bin_fHigh" ~ "Offer: High vs Middle",
        term == "offer_bin_fLow:stimPost-Stim" ~ "Low Offer × Stimulation",
        term == "offer_bin_fHigh:stimPost-Stim" ~ "High Offer × Stimulation",
        term == "accept_trial" ~ "Accept vs Reject",
        term == "accept_trial:stimPost-Stim" ~ "Accept × Stimulation",
        term == "offer_z" ~ "Offer (z-scored)",
        term == "offer_z:stim1" ~ "Offer × Stimulation",
        
        # Previous offer effects - be more specific
        term == "prev_offer_bin_fLow" ~ "Previous Offer: Low vs Middle",
        term == "prev_offer_bin_fHigh" ~ "Previous Offer: High vs Middle", 
        term == "stimPost-Stim:prev_offer_bin_fLow" ~ "Previous Low Offer × Stimulation",
        term == "stimPost-Stim:prev_offer_bin_fHigh" ~ "Previous High Offer × Stimulation",
        
        # Current × Previous offer interactions
        term == "offer_bin_fHigh:prev_offer_bin_fLow" ~ "High Offer × Previous Low",
        term == "offer_bin_fHigh:prev_offer_bin_fHigh" ~ "High Offer × Previous High",
        term == "offer_bin_fLow:prev_offer_bin_fLow" ~ "Low Offer × Previous Low",
        term == "offer_bin_fLow:prev_offer_bin_fHigh" ~ "Low Offer × Previous High",
        
        # Three-way offer interactions
        term == "offer_bin_fHigh:stimPost-Stim:prev_offer_bin_fLow" ~ "3-way: High×Stim×PrevLow",
        term == "offer_bin_fHigh:stimPost-Stim:prev_offer_bin_fHigh" ~ "3-way: High×Stim×PrevHigh",
        term == "offer_bin_fLow:stimPost-Stim:prev_offer_bin_fLow" ~ "3-way: Low×Stim×PrevLow",
        term == "offer_bin_fLow:stimPost-Stim:prev_offer_bin_fHigh" ~ "3-way: Low×Stim×PrevHigh",
        
        # Offer change effects
        term == "offer_change_typeOffer_Worse" ~ "Offer Change: Worse vs Same",
        term == "offer_change_typeOffer_Better" ~ "Offer Change: Better vs Same",
        term == "offer_change_typeOffer_Worse:stimPost-Stim" ~ "Worse Offer × Stimulation",
        term == "offer_change_typeOffer_Better:stimPost-Stim" ~ "Better Offer × Stimulation",

        # Default for unmatched terms
        TRUE ~ str_replace_all(term, "_", " ") %>% str_to_title()
      ),
      # Clean up p-values and significance
      p_adjusted = p.adjust(p.value, method = "BH"),  # Benjamini-Hochberg FDR
      p_formatted = case_when(
        p.value < 0.001 ~ "< 0.001***",
        p.value < 0.01 ~ paste0(sprintf("%.3f", p.value), "**"),
        p.value < 0.05 ~ paste0(sprintf("%.3f", p.value), "*"),
        p.value < 0.10 ~ paste0(sprintf("%.3f", p.value), "†"),
        TRUE ~ sprintf("%.3f", p.value)
      ),
      p_adj_formatted = case_when(
        p_adjusted < 0.001 ~ "< 0.001***",
        p_adjusted < 0.01 ~ paste0(sprintf("%.3f", p_adjusted), "**"),
        p_adjusted < 0.05 ~ paste0(sprintf("%.3f", p_adjusted), "*"),
        p_adjusted < 0.10 ~ paste0(sprintf("%.3f", p_adjusted), "†"),
        TRUE ~ sprintf("%.3f", p_adjusted)
      ),
      # Format estimates and standard errors
      estimate_formatted = sprintf("%.3f", estimate),
      se_formatted = sprintf("%.2f", std.error),
      t_formatted = sprintf("%.2f", statistic),
      # Add effect size categories
      effect_size = case_when(
        abs(estimate) < 0.2 ~ "Small",
        abs(estimate) < 0.5 ~ "Medium",
        abs(estimate) < 0.8 ~ "Large", 
        TRUE ~ "Very Large"
      ),
      # Identify significant effects
      is_significant = grepl("\\*|†", p_formatted)
    ) %>%
    select(neurotransmitter, task, model_type, term_clean, estimate_formatted, 
           se_formatted, t_formatted, p_formatted,p_adj_formatted, effect_size, is_significant)
}

# Compile
rl.DA.lm1.res <- extract_fixed_effects(rl.DA.trial.lm1, "Dopamine", "Reversal learning","Basic effects")
rl.SE.lm1.res <- extract_fixed_effects(rl.SE.trial.lm1, "Serotonin", "Reversal learning","Basic effects")
rl.NE.lm1.res <- extract_fixed_effects(rl.NE.trial.lm1, "Norepinephrine", "Reversal learning","Basic effects")

rl.DA.lm2.res <- extract_fixed_effects(rl.DA.trial.lm2, "Dopamine", "Reversal learning", "Learning dynamics")
rl.SE.lm2.res <- extract_fixed_effects(rl.SE.trial.lm2, "Serotonin", "Reversal learning", "Learning dynamics")
rl.NE.lm2.res <- extract_fixed_effects(rl.NE.trial.lm2, "Norepinephrine", "Reversal learning", "Learning dynamics")

rl.DA.lm3.res <- extract_fixed_effects(rl.DA.trial.lm3, "Dopamine", "Reversal learning", "Reward interactions")
rl.SE.lm3.res <- extract_fixed_effects(rl.SE.trial.lm3, "Serotonin", "Reversal learning", "Reward interactions")
rl.NE.lm3.res <- extract_fixed_effects(rl.NE.trial.lm3, "Norepinephrine", "Reversal learning", "Reward interactions")

rl.DA.lm4.res <- extract_fixed_effects(rl.DA.trial.lm4, "Dopamine", "Reversal learning", "Prediction errors")
rl.SE.lm4.res <- extract_fixed_effects(rl.SE.trial.lm4, "Serotonin", "Reversal learning", "Prediction errors")
rl.NE.lm4.res <- extract_fixed_effects(rl.NE.trial.lm4, "Norepinephrine", "Reversal learning", "Prediction errors")

ug.DA.lm1.res <- extract_fixed_effects(ug.DA.trial.lm1, "Dopamine", "Ultimatum game","Basic effects")
ug.SE.lm1.res <- extract_fixed_effects(ug.SE.trial.lm1, "Serotonin", "Ultimatum game","Basic effects")
ug.NE.lm1.res <- extract_fixed_effects(ug.NE.trial.lm1, "Norepinephrine", "Ultimatum game","Basic effects")

ug.DA.lm2.res <- extract_fixed_effects(ug.DA.trial.lm2, "Dopamine", "Ultimatum game", "Context effects")
ug.SE.lm2.res <- extract_fixed_effects(ug.SE.trial.lm2, "Serotonin", "Ultimatum game", "Context effects")
ug.NE.lm2.res <- extract_fixed_effects(ug.NE.trial.lm2, "Norepinephrine", "Ultimatum game", "Context effects")

ug.DA.lm3.res <- extract_fixed_effects(ug.DA.trial.lm3, "Dopamine", "Ultimatum game", "Offer interactions")
ug.SE.lm3.res <- extract_fixed_effects(ug.SE.trial.lm3, "Serotonin", "Ultimatum game", "Offer interactions")
ug.NE.lm3.res <- extract_fixed_effects(ug.NE.trial.lm3, "Norepinephrine", "Ultimatum game", "Offer interactions")

ug.DA.lm4.res <- extract_fixed_effects(ug.DA.trial.lm4, "Dopamine", "Ultimatum game", "Offer dynamics")
ug.SE.lm4.res <- extract_fixed_effects(ug.SE.trial.lm4, "Serotonin", "Ultimatum game", "Offer dynamics")
ug.NE.lm4.res <- extract_fixed_effects(ug.NE.trial.lm4, "Norepinephrine", "Ultimatum game", "Offer dynamics")

all.lm.res <- bind_rows(ug.DA.lm1.res,ug.DA.lm2.res,ug.DA.lm3.res,ug.DA.lm4.res,
                        ug.SE.lm1.res,ug.SE.lm2.res,ug.SE.lm3.res,ug.SE.lm4.res,
                        ug.NE.lm1.res,ug.NE.lm2.res,ug.NE.lm3.res,ug.NE.lm4.res,
                        rl.DA.lm1.res,rl.DA.lm2.res,rl.DA.lm3.res,rl.DA.lm4.res,
                        rl.SE.lm1.res,rl.SE.lm2.res,rl.SE.lm3.res,rl.SE.lm4.res,
                        rl.NE.lm1.res,rl.NE.lm2.res,rl.NE.lm3.res,rl.NE.lm4.res) 


# Function to create a formatted flextable for a specific task-model combination
create_task_model_table <- function(data, task_name, model_name) {
 # data = all_results

  # Filter data
  subset_data <- data %>%
    filter(task == task_name, model_type == model_name) %>%
    select(
      Neurotransmitter = neurotransmitter,
      Effect = term_clean,
      Est. = estimate_formatted,
      SE = se_formatted,
      `t` = t_formatted,
      `p` = p_formatted,
      `p-adj` = p_adj_formatted
    )
  
  if(nrow(subset_data) == 0) {
    return(NULL)
  }
  
  # Create flextable
 ft <- subset_data %>%
  flextable() %>%
  theme_vanilla() %>%
  bold(part = "header") %>%
  bg(bg = "#f0f0f0", part = "header") %>%
  align(align = "center", part = "header") %>%
  align(j = 3:7, align = "center", part = "body") %>%
  align(j = 1:2, align = "left", part = "body") %>%
  fontsize(size = 9) %>%  # Smaller font
  font(fontname = "Times New Roman") %>%
  # Set specific column widths (in inches)
  width(j = 1, width = 1.5) %>%  # Neurotransmitter
  width(j = 2, width = 2) %>%  # Effect (longest column)
  width(j = 3, width = 0.5) %>%  # Estimate
  width(j = 4, width = 0.5) %>%  # SE
  width(j = 5, width = 0.5) %>%  # t-value
  width(j = 6, width = 0.8) %>%  # p-value
  width(j = 7, width = 0.8) %>%  # p-adj
  # Add text wrapping for long effect names
  valign(valign = "top") %>%
  padding(padding.top = 2, padding.bottom = 2)
  
  # Highlight significant results
  sig_rows <- which(grepl("\\*|†", subset_data$`p-adj`))
  if(length(sig_rows) > 0) {
    ft <- ft %>%
      bg(i = sig_rows, bg = "#ffffcc", part = "body") %>%
      bold(i = sig_rows, j = 7, part = "body")
  }
  
  return(ft)
}

# Main function to create all tables
create_all_task_model_tables <- function(all_results) {
  
  #all_results = all.lm.res
  # Get unique combinations
  combinations <- all_results %>%
    distinct(task, model_type) %>%
    arrange(task, model_type)
  
  # Create master document
  master_doc <- read_docx()
  
  # Add title
  master_doc <- master_doc %>%
    body_add_par("Trial-Level Analysis Results", style = "heading 1") %>%
    body_add_par(paste("Generated on:", Sys.Date())) %>%
    body_add_break()
  
  # Process each combination
  for(i in 1:nrow(combinations)) {
    task_name <- combinations$task[i]
    model_name <- combinations$model_type[i]
    
    cat(paste("Creating table for:", task_name, "-", model_name, "\n"))
    
    # Create table
    ft <- create_task_model_table(all_results, task_name, model_name)
    
    if(!is.null(ft)) {
      # Add to master document
      master_doc <- master_doc %>%
        body_add_par(paste(task_name, ":", model_name, "model"), 
                     style = "heading 2") %>%
        body_add_flextable(ft) %>%
        body_add_break()
      
      # Create individual document
      individual_doc <- read_docx() %>%
        body_add_par(paste(task_name, ":", model_name, "model"), 
                     style = "heading 1") %>%
        body_add_flextable(ft)
      
      # Save individual file
      clean_filename <- paste0(gsub("[^A-Za-z0-9]", "_", task_name), "_", 
                              gsub("[^A-Za-z0-9]", "_", model_name), "_table.docx")
      print(individual_doc, target = here("results",clean_filename))
    }
  }
  
  # Save master document
  print(master_doc, target = here("results","NT_analysis_tables_.docx"))
  
  cat("Tables created successfully!\n")
  cat("Files generated:\n")
  cat("- all_trial_analysis_tables.docx (master document)\n")
  for(i in 1:nrow(combinations)) {
    clean_filename <- paste0(gsub("[^A-Za-z0-9]", "_", combinations$task[i]), "_", 
                            gsub("[^A-Za-z0-9]", "_", combinations$model_type[i]), "_table.docx")
    cat(paste("-", clean_filename, "\n"))
  }
  
  return(combinations)
}

combinations <- create_all_task_model_tables(all.lm.res)

# Function to extract and clean contrast results
extract_contrast_results <- function(contrast_obj, nt_name, task_type, contrast_type) {
  # Convert emmeans pairs result to dataframe
  contrast_df <- as.data.frame(contrast_obj)
  
  # Apply multiple testing correction within this contrast set
  contrast_df <- contrast_df %>%
    mutate(
      neurotransmitter = nt_name,
      task = task_type,
      contrast_type = contrast_type,
      # Clean contrast names for better readability
      contrast_clean = case_when(
        # Block contrasts
        grepl("Mixed.*Negative", contrast) ~ "Mixed vs Negative",
        grepl("Mixed.*Positive", contrast) ~ "Mixed vs Positive", 
        grepl("Negative.*Positive", contrast) ~ "Negative vs Positive",
        
        # Outcome contrasts
        grepl("Win.*Loss|Loss.*Win", contrast) ~ "Win vs Loss",
        
        # Offer contrasts (Ultimatum Game)
        grepl("Low.*Middle", contrast) ~ "Low vs Middle",
        grepl("Middle.*High|High.*Middle", contrast) ~ "Middle vs High", 
        grepl("Low.*High|High.*Low", contrast) ~ "Low vs High",
        
        # Reward contrasts (Reversal Learning)
        grepl("\\-10.*0", contrast) ~ "-10 vs 0",
        grepl("0.*10", contrast) ~ "0 vs 10",
        grepl("\\-10.*10", contrast) ~ "-10 vs 10",
        
        # Prediction error contrasts
        grepl("Negative_PE.*No_PE", contrast) ~ "Negative PE vs No PE",
        grepl("Positive_PE.*No_PE", contrast) ~ "Positive PE vs No PE",
        grepl("Positive_PE.*Negative_PE", contrast) ~ "Positive PE vs Negative PE",
        
        # Offer change contrasts
        grepl("Worse.*Same", contrast) ~ "Worse vs Same",
        grepl("Better.*Same", contrast) ~ "Better vs Same",
        grepl("Better.*Worse", contrast) ~ "Better vs Worse",
        
        # Accept/reject contrasts
        grepl("accept.*reject|reject.*accept", contrast, ignore.case = TRUE) ~ "Accept vs Reject",
        
        # Error trial contrasts
        grepl("error.*correct|correct.*error", contrast, ignore.case = TRUE) ~ "Error vs Correct",
        
        # RT speed contrasts
        grepl("slow.*fast|fast.*slow", contrast, ignore.case = TRUE) ~ "Slow vs Fast RT",
        
        # Default - clean up the contrast name
        TRUE ~ str_replace_all(contrast, c("\\(|\\)" = "", "rew_f" = "", "_f" = "", "_" = " ")) %>% 
               str_trim() %>%
               str_replace_all("  +", " ")  # Remove multiple spaces
      ),
      # Apply Benjamini-Hochberg FDR correction within this contrast set
      p_adjusted = p.adjust(p.value, method = "BH"),
      # Format results for presentation
      estimate_formatted = sprintf("%.3f", estimate),
      se_formatted = sprintf("%.3f", SE),
      t_formatted = sprintf("%.2f", t.ratio),
      # Format original p-values
      p_formatted = case_when(
        p.value < 0.001 ~ "< 0.001***",
        p.value < 0.01 ~ paste0(sprintf("%.3f", p.value), "**"),
        p.value < 0.05 ~ paste0(sprintf("%.3f", p.value), "*"),
        p.value < 0.10 ~ paste0(sprintf("%.3f", p.value), "†"),
        TRUE ~ sprintf("%.3f", p.value)
      ),
      # Format adjusted p-values
      p_adj_formatted = case_when(
        p_adjusted < 0.001 ~ "< 0.001***",
        p_adjusted < 0.01 ~ paste0(sprintf("%.3f", p_adjusted), "**"),
        p_adjusted < 0.05 ~ paste0(sprintf("%.3f", p_adjusted), "*"),
        p_adjusted < 0.10 ~ paste0(sprintf("%.3f", p_adjusted), "†"),
        TRUE ~ sprintf("%.3f", p_adjusted)
      ),
      # Flag significant results
      is_significant = grepl("\\*|†", p_formatted),
      is_adj_significant = grepl("\\*|†", p_adj_formatted)
    ) %>%
    select(neurotransmitter, task, contrast_type, contrast_clean, 
           estimate_formatted, se_formatted, t_formatted, 
           p_formatted, p_adj_formatted, is_significant, is_adj_significant)
  
  return(contrast_df)
}

# Function to extract and format contrast results from emmeans pairs objects
extract_contrast_results <- function(contrast_obj, nt_name, task_type, contrast_type) {
  # Convert emmeans pairs result to dataframe
  contrast_df <- as.data.frame(contrast_obj)
  
  # Apply multiple testing correction within this contrast set
  contrast_df <- contrast_df %>%
    mutate(
      neurotransmitter = nt_name,
      task = task_type,
      contrast_type = contrast_type,
      # Clean contrast names for better readability
      contrast_clean = case_when(
        # Block contrasts
        grepl("Mixed.*Negative", contrast) ~ "Mixed vs Negative",
        grepl("Mixed.*Positive", contrast) ~ "Mixed vs Positive", 
        grepl("Negative.*Positive", contrast) ~ "Negative vs Positive",
        
        # Outcome contrasts
        grepl("Win.*Loss|Loss.*Win", contrast) ~ "Win vs Loss",
        
        # Offer contrasts (Ultimatum Game)
        grepl("Low.*Middle", contrast) ~ "Low vs Middle",
        grepl("Middle.*High|High.*Middle", contrast) ~ "Middle vs High", 
        grepl("Low.*High|High.*Low", contrast) ~ "Low vs High",
        
        # Reward contrasts (Reversal Learning)
        grepl("\\-10.*0", contrast) ~ "-10 vs 0",
        grepl("0.*10", contrast) ~ "0 vs 10",
        grepl("\\-10.*10", contrast) ~ "-10 vs 10",
        
        # Prediction error contrasts
        grepl("Negative_PE.*No_PE", contrast) ~ "Negative PE vs No PE",
        grepl("Positive_PE.*No_PE", contrast) ~ "Positive PE vs No PE",
        grepl("Positive_PE.*Negative_PE", contrast) ~ "Positive PE vs Negative PE",
        
        # Offer change contrasts
        grepl("Worse.*Same", contrast) ~ "Worse vs Same",
        grepl("Better.*Same", contrast) ~ "Better vs Same",
        grepl("Better.*Worse", contrast) ~ "Better vs Worse",
        
        # Accept/reject contrasts
        grepl("accept.*reject|reject.*accept", contrast, ignore.case = TRUE) ~ "Accept vs Reject",
        
        # Error trial contrasts
        grepl("error.*correct|correct.*error", contrast, ignore.case = TRUE) ~ "Error vs Correct",
        
        # RT speed contrasts
        grepl("slow.*fast|fast.*slow", contrast, ignore.case = TRUE) ~ "Slow vs Fast RT",
        
        # Default - clean up the contrast name
        TRUE ~ str_replace_all(contrast, c("\\(|\\)" = "", "rew_f" = "", "_f" = "", "_" = " ")) %>% 
               str_trim() %>%
               str_replace_all("  +", " ")  # Remove multiple spaces
      ),
      # Apply Benjamini-Hochberg FDR correction within this contrast set
      p_adjusted = p.adjust(p.value, method = "BH"),
      # Format results for presentation
      estimate_formatted = sprintf("%.3f", estimate),
      se_formatted = sprintf("%.2f", SE),
      t_formatted = sprintf("%.2f", t.ratio),
      # Format original p-values
      p_formatted = case_when(
        p.value < 0.001 ~ "< 0.001***",
        p.value < 0.01 ~ paste0(sprintf("%.3f", p.value), "**"),
        p.value < 0.05 ~ paste0(sprintf("%.3f", p.value), "*"),
        p.value < 0.10 ~ paste0(sprintf("%.3f", p.value), "†"),
        TRUE ~ sprintf("%.3f", p.value)
      ),
      # Format adjusted p-values
      p_adj_formatted = case_when(
        p_adjusted < 0.001 ~ "< 0.001***",
        p_adjusted < 0.01 ~ paste0(sprintf("%.3f", p_adjusted), "**"),
        p_adjusted < 0.05 ~ paste0(sprintf("%.3f", p_adjusted), "*"),
        p_adjusted < 0.10 ~ paste0(sprintf("%.3f", p_adjusted), "†"),
        TRUE ~ sprintf("%.3f", p_adjusted)
      ),
      # Flag significant results
      is_significant = grepl("\\*|†", p_formatted),
      is_adj_significant = grepl("\\*|†", p_adj_formatted)
    ) %>%
    select(neurotransmitter, task, contrast_type, contrast_clean, 
           estimate_formatted, se_formatted, t_formatted, 
           p_formatted, p_adj_formatted, is_significant, is_adj_significant)
  
  return(contrast_df)
}

# Function to create formatted flextable for specific task-contrast combination
create_contrast_table <- function(data, task_name, contrast_name) {
  
  # Filter data for this specific combination
  subset_data <- data %>%
    filter(task == task_name, contrast_type == contrast_name) %>%
    select(
      Neurotransmitter = neurotransmitter,
      Contrast = contrast_clean,
      Estimate = estimate_formatted,
      SE = se_formatted,
      `t` = t_formatted,
      `p` = p_formatted,
      `p-adj` = p_adj_formatted
    )
  
  # Return NULL if no data
  if(nrow(subset_data) == 0) {
    return(NULL)
  }
  
  # Create flextable with consistent formatting
  ft <- subset_data %>%
    flextable() %>%
    theme_vanilla() %>%
    bold(part = "header") %>%
    bg(bg = "#f0f0f0", part = "header") %>%
    align(align = "center", part = "header") %>%
    align(j = 3:7, align = "center", part = "body") %>%
    align(j = 1:2, align = "left", part = "body") %>%
    fontsize(size = 9) %>%
    font(fontname = "Times New Roman") %>%
    # Set specific column widths to prevent margin overflow
    width(j = 1, width = 1.5) %>%  # Neurotransmitter
    width(j = 2, width = 1.5) %>%  # Contrast (shorter than main effects)
    width(j = 3, width = 1.0) %>%  # Estimate
    width(j = 4, width = 0.6) %>%  # SE
    width(j = 5, width = 0.6) %>%  # t-value
    width(j = 6, width = 0.7) %>%  # p-value
    width(j = 7, width = 0.7) %>%  # p-adj
    valign(valign = "top") %>%
    padding(padding.top = 2, padding.bottom = 2)
  
  # Highlight uncorrected significant results (yellow)
  sig_rows <- which(grepl("\\*|†", subset_data$`p-adj`))
  if(length(sig_rows) > 0) {
    ft <- ft %>%
      bg(i = sig_rows, bg = "#ffffcc", part = "body") %>%
      bold(i = sig_rows, j = 7, part = "body")
  }
  
  # Highlight FDR-corrected significant results (light green)
  adj_sig_rows <- which(grepl("\\*|†", subset_data$`p-adj`))
  if(length(adj_sig_rows) > 0) {
    ft <- ft %>%
      bg(i = adj_sig_rows, bg = "#ccffcc", part = "body") %>%
      bold(i = adj_sig_rows, j = 7, part = "body")
  }
  
  return(ft)
}

# Main function to create all contrast tables
create_all_contrast_tables <- function(all_contrasts) {
  
  # Get unique task-contrast combinations
  combinations <- all_contrasts %>%
    distinct(task, contrast_type) %>%
    arrange(task, contrast_type)
  
  # Create master document
  master_doc <- read_docx()
  
  # Add title and header information
  master_doc <- master_doc %>%
    body_add_par("Trial-Level Contrast Analysis Results", style = "heading 1") %>%
    body_add_par(paste("Generated on:", Sys.Date())) %>%
    body_add_par("Note: Yellow highlighting = uncorrected significant, Green highlighting = FDR-corrected significant") %>%
    body_add_break()
  
  # Process each task-contrast combination
  for(i in 1:nrow(combinations)) {
    task_name <- combinations$task[i]
    contrast_name <- combinations$contrast_type[i]
    
    cat(paste("Creating contrast table for:", task_name, "-", contrast_name, "\n"))
    
    # Create table for this combination
    ft <- create_contrast_table(all_contrasts, task_name, contrast_name)
    
    if(!is.null(ft)) {
      # Add section to master document
      master_doc <- master_doc %>%
        body_add_par(paste(task_name, ":", contrast_name, "contrasts"), 
                     style = "heading 2") %>%
        body_add_flextable(ft) %>%
        body_add_break()
      
      # Create individual document for this contrast set
      individual_doc <- read_docx() %>%
        body_add_par(paste(task_name, ":", contrast_name, "contrasts"), 
                     style = "heading 1") %>%
        body_add_par("Note: Yellow highlighting = uncorrected significant, Green highlighting = FDR-corrected significant") %>%
        body_add_flextable(ft)
      
      # Save individual file with clean filename
      clean_filename <- paste0(gsub("[^A-Za-z0-9]", "_", task_name), "_", 
                              gsub("[^A-Za-z0-9]", "_", contrast_name), "_contrasts.docx")
      print(individual_doc, target = here("results",clean_filename))
    }
  }
  # Save master document with all contrasts
  print(master_doc, target = here("results","all_contrast_results.docx"))
  
  # Print summary of what was created
  cat("\nContrast tables created successfully!\n")
  cat("Files generated:\n")
  cat("- all_contrast_results.docx (master document with all contrasts)\n")
  for(i in 1:nrow(combinations)) {
    clean_filename <- paste0(gsub("[^A-Za-z0-9]", "_", combinations$task[i]), "_", 
                            gsub("[^A-Za-z0-9]", "_", combinations$contrast_type[i]), "_contrasts.docx")
    cat(paste("-", clean_filename, "\n"))
  }
  
  return(combinations)
}

# Collect all contrast results
all_contrasts <- bind_rows(
  extract_contrast_results(pairs(rl_da_block_effects, by = "stim"), "Dopamine", "Reversal learning", "Block type effects"),
  extract_contrast_results(pairs(rl_se_block_effects, by = "stim"), "Serotonin", "Reversal learning", "Block type effects"),
  extract_contrast_results(pairs(rl_ne_block_effects, by = "stim"), "Norepinephrine", "Reversal learning", "Block type effects"),

  extract_contrast_results(pairs(rl_da_pe_effects, by = "stim"), "Dopamine", "Reversal learning", "Prediction error effects"),
  extract_contrast_results(pairs(rl_se_pe_effects, by = "stim"), "Serotonin", "Reversal learning", "Prediction error effects"),
  extract_contrast_results(pairs(rl_ne_pe_effects, by = "stim"), "Norepinephrine", "Reversal learning", "Prediction error effects"),

  extract_contrast_results(pairs(rl_da_rt_effects, by = "stim"), "Dopamine", "Reversal learning", "Response time effects"),
  extract_contrast_results(pairs(rl_se_rt_effects, by = "stim"), "Serotonin", "Reversal learning", "Response time  effects"),
  extract_contrast_results(pairs(rl_ne_rt_effects, by = "stim"), "Norepinephrine", "Reversal learning", "Response time  effects"),

  extract_contrast_results(pairs(rl_da_error_effects, by = "stim"), "Dopamine", "Reversal learning", "Choice effects"),
  extract_contrast_results(pairs(rl_se_error_effects, by = "stim"), "Serotonin", "Reversal learning", "Choice effects"),
  extract_contrast_results(pairs(rl_ne_error_effects, by = "stim"), "Norepinephrine", "Reversal learning", "Choice effects"),

  extract_contrast_results(pairs(rl_da_rew_effects, by = "stim"), "Dopamine", "Reversal learning", "Reward effects"),
  extract_contrast_results(pairs(rl_se_rew_effects, by = "stim"), "Serotonin", "Reversal learning", "Reward effects"),
  extract_contrast_results(pairs(rl_ne_rew_effects, by = "stim"), "Norepinephrine", "Reversal learning", "Reward effects"),

  extract_contrast_results(pairs(rl_da_reversal_effects, by = "stim"), "Dopamine", "Reversal learning", "Reversal effects"),
  extract_contrast_results(pairs(rl_se_reversal_effects, by = "stim"), "Serotonin", "Reversal learning", "Reversal effects"),
  extract_contrast_results(pairs(rl_ne_reversal_effects, by = "stim"), "Norepinephrine", "Reversal learning", "Reversal effects"),

  extract_contrast_results(pairs(rl_da_outcome_effects, by = "stim"), "Dopamine", "Reversal learning", "Outcome effects"),
  extract_contrast_results(pairs(rl_se_outcome_effects, by = "stim"), "Serotonin", "Reversal learning", "Outcome effects"),
  extract_contrast_results(pairs(rl_ne_outcome_effects, by = "stim"), "Norepinephrine", "Reversal learning", "Outcome effects"),

  extract_contrast_results(pairs(ug_da_outcome_effects, by = "stim"), "Dopamine", "Ultimatum game", "Outcome effects"),
  extract_contrast_results(pairs(ug_se_outcome_effects, by = "stim"), "Serotonin", "Ultimatum game", "Outcome effects"),
  extract_contrast_results(pairs(ug_ne_outcome_effects, by = "stim"), "Norepinephrine", "Ultimatum game", "Outcome effects"),

  extract_contrast_results(pairs(ug_da_rt_effects, by = "stim"), "Dopamine", "Ultimatum game", "Resposne time effects"),
  extract_contrast_results(pairs(ug_se_rt_effects, by = "stim"), "Serotonin", "Ultimatum game", "Resposne time effects"),
  extract_contrast_results(pairs(ug_ne_rt_effects, by = "stim"), "Norepinephrine", "Ultimatum game", "Resposne time effects"),

  extract_contrast_results(pairs(ug_da_accept_effects, by = "stim"), "Dopamine", "Ultimatum game", "Choice effects"),
  extract_contrast_results(pairs(ug_se_accept_effects, by = "stim"), "Serotonin", "Ultimatum game", "Choice effects"),
  extract_contrast_results(pairs(ug_ne_accept_effects, by = "stim"), "Norepinephrine", "Ultimatum game", "Choice effects"),

  extract_contrast_results(pairs(ug_da_offer_effects, by = "stim"), "Dopamine", "Ultimatum game", "Offer effects"),
  extract_contrast_results(pairs(ug_se_offer_effects, by = "stim"), "Serotonin", "Ultimatum game", "Offer effects"),
  extract_contrast_results(pairs(ug_ne_offer_effects, by = "stim"), "Norepinephrine", "Ultimatum game", "Offer effects"),

  extract_contrast_results(pairs(ug_da_oc_effects, by = "stim"), "Dopamine", "Ultimatum game", "Offer change effects"),
  extract_contrast_results(pairs(ug_se_oc_effects, by = "stim"), "Serotonin", "Ultimatum game", "Offer change effects"),
  extract_contrast_results(pairs(ug_ne_oc_effects, by = "stim"), "Norepinephrine", "Ultimatum game", "Offer change effects")

)

combinations <- create_all_contrast_tables(all_contrasts)
