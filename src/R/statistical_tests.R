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
# 2025/05/05    Blair Shevlin                           wrote original code
# 2025/07/24    Blair Shevlin                           updated to use new NT data
# 2025/09/16    Blair Shevlin                           updated to use revised NT data
# 2025/09/17    Blair Shevlin                           robustness tests for clinical predictions
# 2025/10/1     Blair Shevlin                           statistical tests comparing post-stim OR to baseline
# 2025/11/14    Blair Shevlin                           adding effect sizes
# 2025/11/18    Blair Shevlin                           added task x sess analysis
# 2026/05/05    Blair Shevlin                           added HDRS change analysis

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

# Paths
dir = path(here())
nt_dir = dir / "data" / "nt" / "processed"
beh_dir = dir / "data" / "behavior"
clin_dir = dir / "data" / "clinical"

# Subject-level NT Data
load(nt_dir / "UG_RL_NT-Continuous_9-23-25.RData")

# Generate summary data - focusing only on reward events (RL) and offer events (UG)
rl.EST.Reward = rl.EST %>%
  filter(event == "Reward") %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))

ug.EST.Offer = ug.EST %>%
  filter(event == "Offer") %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim"))) 

# IDs
ids_final = unique(ug.EST.Offer$idx)

# Load clinical data
cl.df = read.csv(clin_dir / "clinical-data_deid_07-10-25.csv" ) 

# Load behavioral data
rl.beh = read.csv(file = beh_dir / "rl-data_deid_07-10-25.csv")

ug.beh = read.csv(file = beh_dir / "ug-data_deid_07-10-25.csv")

##############
# STIM ON NT #
##############

# Regression
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
         trial_z = scale(trial)[,1]
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

# Summarise NT change
rl.EST.pr %>% 
  mutate(nt = factor(nt, levels = c("DA","SE","NE"))) %>%
  group_by(idx, nt) %>%
  summarise(
    nt_change = mean(Oz[stim == "Post-Stim"]) - mean(Oz[stim == "Pre-Stim"]),
    n = n(),
    se = sd(Oz[stim == "Post-Stim"] - Oz[stim == "Pre-Stim"]) / sqrt(n)
  ) %>% 
  as.data.frame()

ug.EST.pr %>% 
  mutate(nt = factor(nt, levels = c("DA","SE","NE"))) %>%
  group_by(idx, nt) %>%
  summarise(
    nt_change = mean(Oz[stim == "Post-Stim"]) - mean(Oz[stim == "Pre-Stim"]),
    n = n(),
    se = sd(Oz[stim == "Post-Stim"] - Oz[stim == "Pre-Stim"]) / sqrt(n)
  ) %>% 
  as.data.frame()

# Set up contrast coding for stim
contrasts(rl.EST.pr$stim) <- c(-1,1)
contrasts(ug.EST.pr$stim) <- c(-1,1)

# ANOVA
rl.EST.aov = rl.EST.pr %>% group_by(idx,stim,nt,block_type) %>%
  summarise(oz = mean(Oz))
ug.EST.aov = ug.EST.pr %>% group_by(idx,stim,nt,offer_bin_f) %>%
  summarise(oz = mean(Oz))

rl.DA.anova=aov(oz ~ stim * block_type,
                data = rl.EST.aov[rl.EST.aov$nt == "DA",])
rl.SE.anova=aov(oz ~ stim * block_type,
                data = rl.EST.aov[rl.EST.aov$nt == "SE",])

ug.DA.anova=aov(oz ~ stim * offer_bin_f,
                data = ug.EST.aov[ug.EST.aov$nt == "DA",])
ug.SE.anova=aov(oz ~ stim * offer_bin_f,
                data = ug.EST.aov[ug.EST.aov$nt == "SE",])

summary(rl.DA.anova)
summary(rl.SE.anova)

summary(ug.DA.anova)
summary(ug.SE.anova)

# Effect sizes
eta_squared(rl.DA.anova, partial = TRUE)
eta_squared(rl.SE.anova, partial = TRUE)
eta_squared(ug.DA.anova, partial = TRUE)
eta_squared(ug.SE.anova, partial = TRUE)

# Test task-by-stim interaction
rl.premerge = rl.EST.pr %>%
  select(idx,stim,nt,Oz,trial_z) %>%
  mutate(task = "RL")
ug.premerge = ug.EST.pr %>%
  select(idx,stim,nt,Oz,trial_z) %>%
  mutate(task = "UG")

task.merge = rbind(rl.premerge,ug.premerge) %>% 
            mutate(task = factor(task))

task.merge.aov = task.merge %>% group_by(idx,stim,nt,task) %>%
  summarise(oz = mean(Oz))         

DA.anova=aov(oz ~ stim * task,
                data = task.merge.aov[task.merge.aov$nt == "DA",])
SE.anova=aov(oz ~ stim * task,
                data = task.merge.aov[task.merge.aov$nt == "SE",])  

summary(DA.anova)
summary(SE.anova)    

TukeyHSD(DA.anova)
TukeyHSD(SE.anova)


rl.lm.DA = lmer(Oz ~ stim * trial_z +
                    (1+stim|idx),
               data = task.merge[task.merge$nt == "DA" & task.merge$task == "RL",],
               REML=F,
               control = lmerControl(optimizer = "bobyqa"))
rl.lm.SE = lmer(Oz ~ stim * trial_z +
                    (1+stim|idx),
                  data = task.merge[task.merge$nt == "SE" & task.merge$task == "RL",],
                  REML=F,
                  control = lmerControl(optimizer = "bobyqa"))
rl.lm.NE = lmer(Oz ~ stim * trial_z +
                    (1+stim|idx),
                  data = task.merge[task.merge$nt == "NE" & task.merge$task == "RL",],
                  REML=F,
                  control = lmerControl(optimizer = "bobyqa"))

ug.lm.DA = lmer(Oz ~ stim * trial_z +
                    (1+stim|idx),
               data = task.merge[task.merge$nt == "DA" & task.merge$task == "UG",],
               REML=F,
               control = lmerControl(optimizer = "bobyqa"))
ug.lm.SE = lmer(Oz ~ stim * trial_z +
                    (1+stim|idx),
                  data = task.merge[task.merge$nt == "SE" & task.merge$task == "UG",],
                  REML=F,
                  control = lmerControl(optimizer = "bobyqa"))
ug.lm.NE = lmer(Oz ~ stim * trial_z +
                  (1+stim|idx),
                  data = task.merge[task.merge$nt == "NE" & task.merge$task == "UG",],
                  REML=F,
                  control = lmerControl(optimizer = "bobyqa"))  

summary(rl.lm.DA)                  
summary(rl.lm.SE)                  
summary(rl.lm.NE) 

summary(ug.lm.DA)                  
summary(ug.lm.SE)  
summary(ug.lm.NE)  

sess.lm.DA = lmer(Oz ~ stim * task *  trial_z +
                    (1+stim+task|idx),
               data = task.merge[task.merge$nt == "DA",],
               REML=F,
               control = lmerControl(optimizer = "bobyqa"))
sess.lm.SE = lmer(Oz ~ stim * task * trial_z +
                    (1+stim+task|idx),
                  data = task.merge[task.merge$nt == "SE",],
                  REML=F,
                  control = lmerControl(optimizer = "bobyqa"))
sess.lm.NE = lmer(Oz ~ stim * task * trial_z +
                    (1+stim+task|idx),
                  data = task.merge[task.merge$nt == "NE",],
                  REML=F,
                  control = lmerControl(optimizer = "bobyqa"))

summary(sess.lm.DA)
summary(sess.lm.SE)
summary(sess.lm.NE)

sess.lm.SE %>%
    tidy(effects = "fixed") %>%
    mutate(p_adjusted = p.adjust(p.value) )

# Simple slopes analysis for stim effects
da_stim <- emmeans(sess.lm.DA, ~ stim * task)
se_stim <- emmeans(sess.lm.SE, ~ stim * task)
ne_stim <- emmeans(sess.lm.NE, ~ stim * task)

pairs(da_stim, by = "task")
pairs(se_stim, by = "task")
pairs(ne_stim, by = "task")

###
# Test Order effects #
task.merge = task.merge %>%
  mutate(order = ifelse(stim == "Pre-Stim", ifelse(task == "RL", "Task 1", "Task 2"),
                        ifelse(task == "RL", "Task 2", "Task 1")),
         order = factor(order))

order.lm.DA = lmer(Oz ~ stim * order  +
                    (1+stim+order|idx),
               data = task.merge[task.merge$nt == "DA",],
               REML=F,
               control = lmerControl(optimizer = "bobyqa"))
order.lm.SE = lmer(Oz ~ stim * order * trial_z +
                    (1+stim+order|idx),
                  data = task.merge[task.merge$nt == "SE",],
                  REML=F,
                  control = lmerControl(optimizer = "bobyqa"))
order.lm.NE = lmer(Oz ~ stim * order * trial_z +
                    (1+stim+order|idx),
                  data = task.merge[task.merge$nt == "NE",],
                  REML=F,
                  control = lmerControl(optimizer = "bobyqa"))

summary(order.lm.DA) # No main effect of order or interactions with stim
summary(order.lm.SE) # No main effect of order or interactions with stim (but interaction with trial)
summary(order.lm.NE) # No main effect of order or interaction with stim



##################
# BEH PROCESSING #
##################

# Get session-level averages
ug.beh.means = ug.beh %>% group_by(idx,sess) %>%
  filter(rt < 10) %>%
  summarise(mChoice = mean(rej==0),
            mRT = mean(rt),
            mLogRT = mean(log(rt))) %>%   
  mutate(idx = factor(idx),
         sess = recode(sess,"fMRI" = "Baseline")) %>%
  as.data.frame()

ug.mood.means = ug.beh %>% group_by(idx,sess) %>%
  filter(rt_mood < 10) %>%
  filter(!is.na(rt_mood)) %>%
  summarise(mMood = mean(mood),
            mRT = mean(rt_mood ),
            mLogRT = mean(log(rt))) %>% 
  mutate(idx = factor(idx),
         sess = recode(sess,"fMRI" = "Baseline")) %>%
  as.data.frame()

rl.beh.means = rl.beh %>%
  filter(rt < 10) %>%
  group_by(idx,sess) %>%
  summarise(mChoice = mean(opt),
            mRT = mean(rt),
            mLogRT = mean(log(rt))) %>% 
  mutate(idx = factor(idx),
         sess = recode(sess,"fMRI" = "Baseline")) %>%
  as.data.frame()


#########################
# LONGITUDINAL BEHAVIOR #
#########################

# Function to calculate t.test

calculate_t <- function(data, value_col, time_col, baseline_level, 
                        comparison_level, subject_col = "idx", 
                        paired = TRUE,
                        alternative = "two.sided") {
  # Extract baseline and comparison data
  baseline_data <- data[data[[time_col]] == baseline_level, ]
  comparison_data <- data[data[[time_col]] == comparison_level, ]
  
  # Ensure data is properly paired by subject
  merged_data <- merge(
    baseline_data[, c(subject_col, value_col)],
    comparison_data[, c(subject_col, value_col)],
    by = subject_col
  )
  
  # Perform paired t-test
  t_test <- t.test(merged_data[[paste0(value_col, ".y")]],
                    merged_data[[paste0(value_col, ".x")]],
                   paired = paired,
                   alternative = alternative)
  
  return(t_test)
}

# RL - RT
calculate_t(rl.beh.means,value_col = "mLogRT", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Post-Stim")

calculate_t(rl.beh.means,value_col = "mLogRT", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 1")

calculate_t(rl.beh.means,value_col = "mLogRT", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 3")

calculate_t(rl.beh.means,value_col = "mLogRT", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 6")
            
rl.beh.means %>% group_by(sess) %>%
  summarise(meanRT = mean(mRT),
            sdRT = sd(mRT),
            n = n(),
            seRT = sdRT/sqrt(n),
            ciLower = meanRT - qt(0.975, df=n-1)*seRT,
            ciUpper = meanRT + qt(0.975, df=n-1)*seRT)
# UG - RT
calculate_t(ug.beh.means,value_col = "mLogRT", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Post-Stim")

calculate_t(ug.beh.means,value_col = "mLogRT", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 1")

calculate_t(ug.beh.means,value_col = "mLogRT", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 3")

calculate_t(ug.beh.means,value_col = "mLogRT", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 6")

ug.beh.means %>% group_by(sess) %>%
  summarise(meanRT = mean(mRT),
            sdRT = sd(mRT),
            n = n(),
            seRT = sdRT/sqrt(n),
            ciLower = meanRT - qt(0.975, df=n-1)*seRT,
            ciUpper = meanRT + qt(0.975, df=n-1)*seRT)

# RL - Choice
calculate_t(rl.beh.means,value_col = "mChoice", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Post-Stim")

calculate_t(rl.beh.means,value_col = "mChoice", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 1")

calculate_t(rl.beh.means,value_col = "mChoice", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 3")

calculate_t(rl.beh.means,value_col = "mChoice", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 6")

rl.beh.means %>% group_by(sess) %>%
  summarise(mOpt = mean(mChoice),
            sdOpt = sd(mChoice),
            n = n(),
            seAccept = sdOpt/sqrt(n),
            ciLower = mOpt - qt(0.975, df=n-1)*sdOpt,
            ciUpper = mOpt + qt(0.975, df=n-1)*sdOpt)

# UG - Choice
calculate_t(ug.beh.means,value_col = "mChoice", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Post-Stim")

calculate_t(ug.beh.means,value_col = "mChoice", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 1")

calculate_t(ug.beh.means,value_col = "mChoice", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 3")

calculate_t(ug.beh.means,value_col = "mChoice", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 6")

ug.beh.means %>% group_by(sess) %>%
  summarise(mAccept = mean(mChoice),
            sdAccept = sd(mChoice),
            n = n(),
            seAccept = sdAccept/sqrt(n),
            ciLower = mAccept - qt(0.975, df=n-1)*seAccept,
            ciUpper = mAccept + qt(0.975, df=n-1)*seAccept)

# UG - Mood
calculate_t(ug.mood.means,value_col = "mMood", time_col = "sess",
            baseline_level = "Pre-Stim",
            comparison_level = "Post-Stim")
calculate_t(ug.mood.means,value_col = "mMood", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Post-Stim")
calculate_t(ug.mood.means,value_col = "mMood", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 6")

round(sd(ug.mood.means$mMood[ug.mood.means$sess == "Baseline"]),2)
round(sd(ug.mood.means$mMood[ug.mood.means$sess == "Pre-Stim"]),2)
round(sd(ug.mood.means$mMood[ug.mood.means$sess == "Post-Stim"]),2)
round(sd(ug.mood.means$mMood[ug.mood.means$sess == "Week 1"]),2)
round(sd(ug.mood.means$mMood[ug.mood.means$sess == "Month 1"]),2)
round(sd(ug.mood.means$mMood[ug.mood.means$sess == "Month 2"]),2)
round(sd(ug.mood.means$mMood[ug.mood.means$sess == "Month 3"]),2)
round(sd(ug.mood.means$mMood[ug.mood.means$sess == "Month 4"]),2)
round(sd(ug.mood.means$mMood[ug.mood.means$sess == "Month 5"]),2)
round(sd(ug.mood.means$mMood[ug.mood.means$sess == "Month 6"]),2)

round(mean(ug.mood.means$mMood[ug.mood.means$sess == "Baseline"]),2)
round(mean(ug.mood.means$mMood[ug.mood.means$sess == "Pre-Stim"]),2)
round(mean(ug.mood.means$mMood[ug.mood.means$sess == "Post-Stim"]),2)
round(mean(ug.mood.means$mMood[ug.mood.means$sess == "Week 1"]),2)
round(mean(ug.mood.means$mMood[ug.mood.means$sess == "Month 1"]),2)
round(mean(ug.mood.means$mMood[ug.mood.means$sess == "Month 2"]),2)
round(mean(ug.mood.means$mMood[ug.mood.means$sess == "Month 3"]),2)
round(mean(ug.mood.means$mMood[ug.mood.means$sess == "Month 4"]),2)
round(mean(ug.mood.means$mMood[ug.mood.means$sess == "Month 5"]),2)
round(mean(ug.mood.means$mMood[ug.mood.means$sess == "Month 6"]),2)


ug.mood.means %>% group_by(sess) %>%
  dplyr::summarise(mMood = round(mean(mMood),3),
           # sdMood = sd(mMood, na.rm = T),
            n = n(),
            #seMood = sdMood/sqrt(n),
            #ciLower = mMood - qt(0.975, df=n-1)*seMood,
            #ciUpper = mMood + qt(0.975, df=n-1)*seMood
            )

#################
# NT ~ BEHAVIOR #
################

# Simplify changes in NT across sessions and task (by NT type)
rl.nt = rl.EST.Reward %>%
  filter(nt %in% c("DA","SE")) %>%
  group_by(idx,nt,stim) %>%
  summarise(est = mean(Oz)) %>%
  group_by(idx,nt) %>%
  summarise(dff = est[stim == "Post-Stim"] - est[stim == "Pre-Stim"]) %>%
  group_by(idx) %>%
  reframe(DA = dff[nt == "DA"],
          SE = dff[nt == "SE"])

ug.nt = ug.EST.Offer %>%
  filter(nt %in% c("DA","SE")) %>%
  group_by(idx,nt,stim) %>%
  summarise(est = mean(Oz)) %>%
  group_by(idx,nt) %>%
  summarise(dff = est[stim == "Post-Stim"] - est[stim == "Pre-Stim"]) %>%
  group_by(idx) %>%
  reframe(DA = dff[nt == "DA"],
          SE = dff[nt == "SE"])

rl.beh.change = 
  rl.beh.means %>%
  mutate(mRT = mLogRT) %>%
  select(!mLogRT) %>%
  pivot_wider(names_from = sess, values_from = c(mChoice, mRT)) %>%
  mutate(# Use Baseline as baseline
    mChoice_PostStim_Change = `mChoice_Post-Stim` - mChoice_Baseline,
    mChoice_W1_Change = `mChoice_Week 1` - mChoice_Baseline,
    mChoice_M1_Change = `mChoice_Month 1` - mChoice_Baseline,
    mChoice_M2_Change = `mChoice_Month 2` - mChoice_Baseline,
    mChoice_M3_Change = `mChoice_Month 3` - mChoice_Baseline,
    mChoice_M4_Change = `mChoice_Month 4` - mChoice_Baseline,
    mChoice_M5_Change = `mChoice_Month 5` - mChoice_Baseline,
    mChoice_M6_Change = `mChoice_Month 6` - mChoice_Baseline,
    mRT_PostStim_Change = `mRT_Post-Stim` - mRT_Baseline,
    mRT_W1_Change = `mRT_Week 1` - mRT_Baseline,
    mRT_M1_Change = `mRT_Month 1` - mRT_Baseline,
    mRT_M2_Change = `mRT_Month 2` - mRT_Baseline,
    mRT_M3_Change = `mRT_Month 3` - mRT_Baseline,
    mRT_M4_Change = `mRT_Month 4` - mRT_Baseline,
    mRT_M5_Change = `mRT_Month 5` - mRT_Baseline,
    mRT_M6_Change = `mRT_Month 6` - mRT_Baseline,
    mChoice_OR_Change = `mChoice_Post-Stim` - `mChoice_Pre-Stim`,
    mRT_OR_Change = `mRT_Post-Stim` - `mRT_Pre-Stim`
  ) %>%
  select(c(idx,mChoice_PostStim_Change,mChoice_W1_Change,mChoice_M1_Change,mChoice_M2_Change,
           mChoice_M3_Change,mChoice_M4_Change,mChoice_M5_Change,mChoice_M6_Change,
           mRT_PostStim_Change,mRT_W1_Change,mRT_M1_Change,mRT_M2_Change,
           mRT_M3_Change,mRT_M4_Change,mRT_M5_Change,mRT_M6_Change,mChoice_OR_Change,mRT_OR_Change
  ))

rl.beh.nt = merge(rl.beh.change, rl.nt) 

rl.beh.nt.long = rl.beh.nt %>% 
  pivot_longer(cols = c(DA,SE), names_to = "NT", values_to = "ests") %>%
  pivot_longer(cols = c(mChoice_PostStim_Change,mChoice_W1_Change,mChoice_M1_Change,mChoice_M2_Change,
                        mChoice_M3_Change,mChoice_M4_Change,mChoice_M5_Change,mChoice_M6_Change,
                        mRT_PostStim_Change,mRT_W1_Change,mRT_M1_Change,mRT_M2_Change,
                        mRT_M3_Change,mRT_M4_Change,mRT_M5_Change,mRT_M6_Change,mChoice_OR_Change,mRT_OR_Change), names_to = "change") %>%
  mutate(month = str_extract(change, "(?<=_)[^_]+(?=_)"),
         beh = ifelse(grepl("mChoice", change), "choice", "rt"))

ug.beh.change = 
  ug.beh.means %>%
  mutate(mRT = mLogRT) %>%
  select(!mLogRT) %>%
  pivot_wider(names_from = sess, values_from = c(mChoice, mRT)) %>%
  mutate(# Use fMRI as baseline
    mChoice_PostStim_Change = `mChoice_Post-Stim` - mChoice_Baseline,
    mChoice_W1_Change = `mChoice_Week 1` - mChoice_Baseline,
    mChoice_M1_Change = `mChoice_Month 1` - mChoice_Baseline,
    mChoice_M2_Change = `mChoice_Month 2` - mChoice_Baseline,
    mChoice_M3_Change = `mChoice_Month 3` - mChoice_Baseline,
    mChoice_M4_Change = `mChoice_Month 4` - mChoice_Baseline,
    mChoice_M5_Change = `mChoice_Month 5` - mChoice_Baseline,
    mChoice_M6_Change = `mChoice_Month 6` - mChoice_Baseline,
    mRT_PostStim_Change = `mRT_Post-Stim` - mRT_Baseline,
    mRT_W1_Change = `mRT_Week 1` - mRT_Baseline,
    mRT_M1_Change = `mRT_Month 1` - mRT_Baseline,
    mRT_M2_Change = `mRT_Month 2` - mRT_Baseline,
    mRT_M3_Change = `mRT_Month 3` - mRT_Baseline,
    mRT_M4_Change = `mRT_Month 4` - mRT_Baseline,
    mRT_M5_Change = `mRT_Month 5` - mRT_Baseline,
    mRT_M6_Change = `mRT_Month 6` - mRT_Baseline,
    mChoice_OR_Change = `mChoice_Post-Stim` - `mChoice_Pre-Stim`,
    mRT_OR_Change = `mRT_Post-Stim` - `mRT_Pre-Stim`
  ) %>%
  select(c(idx,mChoice_PostStim_Change,mChoice_W1_Change,mChoice_M1_Change,mChoice_M2_Change,
           mChoice_M3_Change,mChoice_M4_Change,mChoice_M5_Change,mChoice_M6_Change,
           mRT_PostStim_Change,mRT_W1_Change,mRT_M1_Change,mRT_M2_Change,
           mRT_M3_Change,mRT_M4_Change,mRT_M5_Change,mRT_M6_Change,mChoice_OR_Change,mRT_OR_Change ))

ug.beh.nt = merge(ug.beh.change, ug.nt) 

ug.beh.nt.long = ug.beh.nt %>% 
  pivot_longer(cols = c(DA,SE), names_to = "NT", values_to = "ests") %>%
  pivot_longer(cols = c(mChoice_PostStim_Change,mChoice_W1_Change,mChoice_M1_Change,mChoice_M2_Change,
                        mChoice_M3_Change,mChoice_M4_Change,mChoice_M5_Change,mChoice_M6_Change,
                        mRT_PostStim_Change,mRT_W1_Change,mRT_M1_Change,mRT_M2_Change,
                        mRT_M3_Change,mRT_M4_Change,mRT_M5_Change,mRT_M6_Change,mChoice_OR_Change,mRT_OR_Change), 
               names_to = "change") %>%
  mutate(month = str_extract(change, "(?<=_)[^_]+(?=_)"),
         beh = ifelse(grepl("mChoice", change), "choice", "rt"),
         baseline = ifelse(grepl("Change2",change),"w1","fmri"))


# UG Mood
ug.mood.change = 
  ug.mood.means %>%
  mutate(mRT = mLogRT) %>%
  select(!mLogRT) %>%
  pivot_wider(names_from = sess, values_from = c(mMood, mRT)) %>%
  mutate(# Use fMRI as baseline
    mMood_PostStim_Change = `mMood_Post-Stim` - mMood_Baseline,
    mMood_W1_Change = `mMood_Week 1` - mMood_Baseline,
    mMood_M1_Change = `mMood_Month 1` - mMood_Baseline,
    mMood_M2_Change = `mMood_Month 2` - mMood_Baseline,
    mMood_M3_Change = `mMood_Month 3` - mMood_Baseline,
    mMood_M4_Change = `mMood_Month 4` - mMood_Baseline,
    mMood_M5_Change = `mMood_Month 5` - mMood_Baseline,
    mMood_M6_Change = `mMood_Month 6` - mMood_Baseline,
    mRT_PostStim_Change = `mRT_Post-Stim` - mRT_Baseline,
    mRT_W1_Change = `mRT_Week 1` - mRT_Baseline,
    mRT_M1_Change = `mRT_Month 1` - mRT_Baseline,
    mRT_M2_Change = `mRT_Month 2` - mRT_Baseline,
    mRT_M3_Change = `mRT_Month 3` - mRT_Baseline,
    mRT_M4_Change = `mRT_Month 4` - mRT_Baseline,
    mRT_M5_Change = `mRT_Month 5` - mRT_Baseline,
    mRT_M6_Change = `mRT_Month 6` - mRT_Baseline,
    mMood_OR_Change = `mMood_Post-Stim` - `mMood_Pre-Stim`,
    mRT_OR_Change = `mRT_Post-Stim` - `mRT_Pre-Stim`
  ) %>%
  select(c(idx,mMood_PostStim_Change,mMood_W1_Change,mMood_M1_Change,mMood_M2_Change,
           mMood_M3_Change,mMood_M4_Change,mMood_M5_Change,mMood_M6_Change,
           mRT_PostStim_Change,mRT_W1_Change,mRT_M1_Change,mRT_M2_Change,
           mRT_M3_Change,mRT_M4_Change,mRT_M5_Change,mRT_M6_Change,mMood_OR_Change,mRT_OR_Change ))

ug.mood.nt = merge(ug.mood.change, ug.nt) 

ug.mood.nt.long = ug.mood.nt %>% 
  pivot_longer(cols = c(DA,SE), names_to = "NT", values_to = "ests") %>%
  pivot_longer(cols = c(mMood_PostStim_Change,mMood_W1_Change,mMood_M1_Change,mMood_M2_Change,
                        mMood_M3_Change,mMood_M4_Change,mMood_M5_Change,mMood_M6_Change,
                        mRT_PostStim_Change,mRT_W1_Change,mRT_M1_Change,mRT_M2_Change,
                        mRT_M3_Change,mRT_M4_Change,mRT_M5_Change,mRT_M6_Change,mMood_OR_Change,mRT_OR_Change), 
               names_to = "change") %>%
  mutate(month = str_extract(change, "(?<=_)[^_]+(?=_)"),
         beh = ifelse(grepl("mMood", change), "mood", "rt"),
         baseline = ifelse(grepl("Change2",change),"w1","fmri"))


# Run correlations and tidy the results

# Longitudinal changes
rl.beh.nt.long %>%
  group_by(NT,beh,month) %>%
  do(tidy(cor.test(.$value, .$ests, method = "spearman"))) %>%
  # Filter for significant results (p < 0.05)
  filter(p.value < .05)

# All SE~RT for reporting
rl.beh.nt.long %>%
  group_by(NT,beh,month) %>%
  do(tidy(cor.test(.$value, .$ests, method = "spearman"))) %>%
  # Filter for significant results (p < 0.05)
  filter(NT == "SE", beh == "rt")

ug.beh.nt.long %>%
  group_by(NT,beh,month) %>%
  do(tidy(cor.test(.$value, .$ests, method = "spearman"))) %>%
  # Filter for significant results (p < 0.05)
  filter(p.value < .05)

ug.mood.nt.long %>%
  group_by(NT,beh,month) %>%
  do(tidy(cor.test(.$value, .$ests, method = "spearman"))) %>%
  # Filter for significant results (p < 0.05)
  filter(p.value < .05)

# All DA ~ mood
ug.mood.nt.long %>%
  group_by(NT,beh,month) %>%
  do(tidy(cor.test(.$value, .$ests, method = "spearman"))) %>%
  # Filter for significant results (p < 0.05)
  filter(NT == "DA", beh == "mood")

#################
# NT ~ CLINICAL #
#################

ug.nt.cl = 
  ug.EST.Offer %>% 
  pivot_longer(cols = c("Oz","Rz","Pz","Mz"),
               names_to = "nt_metric", values_to = "nt_val") %>%
  group_by(idx,stim,nt,nt_metric) %>%
  summarise(mTrial = mean(nt_val)) %>%
  filter(nt_metric == "Oz") %>%
  select(idx,stim,nt,mTrial) %>%
  pivot_wider(values_from = mTrial, names_from = c("stim","nt")) %>%
  mutate(DA_Pre_UG = `Pre-Stim_DA`,
         DA_Post_UG = `Post-Stim_DA`,
         deltaDA_UG = DA_Post_UG - DA_Pre_UG,
         SE_Pre_UG = `Pre-Stim_SE`,
         SE_Post_UG = `Post-Stim_SE`,
         deltaSE_UG = SE_Post_UG - SE_Pre_UG,
         NE_Pre_UG = `Pre-Stim_NE`,
         NE_Post_UG = `Post-Stim_NE`,
         deltaNE_UG = NE_Post_UG - NE_Pre_UG,
  ) %>%
  select(idx,
         DA_Pre_UG,DA_Post_UG,deltaDA_UG,
         SE_Pre_UG,SE_Post_UG,deltaSE_UG,
         NE_Pre_UG,NE_Post_UG,deltaNE_UG
  ) 

rl.nt.cl = 
  rl.EST.Reward %>%
  pivot_longer(cols = c("Oz","Rz","Pz","Mz"),
               names_to = "nt_metric", values_to = "nt_val") %>%
  group_by(idx,stim,nt,nt_metric) %>%
  summarise(mTrial = mean(nt_val)) %>%
  filter(nt_metric == "Oz") %>%
  select(idx,stim,nt,mTrial) %>%
  pivot_wider(values_from = mTrial, names_from = c("stim","nt")) %>%
  mutate(DA_Pre_RL = `Pre-Stim_DA`,
         DA_Post_RL = `Post-Stim_DA`,
         deltaDA_RL = DA_Post_RL - DA_Pre_RL,
         SE_Pre_RL = `Pre-Stim_SE`,
         SE_Post_RL = `Post-Stim_SE`,
         deltaSE_RL = SE_Post_RL - SE_Pre_RL,
         NE_Pre_RL = `Pre-Stim_NE`,
         NE_Post_RL = `Post-Stim_NE`,
         deltaNE_RL = NE_Post_RL - NE_Pre_RL,
  ) %>%
  select(idx,
         DA_Pre_RL,DA_Post_RL,deltaDA_RL,
         SE_Pre_RL,SE_Post_RL,deltaSE_RL,
         NE_Pre_RL,NE_Post_RL,deltaNE_RL
  ) 

nt.cl = merge(ug.nt.cl,rl.nt.cl,all.x = TRUE)

cl.Oz = merge(cl.df,nt.cl)

cl.Oz.lme = cl.Oz %>%
  filter(!session %in% c("post stim")) %>%
  mutate(sess_f = factor(session, levels  = c("fmri","pre stim", "week 1", "month 1", "month 2", "month 3", "month 4", "month 5", "month 6")) ) %>%
  group_by(idx) %>%
  mutate(baseline_HDRS = HDRS[session == "fmri"],#HDRS[session == "pre stim"],
         baseline_MADRS = MADRS[session == "fmri"],#MADRS[session == "pre stim"],
         m6_HDRS = HDRS[session == "month 6"],
         m6_MADRS = MADRS[session == "month 6"],
         deltaPerDA_UG = (DA_Post_UG - DA_Pre_UG)/DA_Pre_UG,
         deltaPerSE_UG = (SE_Post_UG - SE_Pre_UG)/SE_Pre_UG,
         deltaPerNE_UG = (NE_Post_UG - NE_Pre_UG)/NE_Pre_UG,
         deltaPerDA_RL = (DA_Post_RL - DA_Pre_RL)/DA_Pre_RL,
         deltaPerSE_RL = (SE_Post_RL - SE_Pre_RL)/SE_Pre_RL,
         deltaPerNE_RL = (NE_Post_RL - NE_Pre_RL)/NE_Pre_RL,
         deltaPerHDRS = (HDRS[session == "fmri"] - HDRS[session == "month 6"])/HDRS[session == "fmri"]
  ) %>%
  filter(!session %in% c("pre stim","fmri")) %>%
  ungroup()

data_m1 = cl.Oz.lme %>%
  filter(session == "month 1") %>%
  select(idx,HDRS,baseline_HDRS,deltaPerHDRS,
         deltaDA_UG,deltaSE_UG,deltaNE_UG,
         deltaDA_RL,deltaSE_RL,deltaNE_RL) 

data_m3 = cl.Oz.lme %>%
  filter(session == "month 3") %>%
  select(idx,HDRS,baseline_HDRS,deltaPerHDRS,
         deltaDA_UG,deltaSE_UG,deltaNE_UG,
         deltaDA_RL,deltaSE_RL,deltaNE_RL) 

data_m6 = cl.Oz.lme %>%
  filter(session == "month 6") %>%
  select(idx,HDRS,baseline_HDRS,deltaPerHDRS,
         deltaDA_UG,deltaSE_UG,deltaNE_UG,
         deltaDA_RL,deltaSE_RL,deltaNE_RL)

# Did HDRS change from baseline to month 6?
t.test(data_m6$HDRS, data_m6$baseline_HDRS, paired = TRUE)
# Effect size of difference
data_m6_wide = data_m6 %>%
  pivot_longer(cols = c(HDRS, baseline_HDRS), names_to = "time", values_to = "score")
cohens_d(score ~ time, data = data_m6_wide)

# Calculate models
model_RL_m1 <- lm(HDRS ~ deltaDA_RL * deltaSE_RL, data = data_m1)
model_RL_m3 <- lm(HDRS ~ deltaDA_RL * deltaSE_RL, data = data_m3)
model_RL_m6 <- lm(HDRS ~ deltaDA_RL * deltaSE_RL, data = data_m6)
model_UG_m1 <- lm(HDRS ~ deltaDA_UG * deltaSE_UG, data = data_m1)
model_UG_m3 <- lm(HDRS ~ deltaDA_UG * deltaSE_UG, data = data_m3)
model_UG_m6 <- lm(HDRS ~ deltaDA_UG * deltaSE_UG, data = data_m6)

# For reporting
summary(model_RL_m1)
summary(model_RL_m3)
summary(model_RL_m6)
summary(model_UG_m1)
summary(model_UG_m3)
summary(model_UG_m6)

# Compare alternative models

# Single NT models
model_DA = lm(HDRS ~ deltaDA_UG, data = data_m6)
model_SE = lm(HDRS ~ deltaSE_UG, data = data_m6)
model_NE = lm(HDRS ~ deltaNE_UG, data = data_m6)

# Additive models
model_DASE = lm(HDRS ~ deltaDA_UG + deltaSE_UG, data = data_m6)
model_DANE = lm(HDRS ~ deltaDA_UG + deltaNE_UG, data = data_m6)
model_NESE = lm(HDRS ~ deltaNE_UG + deltaSE_UG, data = data_m6)
model_DASENE = lm(HDRS ~ deltaDA_UG + deltaSE_UG + deltaNE_UG, data = data_m6)

# Interaction models
model_DAxSE = lm(HDRS ~ deltaDA_UG * deltaSE_UG, data = data_m6)
model_NExSE = lm(HDRS ~ deltaNE_UG * deltaSE_UG, data = data_m6)
model_DAxNE = lm(HDRS ~ deltaDA_UG * deltaNE_UG, data = data_m6)
model_DAxSExNE = lm(HDRS ~ deltaDA_UG * deltaSE_UG * deltaNE_UG, data = data_m6)

models <- list( model_DA, model_SE, model_NE,
                model_DASE,model_NESE,model_DANE,model_DASENE,
                model_DAxSE, model_NExSE, model_DAxNE, model_DAxSExNE)
model_names <- c("DA", "SE", "NE",
"DA+SE", "NE+SE", "DA+NE", "DA+SE+NE",
 "DA×SE", "NE×SE", "DA×NE", "DA×SE×NE")
comparison_table <- data.frame(
  Model = model_names,
  AIC = sapply(models, AIC),
  BIC = sapply(models, BIC),
  R_squared = sapply(models, function(x) summary(x)$r.squared),
  Adj_R_squared = sapply(models, function(x) summary(x)$adj.r.squared)
)

print(comparison_table)


# Calculate baseline models

baseline_model_RL <- lm(baseline_HDRS ~ deltaDA_RL * deltaSE_RL, data = data_m1)
baseline_model_UG <- lm(baseline_HDRS ~ deltaDA_UG * deltaSE_UG, data = data_m1)

# For reporting
summary(baseline_model_RL)
summary(baseline_model_UG)

AIC(baseline_model_RL) # 59.32
AIC(baseline_model_UG) # 57.31

BIC(baseline_model_RL) # 59.72
BIC(baseline_model_UG) # 58.82

