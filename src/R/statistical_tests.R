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
load(nt_dir / "UG_RL_NT-Continuous_7-14-25.RData")

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
         rew_f = factor(rew, levels=c(0,-10,10)))

ug.EST.pr = ug.EST.Offer %>%
  mutate(offer_bin_f = factor(offer_bin, levels = c("Middle","Low","High")),
         offer_z = scale(offer)[,1],
         offer_change = offer - prev_offer_raw,
         offer_change_z = scale(offer_change)[,1],
         offer_change_f = factor(ifelse(offer_change > 0, "improve","worse"),
                                 levels=c("worse","improve"))
  )

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

##################
# BEH PROCESSING #
##################

# Get session-level averages
ug.beh.means = ug.beh %>% group_by(idx,sess) %>%
  summarise(mChoice = mean(rej==0),
            mRT = mean(rt),
            mLogRT = mean(log(rt))) %>%   
  mutate(idx = factor(idx),
         sess = recode(sess,"fMRI" = "Baseline")) %>%
  as.data.frame()

ug.mood.means = ug.beh %>% group_by(idx,sess) %>%
  filter(!is.na(rt_mood)) %>%
  summarise(mMood = mean(mood),
            mRT = mean(rt_mood ),
            mLogRT = mean(log(rt))) %>% 
  mutate(idx = factor(idx),
         sess = recode(sess,"fMRI" = "Baseline")) %>%
  as.data.frame()

rl.beh.means = rl.beh %>%
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
  t_test <- wilcox.test(merged_data[[paste0(value_col, ".x")]], 
                   merged_data[[paste0(value_col, ".y")]], 
                   paired = paired,
                   alternative = alternative)
  
  return(t_test)
}

calculate_t(rl.beh.means,value_col = "mLogRT", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 1")

calculate_t(rl.beh.means,value_col = "mLogRT", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 3")

calculate_t(rl.beh.means,value_col = "mLogRT", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 6")

calculate_t(ug.beh.means,value_col = "mLogRT", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 1")

calculate_t(ug.beh.means,value_col = "mLogRT", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 3")

calculate_t(ug.beh.means,value_col = "mLogRT", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 6")

calculate_t(rl.beh.means,value_col = "mChoice", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 1")

calculate_t(rl.beh.means,value_col = "mChoice", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 3")

calculate_t(rl.beh.means,value_col = "mChoice", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 6")

calculate_t(ug.beh.means,value_col = "mChoice", time_col = "sess",
            baseline_level = "Baseline",
            comparison_level = "Month 6")

calculate_t(ug.mood.means,value_col = "mMood", time_col = "sess",
            baseline_level = "Month 6",
            comparison_level = "Baseline")


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
    mRT_M6_Change = `mRT_Month 6` - mRT_Baseline
  ) %>%
  select(c(idx,mChoice_PostStim_Change,mChoice_W1_Change,mChoice_M1_Change,mChoice_M2_Change,
           mChoice_M3_Change,mChoice_M4_Change,mChoice_M5_Change,mChoice_M6_Change,
           mRT_PostStim_Change,mRT_W1_Change,mRT_M1_Change,mRT_M2_Change,
           mRT_M3_Change,mRT_M4_Change,mRT_M5_Change,mRT_M6_Change,
  ))

rl.beh.nt = merge(rl.beh.change, rl.nt) 

rl.beh.nt.long = rl.beh.nt %>% 
  pivot_longer(cols = c(DA,SE), names_to = "NT", values_to = "ests") %>%
  pivot_longer(cols = c(mChoice_PostStim_Change,mChoice_W1_Change,mChoice_M1_Change,mChoice_M2_Change,
                        mChoice_M3_Change,mChoice_M4_Change,mChoice_M5_Change,mChoice_M6_Change,
                        mRT_PostStim_Change,mRT_W1_Change,mRT_M1_Change,mRT_M2_Change,
                        mRT_M3_Change,mRT_M4_Change,mRT_M5_Change,mRT_M6_Change), names_to = "change") %>%
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
    mRT_M6_Change = `mRT_Month 6` - mRT_Baseline
  ) %>%
  select(c(idx,mChoice_PostStim_Change,mChoice_W1_Change,mChoice_M1_Change,mChoice_M2_Change,
           mChoice_M3_Change,mChoice_M4_Change,mChoice_M5_Change,mChoice_M6_Change,
           mRT_PostStim_Change,mRT_W1_Change,mRT_M1_Change,mRT_M2_Change,
           mRT_M3_Change,mRT_M4_Change,mRT_M5_Change,mRT_M6_Change ))

ug.beh.nt = merge(ug.beh.change, ug.nt) 

ug.beh.nt.long = ug.beh.nt %>% 
  pivot_longer(cols = c(DA,SE), names_to = "NT", values_to = "ests") %>%
  pivot_longer(cols = c(mChoice_PostStim_Change,mChoice_W1_Change,mChoice_M1_Change,mChoice_M2_Change,
                        mChoice_M3_Change,mChoice_M4_Change,mChoice_M5_Change,mChoice_M6_Change,
                        mRT_PostStim_Change,mRT_W1_Change,mRT_M1_Change,mRT_M2_Change,
                        mRT_M3_Change,mRT_M4_Change,mRT_M5_Change,mRT_M6_Change), 
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
    mRT_M6_Change = `mRT_Month 6` - mRT_Baseline
  ) %>%
  select(c(idx,mMood_PostStim_Change,mMood_W1_Change,mMood_M1_Change,mMood_M2_Change,
           mMood_M3_Change,mMood_M4_Change,mMood_M5_Change,mMood_M6_Change,
           mRT_PostStim_Change,mRT_W1_Change,mRT_M1_Change,mRT_M2_Change,
           mRT_M3_Change,mRT_M4_Change,mRT_M5_Change,mRT_M6_Change ))

ug.mood.nt = merge(ug.mood.change, ug.nt) 

ug.mood.nt.long = ug.mood.nt %>% 
  pivot_longer(cols = c(DA,SE), names_to = "NT", values_to = "ests") %>%
  pivot_longer(cols = c(mMood_PostStim_Change,mMood_W1_Change,mMood_M1_Change,mMood_M2_Change,
                        mMood_M3_Change,mMood_M4_Change,mMood_M5_Change,mMood_M6_Change,
                        mRT_PostStim_Change,mRT_W1_Change,mRT_M1_Change,mRT_M2_Change,
                        mRT_M3_Change,mRT_M4_Change,mRT_M5_Change,mRT_M6_Change), 
               names_to = "change") %>%
  mutate(month = str_extract(change, "(?<=_)[^_]+(?=_)"),
         beh = ifelse(grepl("mMood", change), "mood", "rt"),
         baseline = ifelse(grepl("Change2",change),"w1","fmri"))


# Run correlations and tidy the results

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
  pivot_longer(cols = c("Oz","Rz","Pz","Mz","Totz"),
               names_to = "nt_metric", values_to = "nt_val") %>%
  group_by(idx,stim,nt,nt_metric) %>%
  summarise(mTrial = mean(nt_val)) %>%
  filter(nt_metric == "Totz") %>%
  select(idx,stim,nt,mTrial) %>%
  pivot_wider(values_from = mTrial, names_from = c("stim","nt")) %>%
  mutate(DA_Pre_UG = `Pre-Stim_DA`,
         DA_Post_UG = `Post-Stim_DA`,
         deltaDA_UG = DA_Post_UG - DA_Pre_UG,
         SE_Pre_UG = `Pre-Stim_SE`,
         SE_Post_UG = `Post-Stim_SE`,
         deltaSE_UG = SE_Post_UG - SE_Pre_UG,
         NE_Pre_UG = `Pre-Stim_NE`,
         NE_Post_UG = `Post-Stim_DA`,
         deltaNE_UG = NE_Post_UG - NE_Pre_UG,
  ) %>%
  select(idx,
         DA_Pre_UG,DA_Post_UG,deltaDA_UG,
         SE_Pre_UG,SE_Post_UG,deltaSE_UG,
         NE_Pre_UG,NE_Post_UG,deltaNE_UG
  ) 

rl.nt.cl = 
  rl.EST.Reward %>%
  pivot_longer(cols = c("Oz","Rz","Pz","Mz","Totz"),
               names_to = "nt_metric", values_to = "nt_val") %>%
  group_by(idx,stim,nt,nt_metric) %>%
  summarise(mTrial = mean(nt_val)) %>%
  filter(nt_metric == "Totz") %>%
  select(idx,stim,nt,mTrial) %>%
  pivot_wider(values_from = mTrial, names_from = c("stim","nt")) %>%
  mutate(DA_Pre_RL = `Pre-Stim_DA`,
         DA_Post_RL = `Post-Stim_DA`,
         deltaDA_RL = DA_Post_RL - DA_Pre_RL,
         SE_Pre_RL = `Pre-Stim_SE`,
         SE_Post_RL = `Post-Stim_SE`,
         deltaSE_RL = SE_Post_RL - SE_Pre_RL,
         NE_Pre_RL = `Pre-Stim_NE`,
         NE_Post_RL = `Post-Stim_DA`,
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

data = cl.Oz.lme %>%
  filter(session == "month 6") %>%
  select(idx,HDRS,baseline_HDRS,deltaPerHDRS,
         deltaDA_UG,deltaSE_UG,deltaNE_UG,
         deltaDA_RL,deltaSE_RL,deltaNE_RL) 

# Create synergy score and additional metrics
data_enhanced <- data %>%
  mutate(
    # Synergy when both change in beneficial direction (DA+, SE+)
    synergy_score = case_when(
      deltaDA_UG > 0 & deltaSE_UG > 0 ~ deltaDA_UG * deltaSE_UG,
      deltaDA_UG < 0 & deltaSE_UG < 0 ~ abs(deltaDA_UG * deltaSE_UG),
      TRUE ~ -abs(deltaDA_UG * deltaSE_UG)  # Opposing changes
    ),
    
    # Alternative synergy: normalized product
    synergy_normalized = deltaDA_UG * deltaSE_UG,
    
    # Change magnitude
    change_magnitude =  HDRS - baseline_HDRS,
    #sqrt(deltaDA_UG^2 + deltaSE_UG^2),
    
    # Response categories
    response_category = case_when(
      HDRS <= 8 ~ "Remission",
      deltaPerHDRS > 0.5 ~ "Responder", 
      TRUE ~ "Non-Responder"
    ),
    response_category = factor(response_category,
                               levels = c("Non-Responder",
                                          "Responder",
                                          "Remission")),
    
    # Change pattern categories
    change_pattern = case_when(
      deltaDA_UG > 0 & deltaSE_UG > 0 ~ "Both Increase",
      deltaDA_UG < 0 & deltaSE_UG < 0 ~ "Both Decrease",
      deltaDA_UG > 0 & deltaSE_UG < 0 ~ "DA↑/5-HT↓",
      TRUE ~ "DA↓/5-HT↑"
    )
  )

# Calculate effect sizes
model_full <- lm(HDRS ~ deltaDA_UG * deltaSE_UG, data = data)
effect_sizes <- eta_squared(model_full, partial = FALSE)
cohens_f_overall <- cohens_f(model_full)

# For reporting
summary(model_full)

