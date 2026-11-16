# Copyright (C) 2026 Blair Shevlin <blairshevlin@gmail.com>
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
# 2026/11/16    Blair Shevlin                          consolidated all SI analyses into one script; console output only
#
# Section map (SI item -> section of this script -> original source)
#   1. Supp Table 1 (participant NT change, mean (SE))         <- supplementary_analyses.R ("Table of Subject-level NT difference values")
#   2. Supp Table 2 (session-level behaviour, N, Mean (SD))    <- statistical_tests.R ("BEH PROCESSING" + session summaries)
#   3. Supp Tables 6-13 (N-RL-M1..M4, N-UG-M1..M4) and the
#      pairwise contrasts in SI Results "Analysis of trial-level
#      neural data"                                            <- generate_tablesSupplement_ntRegressions.R
#   4. Supp Tables 14-20 (Batten-M1, M3, M3-Stim, M3*, M3*-Stim,
#      M4, M4-Stim)                                            <- generate_tablesSupplement_ntRegressions.R
#   5. Supp Fig 5 legend statistics + SI Results "Alternative
#      normalization" paragraph (task x stimulation ANOVAs)    <- supplementary_analyses.R
#   6. Supp Table 23 (stimulation ANOVAs across metrics/windows)<- supplementary_analyses.R
#   7. Supp Table 24 (within-session performance change)       <- supplementary_analyses.R ("Testing Fatigue")
#   8. Supp Table 25 (baseline performance x HDRS, Pearson)    <- supplementary_analyses.R ("HDRS + Beh")
#   9. SI Results "Exploratory analyses of individual
#      variability" (pattern GLMs, behaviour ANOVAs, medication
#      class logistic regressions) + Supp Table 26             <- supplementary_analyses.R
#  10. Supp Tables 34-35 (MDE, Batten et al. 2024 comparison)  <- power_mde_analysis.R
# Not here: Supp Tables 3-5, 21, 27-33 (main_analyses.R), Supp Table 22 (no code),
# figure-legend statistics produced inside figure scripts (Supp Figs 1-4, 7).
#
# Each section runs inside local() so same-named objects in different sections
# cannot overwrite each other. The script prints to the console and writes no files.

rm(list = ls())

library(tidyverse)
library(fs)
library(here)
library(lmerTest)
library(broom)
library(broom.mixed)
library(emmeans)
library(effectsize)
library(pwr)

options(width = 200, pillar.sigfig = 5)

# Optional: set export_docx = TRUE to also write Supplementary Tables 6-20 as Word tables in
# results/si_tables/ (Supplementary Information format, ready to copy). Off by default: no files written.
export_docx = FALSE
if (export_docx) source(here::here("src", "R", "write_si_table_docx.R"))

# Paths
dir = path(here())
nt_dir = dir / "data" / "nt" / "processed"
beh_dir = dir / "data" / "behavior"
clin_dir = dir / "data" / "clinical"

# Subject-level NT Data
load(nt_dir / "UG_RL_NT-Continuous_9-23-25.RData")

tasks.EST.alt = tasks.EST.alt %>%
    mutate( stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
            task = factor(task, levels = c("RL","UG")))

# Load clinical data
cl.df = read.csv(clin_dir / "clinical-data_deid_07-10-25.csv" )

# Load behavioral data
rl.beh = read.csv(file = beh_dir / "rl-data_deid_07-10-25.csv")

ug.beh = read.csv(file = beh_dir / "ug-data_deid_07-10-25.csv")


##########################################################################
# 1. SUPPLEMENTARY TABLE 1                                               #
# Participant-level NT change (post minus pre), mean (SE)                #
##########################################################################
local({
  cat("\n===== Supplementary Table 1 =====\n")
  cat("diff = Delta (post minus pre mean Oz); se_diff = pooled trial SD / sqrt(n_post + n_pre)\n")

  # Generate summary data - focusing only on reward events (RL) and offer events (UG)
  rl.EST.Reward = rl.EST %>%
    filter(event == "Reward") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))

  ug.EST.Offer = ug.EST %>%
    filter(event == "Offer") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))

  rl.premerge = rl.EST.Reward %>%
    select(idx,stim,nt,Oz) %>%
    mutate(task = "RL")
  ug.premerge = ug.EST.Offer %>%
    select(idx,stim,nt,Oz) %>%
    mutate(task = "UG")
  task.merge = rbind(rl.premerge,ug.premerge)

  # Table of Subject-level NT difference values
  task.merge %>%
    group_by(idx, task, nt) %>%
    summarise(
      diff    = mean(Oz[stim == "Post-Stim"], na.rm = TRUE) - mean(Oz[stim == "Pre-Stim"], na.rm = TRUE),
      sd_post = sd(Oz[stim == "Post-Stim"], na.rm = TRUE),
      sd_pre  = sd(Oz[stim == "Pre-Stim"],  na.rm = TRUE),
      n_post  = sum(!is.na(Oz[stim == "Post-Stim"])),
      n_pre   = sum(!is.na(Oz[stim == "Pre-Stim"])),
      .groups = "drop"
    ) %>%
    mutate(
      sd_pooled = sqrt(((n_post - 1) * sd_post^2 + (n_pre - 1) * sd_pre^2) / (n_post + n_pre - 2)),
      se_diff   = round( sd_pooled / sqrt(n_post + n_pre), 2),
      diff = round(diff,2),
      task = factor(task, levels = c("RL", "UG")),
      nt   = factor(nt,   levels = c("DA", "SE", "NE"))
    ) %>%
    select(idx, task, nt, diff, se_diff) %>%
    arrange(task, nt) %>%
    as.data.frame() %>%
    print()
})


##########################################################################
# 2. SUPPLEMENTARY TABLE 2                                               #
# Session-level behaviour: N, Mean (SD) of participant session means     #
##########################################################################
local({
  cat("\n===== Supplementary Table 2 =====\n")

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

  sess_order = c("Baseline","Pre-Stim","Post-Stim","Week 1",
                 paste("Month", 1:6))

  cat("\nRL: p(Optimal) = mOpt (sdOpt); RT in seconds = meanRT (sdRT)\n")
  rl.beh.means %>% group_by(sess = factor(sess, levels = sess_order)) %>%
    summarise(n = n(),
              mOpt = mean(mChoice),
              sdOpt = sd(mChoice),
              meanRT = mean(mRT),
              sdRT = sd(mRT)) %>%
    print(n = Inf)

  cat("\nUG: p(Accept) = mAccept (sdAccept); RT in seconds = meanRT (sdRT)\n")
  ug.beh.means %>% group_by(sess = factor(sess, levels = sess_order)) %>%
    summarise(n = n(),
              mAccept = mean(mChoice),
              sdAccept = sd(mChoice),
              meanRT = mean(mRT),
              sdRT = sd(mRT)) %>%
    print(n = Inf)

  cat("\nUG: Mood = meanMood (sdMood)\n")
  ug.mood.means %>% group_by(sess = factor(sess, levels = sess_order)) %>%
    summarise(n = n(),
              meanMood = mean(mMood),
              sdMood = sd(mMood)) %>%
    print(n = Inf)
})


##########################################################################
# 3-4. SUPPLEMENTARY TABLES 6-13 (trial-level neural models), THE        #
# PAIRWISE CONTRASTS QUOTED IN THE SI RESULTS, AND SUPPLEMENTARY TABLES  #
# 14-20 (Batten et al. 2024 models).                                     #
##########################################################################
local({
  cl.df.baseline = cl.df %>% filter(session == "pre stim") %>% select(idx,HDRS)

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

  ## RL ##
  ########
  # Session (stim) is coded pre = -1, post = +1 in every model below (set on each model's data frame,
  # because merge() drops contrasts). Contrasts of estimated marginal means do not depend on this coding.
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

  contrasts(rl_trial_data$stim) <- c(-1, 1)  # session coding: pre-stimulation = -1, post-stimulation = +1

  # Model 1: Basic trial-level effects (N-RL-M1, Supp Table 6)
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

  #Model 2: Error and reversal effects (N-RL-M2, Supp Table 7)
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

  #Model 3: Current and previous reward effects (N-RL-M3, Supp Table 8)
  rl.DA.trial.lm3 <- lmer(Oz ~ rew_f * stim +
                      prev_rew_f * stim +
                      rew_f * prev_rew_f * stim +  # interaction between current and previous reward
                      trial_z * stim +
                      (1 + stim | idx),
                    data = rl_trial_data %>% filter(nt == "DA", !is.na(prev_rew)),
                    REML = FALSE,
                    control = lmerControl(optimizer = "nmkbw"))
  # NOTE: unlike DA/NE, the 5-HT model interacts rew_f with numeric prev_rew (not prev_rew_f);
  # kept as in the original script so Supp Table 8 is reproduced as published.
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

  # Model 4: Prediction error effects (N-RL-M4, Supp Table 9)
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

  # Estimated marginal means for the pairwise contrasts

  # Win vs Loss effects by stimulation condition
  rl_da_outcome_effects <- emmeans(rl.DA.trial.lm1, ~ outcome_f * stim)
  rl_se_outcome_effects <- emmeans(rl.SE.trial.lm1, ~ outcome_f * stim)
  rl_ne_outcome_effects <- emmeans(rl.NE.trial.lm1, ~ outcome_f * stim)

  # Block type effects
  rl_da_block_effects <- emmeans(rl.DA.trial.lm1, ~ block_type_factor * stim)
  rl_se_block_effects <- emmeans(rl.SE.trial.lm1, ~ block_type_factor * stim)
  rl_ne_block_effects <- emmeans(rl.NE.trial.lm1, ~ block_type_factor * stim)

  # Slow trial effects
  rl_da_rt_effects <- emmeans(rl.DA.trial.lm2, ~ rt_slow * stim)
  rl_se_rt_effects <- emmeans(rl.SE.trial.lm2, ~ rt_slow * stim)
  rl_ne_rt_effects <- emmeans(rl.NE.trial.lm2, ~ rt_slow * stim)

  # Error trial effects
  rl_da_error_effects <- emmeans(rl.DA.trial.lm2, ~ error_trial * stim)
  rl_se_error_effects <- emmeans(rl.SE.trial.lm2, ~ error_trial * stim)
  rl_ne_error_effects <- emmeans(rl.NE.trial.lm2, ~ error_trial * stim)

  # Reversal trial effects
  rl_da_reversal_effects <- emmeans(rl.DA.trial.lm2, ~ post_reversal * stim)
  rl_se_reversal_effects <- emmeans(rl.SE.trial.lm2, ~ post_reversal * stim)
  rl_ne_reversal_effects <- emmeans(rl.NE.trial.lm2, ~ post_reversal * stim)

  # Reward effects
  rl_da_rew_effects <- emmeans(rl.DA.trial.lm3, ~ prev_rew_f * rew_f * stim)
  rl_se_rew_effects <- emmeans(rl.SE.trial.lm3, ~ prev_rew_f * rew_f * stim)
  rl_ne_rew_effects <- emmeans(rl.NE.trial.lm3, ~ prev_rew_f * rew_f * stim)

  # Prediction error effects
  rl_da_pe_effects <- emmeans(rl.DA.trial.lm4, ~ pred_error_type * stim)
  rl_se_pe_effects <- emmeans(rl.SE.trial.lm4, ~ pred_error_type * stim)
  rl_ne_pe_effects <- emmeans(rl.NE.trial.lm4, ~ pred_error_type * stim)

  ## UG ##
  ########
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

  contrasts(ug_trial_data$stim) <- c(-1, 1)  # session coding: pre-stimulation = -1, post-stimulation = +1

  # Model 1: Basic trial-level effects (N-UG-M1, Supp Table 10)
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

  #Model 2: Choice and response-time effects (N-UG-M2, Supp Table 11)
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

  #Model 3: Current and previous offer effects (N-UG-M3, Supp Table 12)
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

  # Model 4: Prediction error/offer change effects (N-UG-M4, Supp Table 13)
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

  # Offer Bin by stimulation condition
  ug_da_outcome_effects <- emmeans(ug.DA.trial.lm1, ~ offer_bin_f* stim)
  ug_se_outcome_effects <- emmeans(ug.SE.trial.lm1, ~ offer_bin_f * stim)
  ug_ne_outcome_effects <- emmeans(ug.NE.trial.lm1, ~ offer_bin_f * stim)

  # Slow trial effects
  # FIXED 2026/11/16: the original script assigned the SE model to ug_ne_rt_effects and
  # the NE model to ug_se_rt_effects (labels swapped); each object now uses its own model.
  ug_da_rt_effects <- emmeans(ug.DA.trial.lm2, ~ rt_slow * stim)
  ug_se_rt_effects <- emmeans(ug.SE.trial.lm2, ~ rt_slow * stim)
  ug_ne_rt_effects <- emmeans(ug.NE.trial.lm2, ~ rt_slow * stim)

  # Accept trial effects
  ug_da_accept_effects <- emmeans(ug.DA.trial.lm2, ~ accept_trial * stim)
  ug_se_accept_effects <- emmeans(ug.SE.trial.lm2, ~ accept_trial * stim)
  ug_ne_accept_effects <- emmeans(ug.NE.trial.lm2, ~ accept_trial * stim)

  # Offer effects
  ug_da_offer_effects <- emmeans(ug.DA.trial.lm3, ~ prev_offer_bin_f * offer_bin_f * stim)
  ug_se_offer_effects <- emmeans(ug.SE.trial.lm3, ~ prev_offer_bin_f * offer_bin_f * stim)
  ug_ne_offer_effects <- emmeans(ug.NE.trial.lm3, ~ prev_offer_bin_f * offer_bin_f * stim)

  # Offer change effects
  ug_da_oc_effects <- emmeans(ug.DA.trial.lm4, ~ offer_change_type * stim)
  ug_se_oc_effects <- emmeans(ug.SE.trial.lm4, ~ offer_change_type * stim)
  ug_ne_oc_effects <- emmeans(ug.NE.trial.lm4, ~ offer_change_type * stim)

  # Batten et al. (2024) models
  # Model 1: OZ ~ condition x choice + (1 + condition × choice | dataset)
  # Model 3: RZ ~ value + value difference + (1 + value + value difference | dataset)
  # Model 4: RZ ~ value + value difference + absolute value difference + (1 + value + value difference + absolute value difference | dataset)

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
  contrasts(ug_trial_data_n3_star$stim) <- c(-1, 1)  # session coding: pre-stimulation = -1, post-stimulation = +1
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

  # Function to extract and clean fixed effects - flexible for any model structure
  extract_fixed_effects <- function(model, nt_name, task_type, model_type = "basic") {
    tidy(model, effects = "fixed") %>%
      mutate(term = str_replace_all(term, fixed("stim1"), "stimPost-Stim")) %>%  # -1/+1 session term
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
          term == "offer_z:stimPost-Stim" ~ "Offer × Stimulation",

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

          # Batten et al. (2024) models (Supplementary Tables 14-20)
          term == "stimPost-Stim:accept_trial" ~ "Accept × Stimulation",
          term == "offer" ~ "Offer Value",
          term == "offer_change" ~ "Offer Value Difference",
          term == "offer_change_resid" ~ "Offer Value Difference (Residual)",
          term == "abs(offer_change)" ~ "Offer Value Difference (Unsigned)",
          term %in% c("offer:stimPost-Stim", "stimPost-Stim:offer") ~ "Offer Value × Stimulation",
          term %in% c("offer_change:stimPost-Stim", "stimPost-Stim:offer_change") ~ "Value Difference × Stimulation",
          term %in% c("offer_change_resid:stimPost-Stim", "stimPost-Stim:offer_change_resid") ~ "Residual Difference × Stimulation",
          term %in% c("abs(offer_change):stimPost-Stim", "stimPost-Stim:abs(offer_change)") ~ "Unsigned Difference × Stimulation",

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

  # Print Supp Tables 6-13 (SI reports unadjusted P = p_formatted)
  si_nt_tables = tribble(
    ~si, ~task,               ~model_type,           ~label,
    6,   "Reversal learning", "Basic effects",       "N-RL-M1",
    7,   "Reversal learning", "Learning dynamics",   "N-RL-M2",
    8,   "Reversal learning", "Reward interactions", "N-RL-M3",
    9,   "Reversal learning", "Prediction errors",   "N-RL-M4",
    10,  "Ultimatum game",    "Basic effects",       "N-UG-M1",
    11,  "Ultimatum game",    "Context effects",     "N-UG-M2",
    12,  "Ultimatum game",    "Offer interactions",  "N-UG-M3",
    13,  "Ultimatum game",    "Offer dynamics",      "N-UG-M4"
  )
  for (i in seq_len(nrow(si_nt_tables))) {
    cat(sprintf("\n===== Supplementary Table %d (%s) =====\n", si_nt_tables$si[i], si_nt_tables$label[i]))
    si_table = all.lm.res %>%
      filter(task == si_nt_tables$task[i], model_type == si_nt_tables$model_type[i]) %>%
      select(Neurotransmitter = neurotransmitter,
             Effect = term_clean,
             Est. = estimate_formatted,
             SE = se_formatted,
             `t` = t_formatted,
             `p` = p_formatted,
             `p-adj` = p_adj_formatted) %>%
      as.data.frame()
    print(si_table)
    if (export_docx) write_si_table_docx(
      si_table %>% select(Neurotransmitter, Effect, Est., SE, T = t, P = p),
      sprintf("Supplementary Table %d (%s)", si_nt_tables$si[i], si_nt_tables$label[i]),
      sprintf("Supplementary_Table_%02d.docx", si_nt_tables$si[i]))
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
      select(neurotransmitter, task, contrast_type, stim, contrast, contrast_clean, # stim/contrast kept so the session is visible
             estimate_formatted, se_formatted, t_formatted,
             p_formatted, p_adj_formatted, is_significant, is_adj_significant)

    return(contrast_df)
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

  cat("\n===== SI Results, trial-level neural data: pairwise contrasts of estimated marginal means (by session) =====\n")
  cat("p = emmeans (Tukey-adjusted within session); p-adj = additional BH correction within each contrast set\n")
  all_contrasts %>%
    select(neurotransmitter, task, contrast_type, stim, contrast,
           estimate_formatted, se_formatted, t_formatted, p_formatted, p_adj_formatted) %>%
    print()

  # Contrasts conditioned on the previous trial (quoted in SI Results for N-RL-M3 DA and N-UG-M3 NE)
  cat("\n===== SI Results, N-RL-M3 DA: current-reward contrasts by previous reward and session =====\n")
  print(pairs(rl_da_rew_effects, by = c("prev_rew_f","stim")))
  cat("\n===== SI Results, N-UG-M3 NE: current-offer contrasts by previous offer and session =====\n")
  print(pairs(ug_ne_offer_effects, by = c("prev_offer_bin_f","stim")))

  # Print Supp Tables 14-20 (Batten et al. 2024 models)
  batten_tables = list(
    "14 (Batten-M1)"       = list(Dopamine = ug.DA.trial.s1,           Serotonin = ug.SE.trial.s1,           Norepinephrine = ug.NE.trial.s1),
    "15 (Batten-M3)"       = list(Dopamine = ug.DA.trial.s3,           Serotonin = ug.SE.trial.s3,           Norepinephrine = ug.NE.trial.s3),
    "16 (Batten-M3-Stim)"  = list(Dopamine = ug.DA.trial.s3.stim,      Serotonin = ug.SE.trial.s3.stim,      Norepinephrine = ug.NE.trial.s3.stim),
    "17 (Batten-M3*)"      = list(Dopamine = ug.DA.trial.s3_star,      Serotonin = ug.SE.trial.s3_star,      Norepinephrine = ug.NE.trial.s3_star),
    "18 (Batten-M3*-Stim)" = list(Dopamine = ug.DA.trial.s3_star.stim, Serotonin = ug.SE.trial.s3_star.stim, Norepinephrine = ug.NE.trial.s3_star.stim),
    "19 (Batten-M4)"       = list(Dopamine = ug.DA.trial.s4,           Serotonin = ug.SE.trial.s4,           Norepinephrine = ug.NE.trial.s4),
    "20 (Batten-M4-Stim)"  = list(Dopamine = ug.DA.trial.s4.stim,      Serotonin = ug.SE.trial.s4.stim,      Norepinephrine = ug.NE.trial.s4.stim)
  )
  for (tab in names(batten_tables)) {
    cat(sprintf("\n===== Supplementary Table %s =====\n", tab))
    map_dfr(batten_tables[[tab]], ~ .x %>%
              tidy(effects = "fixed") %>%
              mutate(p_adjusted = p.adjust(p.value, method = "BH") ),
            .id = "neurotransmitter") %>%
      select(neurotransmitter, term, estimate, std.error, statistic, df, p.value, p_adjusted) %>%
      as.data.frame() %>%
      print(digits = 3)
    if (export_docx) write_si_table_docx(
      imap_dfr(batten_tables[[tab]], ~ extract_fixed_effects(.x, .y, "Ultimatum game", tab)) %>%
        transmute(Neurotransmitter = neurotransmitter, Effect = term_clean,
                  Est. = sprintf("%.2f", as.numeric(estimate_formatted)), SE = se_formatted,
                  T = t_formatted, P = p_formatted),
      paste("Supplementary Table", tab),
      sprintf("Supplementary_Table_%s.docx", sub(" .*", "", tab)))
  }
})


##################################################################
# 5. SUPPLEMENTARY FIGURE 5 STATISTICS + SI RESULTS "ALTERNATIVE #
# NORMALIZATION" PARAGRAPH                                       #
# Cross-task ("shared") z-scoring, Oz metric, 500 ms post-event. #
##################################################################
local({
  cat("\n===== Supplementary Figure 5 statistics (and SI Results, alternative normalization) =====\n")

  shared = function(task_name) tasks.EST.alt %>%
    filter(task == task_name, event == "Shared") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))

  rl.OR = rl.beh %>% filter(sess %in% c("Pre-Stim","Post-Stim")) %>%
    mutate(stim = factor(sess, levels = c("Pre-Stim","Post-Stim")))
  ug.OR = ug.beh %>% filter(sess %in% c("Pre-Stim","Post-Stim")) %>%
    mutate(stim = factor(sess, levels = c("Pre-Stim","Post-Stim")),
           offer_bin = ifelse(offer < 4, "Low", ifelse(offer > 6, "High", "Middle")))

  rl = merge(shared("RL"), rl.OR, by = c("idx","stim","trial","rt")) %>%
    mutate(block_type = factor(cond, levels = c("Mixed","Negative","Positive")))
  ug = merge(shared("UG"), ug.OR, by = c("idx","stim","trial","rt")) %>%
    mutate(offer_bin_f = factor(offer_bin, levels = c("Middle","Low","High")))
  contrasts(rl$stim) <- c(-1,1)
  contrasts(ug$stim) <- c(-1,1)

  # Panels c-f: stimulation x block type (RL) and stimulation x offer amount (UG), participant-level means
  rl.aov = rl %>% group_by(idx,stim,nt,block_type) %>% summarise(oz = mean(Oz), .groups = "drop")
  ug.aov = ug %>% group_by(idx,stim,nt,offer_bin_f) %>% summarise(oz = mean(Oz), .groups = "drop")
  for (n in c("DA","SE")) {
    cat("\nSupp Fig 5, RL", n, "\n"); print(summary(aov(oz ~ stim * block_type, data = filter(rl.aov, nt == n))))
    cat("\nSupp Fig 5, UG", n, "\n"); print(summary(aov(oz ~ stim * offer_bin_f, data = filter(ug.aov, nt == n))))
  }

  # Task x stimulation interaction across both tasks
  combo = bind_rows(
    RL = rl %>% group_by(idx,stim,nt) %>% summarise(oz = mean(Oz), .groups = "drop"),
    UG = ug %>% group_by(idx,stim,nt) %>% summarise(oz = mean(Oz), .groups = "drop"),
    .id = "task") %>%
    mutate(task = factor(task, levels = c("RL","UG")))
  for (n in c("DA","SE")) {
    cat("\nSupp Fig 5 text, task x stimulation,", n, "\n"); print(summary(aov(oz ~ stim * task, data = filter(combo, nt == n))))
  }
})


##########################################################################
# 6. SUPPLEMENTARY TABLE 23                                              #
# Stimulation ANOVAs across quantification metrics and windows           #
##########################################################################
local({
  cat("\n===== Supplementary Table 23 =====\n")

  rl.EST.Reward = rl.EST %>%
    filter(event == "Reward") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))

  ug.EST.Offer = ug.EST %>%
    filter(event == "Offer") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))

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

  # Set up contrast coding for stim
  contrasts(rl.EST.pr$stim) <- c(-1,1)
  contrasts(ug.EST.pr$stim) <- c(-1,1)

  # ANOVA
  rl.EST.aov = rl.EST.pr %>% group_by(idx,stim,nt,block_type) %>%
    summarise(oz = mean(Oz), mz = mean(Mz), pz = mean(Pz), rz = mean(Rz),
              oz10 = mean(Oz10), mz10 = mean(Mz10), pz10 = mean(Pz10), rz10 = mean(Rz10))

  ug.EST.aov = ug.EST.pr %>% group_by(idx,stim,nt,offer_bin_f) %>%
    summarise(oz = mean(Oz), mz = mean(Mz), pz = mean(Pz), rz = mean(Rz),
              oz10 = mean(Oz10), mz10 = mean(Mz10), pz10 = mean(Pz10), rz10 = mean(Rz10))

  # Define metrics and neurotransmitter types
  metrics <- c("oz", "mz", "pz", "rz","oz10", "mz10", "pz10", "rz10")
  nt_types <- c("DA", "SE", "NE")

  # Function to run ANOVA for a given dataset, metric, and NT type
  run_anova <- function(data, metric, nt_type, formula_str) {
    # Filter data for specific NT type
    filtered_data <- data[data$nt == nt_type, ]

    # Create formula with the specific metric
    formula_obj <- as.formula(paste(metric, formula_str))

    # Run ANOVA
    aov_result <- aov(formula_obj, data = filtered_data)

    return(aov_result)
  }

  rl_combinations <- expand.grid(
    metric = metrics,
    nt_type = nt_types,
    stringsAsFactors = FALSE
  )
  ug_combinations <- expand.grid(
    metric = metrics,
    nt_type = nt_types,
    stringsAsFactors = FALSE
  )

  rl_anovas <- map2(
    rl_combinations$metric,
    rl_combinations$nt_type,
    ~ run_anova(rl.EST.aov, .x, .y, "~ stim * block_type")
  )
  ug_anovas <- map2(
    ug_combinations$metric,
    ug_combinations$nt_type,
    ~ run_anova(ug.EST.aov, .x, .y, "~ stim * offer_bin_f")
  )

  names(rl_anovas) <- paste("rl", rl_combinations$nt_type, rl_combinations$metric, sep = " ")
  names(ug_anovas) <- paste("ug", ug_combinations$nt_type, ug_combinations$metric, sep = " ")

  extract_anova_summary <- function(anova_result, task, nt_type, metric) {
    summary_tidy <- tidy(anova_result)
    summary_tidy$task <- task
    summary_tidy$nt <- nt_type
    summary_tidy$metric <- metric
    return(summary_tidy)
  }

  # Extract summaries for all ANOVAs
  rl_summaries <- pmap_dfr(
    list(
      rl_anovas,
      rep("RL", length(rl_anovas)),
      rep(rl_combinations$nt_type, 1),
      rep(rl_combinations$metric, 1)
    ),
    extract_anova_summary
  )

  ug_summaries <- pmap_dfr(
    list(
      ug_anovas,
      rep("UG", length(ug_anovas)),
      rep(ug_combinations$nt_type, 1),
      rep(ug_combinations$metric, 1)
    ),
    extract_anova_summary
  )

  anova_summaries <- bind_rows(rl_summaries, ug_summaries)

  # Format the ANOVA summary data
  format_anova_results <- function(anova_summaries) {

    # Define metric labels for clarity
    metric_labels <- c(
      "oz" = "Sum (500ms)",
      "mz" = "Mean (500ms)",
      "pz" = "Peak (500ms)",
      "rz" = "Relative (500ms)",
      "oz10" = "Sum (1000ms)",
      "mz10" = "Mean (1000ms)",
      "pz10" = "Peak (1000ms)",
      "rz10" = "Relative (1000ms)"
    )

    # Process the data
    formatted_data <- anova_summaries %>%
      filter(term == "stim") %>%  # Focus on stimulation main effect
      mutate(
        # Rename metrics for clarity
        metric_label = recode(metric, !!!metric_labels),
        # Format p-values with significance stars
        p_formatted = case_when(
          p.value < 0.001 ~ sprintf("%.4f***", p.value),
          p.value < 0.01 ~ sprintf("%.3f**", p.value),
          p.value < 0.05 ~ sprintf("%.3f*", p.value),
          TRUE ~ sprintf("%.3f", p.value)
        ),
        # Format F-statistics
        f_formatted = sprintf("%.2f", statistic),
        # Organize neurotransmitter names
        nt_full = case_when(
          nt == "DA" ~ "Dopamine",
          nt == "SE" ~ "Serotonin",
          nt == "NE" ~ "Norepinephrine",
          TRUE ~ nt
        ),
        # Organize task names
        task_full = case_when(
          task == "RL" ~ "Reversal Learning",
          task == "UG" ~ "Ultimatum Game",
          TRUE ~ task
        )
      ) %>%
      select(
        Task = task_full,
        Neurotransmitter = nt_full,
        Metric = metric_label,
        F = f_formatted,
        df1 = df,
        p = p_formatted
      )

    return(formatted_data)
  }

  print(format_anova_results(anova_summaries), n = Inf)

  cat("\nResidual df (Supp Table 23 note)\n")
  anova_summaries %>% filter(term == "Residuals") %>% distinct(task, df) %>% print()
})


##########################################################################
# 7. SUPPLEMENTARY TABLE 24                                              #
# Within-session performance change (trial slopes per session)          #
##########################################################################
local({
  cat("\n===== Supplementary Table 24 =====\n")

  cl.df.baseline = cl.df %>% filter(session == "pre stim") %>% select(idx,HDRS)

  # Generate summary data - focusing only on reward events (RL) and offer events (UG)
  rl.EST.Reward = rl.EST %>%
    filter(event == "Reward") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))

  ug.EST.Offer = ug.EST %>%
    filter(event == "Offer") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))

  rl.EST.pr = merge(rl.EST.Reward,cl.df.baseline, by = "idx") %>%
    filter(rt < 10) %>%
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
           HDRS_z = scale(HDRS)[,1],
            trial_z = scale(trial)[,1]
           )

  ug.EST.pr =  merge(ug.EST.Offer,cl.df.baseline, by = "idx") %>%
    filter(rt < 10) %>%
    mutate(offer_bin_f = factor(offer_bin, levels = c("Middle","Low","High")),
           offer_z = scale(offer)[,1],
           offer_change = offer - prev_offer_raw,
           offer_change_z = scale(offer_change)[,1],
           offer_change_f = factor(ifelse(offer_change > 0, "improve","worse"),
                                   levels=c("worse","improve")),
           HDRS_z = scale(HDRS)[,1],
           trial_z = scale(trial)[,1]
    )

  # Set up contrast coding for stim
  contrasts(rl.EST.pr$stim) <- c(-1,1)
  contrasts(ug.EST.pr$stim) <- c(-1,1)

  # Reversal Learning RT analysis
  rl_rt_analysis = rl.EST.pr %>%
    mutate(trial_z = scale(trial)[,1],
           log_rt = log(rt)) %>%
    group_by(idx, stim, trial, trial_z, block_type) %>%
    summarise(mean_rt = mean(rt, na.rm = TRUE),
              log_rt = mean(log_rt, na.rm = TRUE))

  # Ultimatum Game RT analysis
  ug_rt_analysis = ug.EST.pr %>%
    mutate(trial_z = scale(trial)[,1],
           log_rt = log(rt)) %>%
    group_by(idx, stim, trial, trial_z, offer_bin) %>%
    summarise(mean_rt = mean(rt, na.rm = TRUE),
              log_rt = mean(log_rt, na.rm = TRUE))

  contrasts(rl_rt_analysis$stim) <- c(-1, 1)  # session coding: pre-stimulation = -1, post-stimulation = +1
  contrasts(ug_rt_analysis$stim) <- c(-1, 1)  # session coding: pre-stimulation = -1, post-stimulation = +1

  # Reversal Learning RT slopes
  rl_rt_model = lmer(log_rt ~ trial_z * stim + (1 + stim|idx),
                     data = rl_rt_analysis,
                     REML = FALSE)

  # Ultimatum Game RT slopes
  ug_rt_model = lmer(log_rt ~ trial_z * stim + (1 + stim|idx),
                     data = ug_rt_analysis,
                     REML = FALSE)

  # Choice accuracy analysis (for reversal learning)
  rl_choice_analysis = rl.EST.pr %>%
    mutate(trial_z = scale(trial)[,1]) %>%
    group_by(idx, stim, trial, trial_z, rew) %>%
    summarise(optimal_choice = mean(opt, na.rm = TRUE),
              .groups = 'drop')

  contrasts(rl_choice_analysis$stim) <- c(-1, 1)  # session coding: pre-stimulation = -1, post-stimulation = +1

  # Test choice accuracy over time
  rl_choice_model = glmer(optimal_choice ~ trial_z * stim + (1 + stim|idx),
                          data = rl_choice_analysis,
                          family = binomial,
                          control = glmerControl(optimizer = "bobyqa"))

  # Acceptance rate analysis (for ultimatum game)
  ug_choice_analysis = ug.EST.pr %>%
    mutate(trial_z = scale(trial)[,1],
           offer_z = scale(offer)[,1],
           accept = 1*(rej==0)) %>%
    group_by(idx, stim, trial, trial_z, offer_z) %>%
    summarise(accept_rate = mean(accept, na.rm = TRUE),
              .groups = 'drop')

  contrasts(ug_choice_analysis$stim) <- c(-1, 1)  # session coding: pre-stimulation = -1, post-stimulation = +1

  # Test acceptance rate over time
  ug_choice_model = glmer(accept_rate ~  trial_z * stim + (1 + stim|idx),
                          data = ug_choice_analysis,
                          family = binomial,
                          control = glmerControl(optimizer = "bobyqa"))

  # Extract slope estimates
  rl_rt_slopes <- emtrends(rl_rt_model, ~ stim, var = "trial_z")
  ug_rt_slopes <- emtrends(ug_rt_model, ~ stim, var = "trial_z")
  rl_c_slopes <- emtrends(rl_choice_model, ~ stim, var = "trial_z")
  ug_c_slopes <- emtrends(ug_choice_model, ~ stim, var = "trial_z")

  # Test if slopes are significantly different from zero
  cat("\nRL, response time (log s)\n");  print(test(rl_rt_slopes, null = 0))
  cat("\nRL, choice quality (optimal)\n"); print(test(rl_c_slopes, null = 0))
  cat("\nUG, response time (log s)\n");  print(test(ug_rt_slopes, null = 0))
  cat("\nUG, choice quality (accept)\n"); print(test(ug_c_slopes, null = 0))
})


##########################################################################
# 8. SUPPLEMENTARY TABLE 25                                              #
# Pearson correlations, baseline (fMRI session) performance x HDRS       #
##########################################################################
local({
  cat("\n===== Supplementary Table 25 =====\n")

  cl.df.baseline = cl.df %>% filter(session == "pre stim") %>% select(idx,HDRS)
  cl.df.baseline$HDRS_z = scale(cl.df.baseline$HDRS)[,1]

  # Baseline session
  rl.task.baseline =
      rl.beh %>%
      filter(sess %in% c("fMRI"), rt < 10) %>%
      select(!c(session,sess,score))

  ug.task.baseline =
      ug.beh %>%
      filter(sess %in% c("fMRI"), rt < 10) %>%
      select(!c(session,sess))

  ug.cl.baseline = merge(ug.task.baseline, cl.df.baseline) %>% mutate(offer_z = scale(offer)[,1])
  rl.cl.baseline = merge(rl.task.baseline, cl.df.baseline)

  ug.cl.baseline.mean = ug.cl.baseline %>%
    group_by(idx, HDRS,HDRS_z) %>%
    summarise(mAccept = mean(rej==0),
              mRT = mean(log(rt)),
              mMood = mean(mood, na.rm = T))

  rl.cl.baseline.mean = rl.cl.baseline %>%
    group_by(idx,HDRS,HDRS_z) %>%
    summarise(mOpt = mean(opt==1),
              mRT = mean(log(rt)))

  # RL #
  cat("\nRL, response time\n")
  print(cor.test(rl.cl.baseline.mean$HDRS_z,
                     rl.cl.baseline.mean$mRT,
                     method = "pearson"))
  cat("\nRL, choice quality\n")
  print(cor.test(rl.cl.baseline.mean$HDRS_z,
                     rl.cl.baseline.mean$mOpt,
                     method = "pearson"))

  # UG #
  cat("\nUG, response time\n")
  # FIXED 2026/11/16: original passed `pmethod = "pearson"` (typo, silently ignored);
  # Pearson is cor.test()'s default, so the result is unchanged.
  print(cor.test(ug.cl.baseline.mean$HDRS_z,
                     ug.cl.baseline.mean$mRT,
                     method = "pearson"))
  cat("\nUG, choice quality\n")
  print(cor.test(ug.cl.baseline.mean$HDRS_z,
                     ug.cl.baseline.mean$mAccept,
                     method = "pearson"))
  cat("\nUG, mood\n")
  print(cor.test(ug.cl.baseline.mean$HDRS_z,
                     ug.cl.baseline.mean$mMood,
                     method = "pearson"))
})


##########################################################################
# 9. SI RESULTS "EXPLORATORY ANALYSES OF INDIVIDUAL VARIABILITY IN       #
# NEUROCHEMICAL RESPONSES" + SUPPLEMENTARY TABLE 26                      #
##########################################################################
local({
  cl.df.baseline = cl.df %>% filter(session == "pre stim") %>% select(idx,HDRS)
  cl.df.m6 = cl.df %>% filter(session == "month 6") %>% select(idx,HDRS)

  # Prepare data
  rl.EST.Reward = rl.EST %>%
    filter(event == "Reward") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))

  ug.EST.Offer = ug.EST %>%
    filter(event == "Offer") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))

  ug.EST.pr = ug.EST.Offer %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
           offer_bin_f = factor(offer_bin, levels = c("Middle","Low","High")),
           offer_z = scale(offer)[,1],
           offer_change = offer - prev_offer_raw,
           offer_change_z = scale(offer_change)[,1],
           offer_change_f = factor(ifelse(offer_change > 0, "improve","worse"),
                                   levels=c("worse","improve")),
           trial_z = scale(trial)[,1]
    )

  rl.EST.pr = rl.EST.Reward %>%
      mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
           block_type = factor(cond, levels = c("Mixed","Negative","Positive")),
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

  rl.premerge = rl.EST.pr %>%
    select(idx,stim,nt,Oz,O,trial_z) %>%
    mutate(task = "RL")
  ug.premerge = ug.EST.pr %>%
    select(idx,stim,nt,Oz,O,trial_z) %>%
    mutate(task = "UG")

  task.merge = rbind(rl.premerge,ug.premerge) %>%
              mutate(task = factor(task))

  task.merge.baseline = merge(task.merge, cl.df.baseline)
  task.merge.baseline$HDRS_z = scale(task.merge.baseline$HDRS)[,1]

  # What differentiates participants who show the expected pattern vs. those who do not?
  task.means = task.merge.baseline %>%
    group_by(idx, task, nt, stim, HDRS_z, HDRS) %>%
    summarize(Oz_mean = mean(Oz, na.rm = TRUE), O_mean = mean(O, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(HDRS_M6 = cl.df.m6$HDRS[match(idx, cl.df.m6$idx)],
           HDRS_M6_Change = HDRS_M6 - HDRS,
           HDRS_M6_z = scale(HDRS_M6)[,1])

  # Pull out the post-stim Oz_mean for each subject
  # (pattern = sign of the post-stimulation mean, not of the post-minus-pre change)
  post_stim_vals <- task.means %>%
    filter(stim == "Post-Stim") %>%
    select(idx, task, nt, Oz_mean) %>%
    mutate(change_pattern = ifelse(Oz_mean > 0, "Increased", "Decreased")) %>%
    group_by(idx, task) %>%
    mutate(expected_pattern = ifelse(
      any(nt == "DA" & change_pattern == "Increased") &
      any(nt == "SE" & change_pattern == "Increased"),
      "Expected", "Unexpected"
    ),
    DA_change = change_pattern[nt == "DA"],
    SE_change = change_pattern[nt == "SE"]
    ) %>%
    ungroup() %>% select(!Oz_mean)

  # Join and apply logic
  task.means <- task.means %>%
    left_join(post_stim_vals, by = c("idx", "task", "nt"))

  # Medication list (pre-operative prescriptions)
  cl.alt = data.frame(idx = c(1:10),
    meds = c(NA,
            paste0("quetiapine,", "propranolol,", "losartan,", "amlodipine,", "hydrochlorothiazide,","gabapentin,"),
            paste0("fetzima,","lamotrigine,","allopurinol,","lisinopril"),
            paste0("desvenlafaxine,","dextroamphetamine,","buprenorphine,","lamotrigine,","mirtazapine,","diazepam,","levothyroxine,","semaglutide,","testosterone,"),
            paste0("wellbutrin,","junel fe,","levothyroxine,","mirtazapine,","lithium,","timolol,","aripiprazole,"),
            paste0("wellbutrin,","desvenlafaxine,","NAC powder,","methylphenidate,","diazepam,","hydroxyzine,","levothyroxine,","liothyronine,"),
            paste0("venlafaxine,","modafinil,","clonazepam,","lisinopril-hydrochlorothiazide,"),
            paste0("venlafaxine,","lisdexamfetamine,"),
            paste0("nortriptyline,","lamotrigine,","lithium carbonate,","doxepin,", "busprione,", "hydroxyzine,", "nadolol,"),
            paste0("mirtazapine,","clonazepam,","quetiapine,","crestor,","metoprolol,")
            )
  )

  cat("\n===== Supplementary Table 26 =====\n")
  cat("UG post-stimulation DA/5-HT pattern (sign of post-stim mean, used for 'expected pattern'), RL pattern, Month-6 HDRS change, medication list\n")
  post_stim_vals %>%
    filter(task == "UG", nt == "DA") %>%
    select(idx, DA_change, SE_change, expected_pattern) %>%
    left_join(post_stim_vals %>% filter(task == "RL", nt == "DA") %>%
                select(idx, RL_expected_pattern = expected_pattern), by = "idx") %>%
    left_join(task.means %>% distinct(idx, HDRS, HDRS_M6, HDRS_M6_Change), by = "idx") %>%
    left_join(cl.alt, by = "idx") %>%
    arrange(idx) %>%
    as.data.frame() %>%
    print()

  cat("\nPer-participant post-minus-pre Oz change (DA, 5-HT) with Month-6 HDRS change\n")
  task.means %>% select(idx, task, nt, stim, Oz_mean, HDRS_M6_Change) %>% distinct() %>%
    filter(nt!="NE") %>%
    group_by(idx, task, nt) %>%
    mutate(Post_minus_Pre = Oz_mean[stim == "Post-Stim"] - Oz_mean[stim == "Pre-Stim"]) %>%
    select(!c(Oz_mean,stim)) %>%
    distinct() %>%
    pivot_wider(names_from = c(task, nt), values_from = Post_minus_Pre, names_prefix = "Oz_") %>%
    arrange(idx) %>%
    print(n = Inf)

  cat("\n===== SI Results, expected vs unexpected pattern (UG): baseline HDRS, Gaussian GLM =====\n")
  task.means %>%
    select(idx, task, expected_pattern, HDRS) %>%
    mutate(expected_pattern = factor(expected_pattern, levels = c("Unexpected", "Expected"))) %>%
    distinct() %>%
    filter(task == "UG") %>%
    glm(HDRS ~ expected_pattern, data = ., family = gaussian) %>%
    summary() %>% print()

  cat("\n===== SI Results, expected vs unexpected pattern (UG): Month-6 HDRS change, Gaussian GLM =====\n")
  # Reference level here is "Expected" (character, alphabetical), so the coefficient is
  # Unexpected - Expected; the SI reports the same effect as Expected - Unexpected (sign flipped).
  task.means %>%
    select(idx, task, expected_pattern, HDRS, HDRS_M6) %>%
    distinct() %>%
    mutate(HDRS_change = HDRS_M6 - HDRS) %>%
    filter(task == "UG") %>%
    glm(HDRS_change ~ expected_pattern, data = ., family = gaussian) %>%
    summary() %>% print()

  # Look at behavior
  ug.EST.pr = merge(ug.EST.pr, post_stim_vals[post_stim_vals$task == "UG",], by = c("idx", "nt"))
  rl.EST.pr = merge(rl.EST.pr, post_stim_vals[post_stim_vals$task == "RL",], by = c("idx", "nt"))

  ug.beh = ug.EST.pr %>%
    filter(nt == "DA") %>%
    group_by(idx,expected_pattern,nt,stim) %>%
    summarise(mAccept = mean(rej==0), mRT = mean(log(rt)))

  # NOTE: mOpt is the proportion of NON-optimal choices (opt==0) as in the original script;
  # for this ANOVA it is a linear transform (1 - p(optimal)) and gives identical F and P.
  rl.beh = rl.EST.pr %>%
    filter(nt == "DA") %>%
    group_by(idx,expected_pattern,nt,stim) %>%
    summarise(mOpt = mean(opt==0), mRT = mean(log(rt)))

  cat("\n===== SI Results, pattern x session ANOVAs on behaviour =====\n")
  cat("\nUG choice (acceptance)\n");  print(summary(aov(mAccept ~ expected_pattern * stim, data = ug.beh)))
  cat("\nUG response time (log)\n");  print(summary(aov(mRT ~ expected_pattern * stim, data = ug.beh)))
  cat("\nRL choice\n");               print(summary(aov(mOpt ~ expected_pattern * stim, data = rl.beh)))
  cat("\nRL response time (log)\n");  print(summary(aov(mRT ~ expected_pattern * stim, data = rl.beh)))

  # Medication classes
  # Create binary columns for each drug
  cl.alt.binary <- cl.alt %>%
    separate_rows(meds, sep = ",") %>%
    mutate(meds = trimws(meds)) %>%
    filter(meds != "" & !is.na(meds)) %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = meds, values_from = value, values_fill = 0) %>%
    right_join(cl.alt %>% select(idx), by = "idx") %>%
    mutate(across(-idx, ~replace_na(., 0))) %>%
    arrange(idx)

  # Let's categorize medications by class and see if any classes are associated with patterns of change
  drug_classes <- list(
    SNRI            = c("venlafaxine", "desvenlafaxine", "fetzima"),
    NDRI            = c("wellbutrin"),
    NaSSA           = c("mirtazapine"),
    TCA             = c("nortriptyline", "doxepin"),
    Antipsychotic   = c("quetiapine", "aripiprazole"),
    Mood_Stabilizer = c("lamotrigine", "lithium", "lithium carbonate"),
    Benzo           = c("diazepam", "clonazepam"),
    Stimulant       = c("methylphenidate", "dextroamphetamine", "lisdexamfetamine", "modafinil"),
    Anticonvulsant  = c("gabapentin"),
    Thyroid         = c("levothyroxine", "liothyronine"),
    Anxiolytic      = c("hydroxyzine", "busprione"),
    Beta_Blocker    = c("propranolol", "metoprolol", "nadolol", "timolol"),
    ACE_ARB         = c("lisinopril", "losartan", "lisinopril-hydrochlorothiazide"),
    Opioid_Adj      = c("buprenorphine"),
    Other_Psych     = c("NAC powder"),
    Hormonal        = c("testosterone", "junel fe"),
    Statin          = c("crestor"),
    Diuretic        = c("hydrochlorothiazide"),
    CCB             = c("amlodipine"),
    GLP1            = c("semaglutide"),
    Xanthine_Ox_Inh = c("allopurinol")
  )

  # Build a class-level binary dataframe
  cl.alt.class <- cl.alt.binary %>%
    mutate(
      SNRI            = as.integer(rowSums(across(any_of(drug_classes$SNRI))) > 0),
      NDRI            = as.integer(rowSums(across(any_of(drug_classes$NDRI))) > 0),
      NaSSA           = as.integer(rowSums(across(any_of(drug_classes$NaSSA))) > 0),
      TCA             = as.integer(rowSums(across(any_of(drug_classes$TCA))) > 0),
      Antipsychotic   = as.integer(rowSums(across(any_of(drug_classes$Antipsychotic))) > 0),
      Mood_Stabilizer = as.integer(rowSums(across(any_of(drug_classes$Mood_Stabilizer))) > 0),
      Benzo           = as.integer(rowSums(across(any_of(drug_classes$Benzo))) > 0),
      Stimulant       = as.integer(rowSums(across(any_of(drug_classes$Stimulant))) > 0),
      Anticonvulsant  = as.integer(rowSums(across(any_of(drug_classes$Anticonvulsant))) > 0),
      Thyroid         = as.integer(rowSums(across(any_of(drug_classes$Thyroid))) > 0),
      Anxiolytic      = as.integer(rowSums(across(any_of(drug_classes$Anxiolytic))) > 0),
      Beta_Blocker    = as.integer(rowSums(across(any_of(drug_classes$Beta_Blocker))) > 0),
      ACE_ARB         = as.integer(rowSums(across(any_of(drug_classes$ACE_ARB))) > 0)
    )

  # Check class counts
  class_counts <- cl.alt.class %>%
    select(SNRI:ACE_ARB) %>%
    summarise(across(everything(), sum)) %>%
    pivot_longer(everything(), names_to = "class", values_to = "n_participants") %>%
    arrange(desc(n_participants))

  cat("\n===== SI Results, medication classes: participants per class =====\n")
  print(class_counts, n = Inf)

  # Define classes with enough participants (n > 1, ideally > 2)
  multi_classes <- class_counts %>%
    filter(n_participants > 1) %>%
    pull(class)

  # Join with stim data
  cl.alt.stim_class <- post_stim_vals %>%
    select(idx, task, expected_pattern) %>%
    left_join(cl.alt.class, by = "idx") %>%
    mutate(expected_pattern_binary = ifelse(expected_pattern == "Expected", 1, 0))

  # --- Expected pattern ~ drug class ---
  class_models_pattern <- list()
  for (task in c("UG", "RL")) {
    for (drug_class in multi_classes) {
      model <- glm(expected_pattern_binary ~ get(drug_class),
                   data = cl.alt.stim_class[cl.alt.stim_class$task == task, ],
                   family = binomial)
      class_models_pattern[[paste0(task, "_", drug_class)]] <- summary(model)
    }
  }

  extract_glm_results <- function(model_list, nt_var = NULL) {
    map_dfr(names(model_list), function(name) {
      coefs <- model_list[[name]]$coefficients

      # Handle models where predictor may have failed (e.g., perfect separation)
      if (nrow(coefs) < 2) {
        return(tibble(
          model     = name,
          task      = str_split(name, "_")[[1]][1],
          predictor = str_split(name, "_", n = 2)[[1]][2],
          nt        = if (!is.null(nt_var)) nt_var else NA_character_,
          estimate  = NA_real_, se = NA_real_, z = NA_real_,
          p_value   = NA_real_
        ))
      }

      pred_row <- coefs[2, , drop = FALSE]
      tibble(
        model     = name,
        task      = str_split(name, "_")[[1]][1],
        predictor = str_split(name, "_", n = 2)[[1]][2],
        nt        = if (!is.null(nt_var)) nt_var else NA_character_,
        estimate  = pred_row[, "Estimate"],
        se        = pred_row[, "Std. Error"],
        z         = pred_row[, "z value"],
        p_value   = pred_row[, "Pr(>|z|)"]
      )
    }) %>%
    mutate(p_holm = p.adjust(p_value, method = "holm")) %>%
    arrange(p_value)
  }

  cat("\n===== SI Results, logistic regressions: expected pattern ~ medication class (per task) =====\n")
  print(extract_glm_results(class_models_pattern), n = Inf)
})


##########################################################################
# 10. SUPPLEMENTARY TABLES 34-35                                         #
# Minimum detectable effects and Batten et al. (2024) comparison         #
##########################################################################
local({
  alpha_uncorrected = 0.05
  alpha_bonferroni  = 0.05 / 3   # Bonferroni within the three primary tests
  target_power = 0.80

  # MDE for tests (a)/(b): paired pre/post comparisons of subject-level means
  mde_paired = function(n, alpha = alpha_uncorrected, power = target_power) {
    fit = pwr.t.test(n = n, sig.level = alpha, power = power, type = "paired",
                      alternative = "two.sided")
    tibble(n = n, alpha = alpha, power = power, mde_dz = fit$d)
  }

  primary_ab_mde = bind_rows(
    mde_paired(n = 8)  %>% mutate(test = "a_DA_RL_task_specificity",   design = "paired, n=8 (RL)"),
    mde_paired(n = 10) %>% mutate(test = "b_5HT_UG_task_specificity",  design = "paired, n=10 (UG)")
  )

  cat("\n===== Supplementary Table 34: MDE, DA x RL and 5-HT x UG (paired) =====\n")
  print(primary_ab_mde)

  # MDE for test (c): interaction term of HDRS ~ deltaDA_UG * deltaSE_UG (u = 1, v = n - 3 - 1)
  assumed_v_test_c = 6   # n=10 - 3 predictors - 1 = 6; re-verified below

  mde_interaction = function(u, v, alpha, power = target_power) {
    fit = pwr.f2.test(u = u, v = v, sig.level = alpha, power = power)
    tibble(u = u, v = v, alpha = alpha, power = power,
           mde_f2 = fit$f2, mde_partial_R2 = fit$f2 / (1 + fit$f2))
  }

  primary_c_mde = bind_rows(
    mde_interaction(u = 1, v = assumed_v_test_c, alpha = alpha_bonferroni)  %>% mutate(correction = "Bonferroni-within-3 (matches main_analyses.R)"),
    mde_interaction(u = 1, v = assumed_v_test_c, alpha = alpha_uncorrected) %>% mutate(correction = "uncorrected, for reference only")
  ) %>%
    mutate(test = "c_DAx5HT_HDRS_key_interaction")

  cat("\n===== Supplementary Table 34: MDE, DA x 5-HT -> Month 6 HDRS interaction =====\n")
  print(primary_c_mde)

  # Observed effect sizes, for comparison against the MDEs
  rl.EST.Reward = rl.EST %>%
    filter(event == "Reward") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))
  ug.EST.Offer = ug.EST %>%
    filter(event == "Offer") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))

  # One-sample Cohen's d (mu = 0) on per-subject change scores, CI from the noncentral t
  observed_d = function(est_data, nt_filter) {
    change_data = est_data %>%
      mutate(nt = factor(nt, levels = c("DA","SE","NE"))) %>%
      filter(nt == nt_filter) %>%
      select(idx, trial, stim, Oz) %>%
      pivot_wider(names_from = stim, values_from = Oz) %>%
      group_by(idx) %>%
      summarise(nt_change = mean(`Post-Stim` - `Pre-Stim`, na.rm = TRUE), .groups = "drop")

    d_result = effectsize::cohens_d(change_data$nt_change, mu = 0)

    tibble(
      mean_change = mean(change_data$nt_change),
      sd_change   = sd(change_data$nt_change),
      n           = nrow(change_data),
      cohen_d     = d_result$Cohens_d,
      ci_lower    = d_result$CI_low,
      ci_upper    = d_result$CI_high
    )
  }

  observed_ab = bind_rows(
    observed_d(rl.EST.Reward, "DA") %>% mutate(test = "a_DA_RL_task_specificity"),
    observed_d(ug.EST.Offer,  "SE") %>% mutate(test = "b_5HT_UG_task_specificity")
  )

  cat("\n===== Supplementary Table 34: observed effect sizes, DA x RL and 5-HT x UG =====\n")
  print(observed_ab)

  # Test (c)'s observed effect size
  ug.nt.cl = ug.EST.Offer %>%
    pivot_longer(cols = c("Oz","Rz","Pz","Mz"), names_to = "nt_metric", values_to = "nt_val") %>%
    group_by(idx, stim, nt, nt_metric) %>%
    summarise(mTrial = mean(nt_val), .groups = "drop") %>%
    filter(nt_metric == "Oz") %>%
    select(idx, stim, nt, mTrial) %>%
    pivot_wider(values_from = mTrial, names_from = c("stim","nt")) %>%
    mutate(deltaDA_UG = `Post-Stim_DA` - `Pre-Stim_DA`,
           deltaSE_UG = `Post-Stim_SE` - `Pre-Stim_SE`) %>%
    select(idx, deltaDA_UG, deltaSE_UG)

  data_m6_power = cl.df %>%
    filter(session == "month 6") %>%
    select(idx, HDRS) %>%
    merge(ug.nt.cl)

  model_DAxSE_check = lm(HDRS ~ deltaDA_UG * deltaSE_UG, data = data_m6_power)
  observed_c = tidy(model_DAxSE_check) %>%
    filter(term == "deltaDA_UG:deltaSE_UG") %>%
    mutate(test = "c_DAx5HT_HDRS_key_interaction")

  # Cross-check against the main-text model_DAxSE (estimate -0.142, P = 0.00282)
  confirmed_estimate_from_main_script = -0.142
  confirmed_pvalue_from_main_script   = 0.00282
  stopifnot(
    "Interaction estimate does not match main_analyses.R's model_DAxSE -- do not proceed until reconciled" =
      isTRUE(all.equal(round(observed_c$estimate, 3), round(confirmed_estimate_from_main_script, 3), tolerance = 1e-6)),
    "Interaction p-value does not match main_analyses.R's model_DAxSE -- do not proceed until reconciled" =
      isTRUE(all.equal(round(observed_c$p.value, 3), round(confirmed_pvalue_from_main_script, 3), tolerance = 1e-6))
  )

  # Partial R2 / f2 for the interaction term via nested-model comparison
  model_DAxSE_noint = lm(HDRS ~ deltaDA_UG + deltaSE_UG, data = data_m6_power)
  sse_reduced = sum(residuals(model_DAxSE_noint)^2)   # model WITHOUT the interaction
  sse_full    = sum(residuals(model_DAxSE_check)^2)   # model WITH the interaction
  observed_partial_R2 = (sse_reduced - sse_full) / sse_reduced
  observed_f2 = observed_partial_R2 / (1 - observed_partial_R2)

  # Re-verify the assumed residual df against the fitted model
  actual_v_test_c = df.residual(model_DAxSE_check)
  if (actual_v_test_c != assumed_v_test_c) {
    warning(sprintf(
      "MDE assumed v=%d residual df, but model_DAxSE_check actually has v=%d. ",
      assumed_v_test_c, actual_v_test_c),
      "Recomputing primary_c_mde at the correct v.")
    primary_c_mde = bind_rows(
      mde_interaction(u = 1, v = actual_v_test_c, alpha = alpha_bonferroni)  %>% mutate(correction = "Bonferroni-within-3 (matches main_analyses.R)"),
      mde_interaction(u = 1, v = actual_v_test_c, alpha = alpha_uncorrected) %>% mutate(correction = "uncorrected, for reference only")
    ) %>%
      mutate(test = "c_DAx5HT_HDRS_key_interaction")
    print(primary_c_mde)
  }

  cat("\n===== Supplementary Table 34: observed partial R2, DA x 5-HT -> Month 6 HDRS interaction =====\n")
  print(observed_c %>% mutate(observed_partial_R2 = observed_partial_R2, observed_f2 = observed_f2))

  # Batten, S. R. et al. Nat. Hum. Behav. 8, 718-728 (2024): published trial-level effects
  # (transcribed, not re-derived), converted to common d and r scales. Effect-magnitude
  # grounding only; not power-comparable to the subject-level MDEs above.
  r_from_d = function(d) d / sqrt(d^2 + 4)
  d_from_r = function(r) 2 * r / sqrt(1 - r^2)

  batten_grounding = tribble(
    ~batten_effect,                                                   ~native_beta, ~native_metric,  ~stat,           ~p,
    "Overall DA level, human vs. computer condition",                  0.85,        "d_equivalent",  "t(454) = 3.04", 0.002,
    "Relative DA tracks trial-by-trial value difference (RPE-like)",   0.21,        "approx_r",      "t(426) = 2.82", 0.005,
    "Relative 5-HT tracks current offer value",                       0.22,        "approx_r",      "t(426) = 3.15", 0.002
  ) %>%
    mutate(
      common_d_scale = ifelse(native_metric == "d_equivalent", native_beta, d_from_r(native_beta)),
      common_r_scale = ifelse(native_metric == "approx_r",     native_beta, r_from_d(native_beta))
    )

  cat("\n===== Supplementary Table 35 =====\n")
  print(batten_grounding)
})
