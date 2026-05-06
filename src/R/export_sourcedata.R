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
# 2026/03/19    Blair Shevlin                         wrote original code
#
# Generates figure source data CSVs for all figures in the manuscript:
#   Shevlin, Fu et al. - "Dopamine and serotonin transients predict depressive
#   symptom relief following deep brain stimulation of human subcallosal cingulate cortex"
#
# Each section below corresponds to one figure and writes a CSV to data/figures/.
# CSVs contain a `panel` column (labels which subplot each row belongs to) and a
# `level` column ("trial" or "subject") to distinguish granularity.
#
# Run this script ONCE before running any generate_figure*.R script.
# Output: data/figures/figure*_source_data.csv (one file per figure)
#
# These CSV files are the source data submitted with the manuscript per
# Nature journal transparency requirements.

library(tidyverse)
library(fs)
library(here)

# Create output directory
dir.create(here::here("data/figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("results/from_source_data"), recursive = TRUE, showWarnings = FALSE)

cat("Exporting figure source data...\n")


# ===== Figure 1 Panel E =====
cat(" Figure 1 Panel E...\n") 
{
    nt_dir   <- here::here("data/nt/processed")
    clin_dir <- here::here("data/clinical")

    load(file.path(nt_dir, "UG_RL_NT-Continuous_9-23-25.RData"))

  rl.EST.Reward <- rl.EST %>%
    filter(event == "Reward") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim", "Post-Stim"))) %>%
    select(!c(O, R, P, M, O10, Oz10, R10, Rz10, P10, Pz10, M10, Mz10))

  ug.EST.Offer <- ug.EST %>%
    filter(event == "Offer") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim", "Post-Stim"))) %>%
    select(!c(O, R, P, M, O10, Oz10, R10, Rz10, P10, Pz10, M10, Mz10))

  cl.df <- read.csv(file.path(clin_dir, "clinical-data_deid_07-10-25.csv"))

  ug.nt.cl <- ug.EST.Offer %>%
    pivot_longer(cols = c("Oz", "Rz", "Pz", "Mz"), names_to = "nt_metric", values_to = "nt_val") %>%
    group_by(idx, stim, nt, nt_metric) %>%
    summarise(mTrial = mean(nt_val)) %>%
    filter(nt_metric == "Oz") %>%
    select(idx, stim, nt, mTrial) %>%
    pivot_wider(values_from = mTrial, names_from = c("stim", "nt")) %>%
    mutate(deltaDA_UG = `Post-Stim_DA` - `Pre-Stim_DA`,
           deltaSE_UG = `Post-Stim_SE` - `Pre-Stim_SE`,
           deltaNE_UG = `Post-Stim_NE` - `Pre-Stim_NE`) %>%
    select(idx, deltaDA_UG, deltaSE_UG, deltaNE_UG)

  rl.nt.cl <- rl.EST.Reward %>%
    pivot_longer(cols = c("Oz", "Rz", "Pz", "Mz"), names_to = "nt_metric", values_to = "nt_val") %>%
    group_by(idx, stim, nt, nt_metric) %>%
    summarise(mTrial = mean(nt_val)) %>%
    filter(nt_metric == "Oz") %>%
    select(idx, stim, nt, mTrial) %>%
    pivot_wider(values_from = mTrial, names_from = c("stim", "nt")) %>%
    mutate(deltaDA_RL = `Post-Stim_DA` - `Pre-Stim_DA`,
           deltaSE_RL = `Post-Stim_SE` - `Pre-Stim_SE`,
           deltaNE_RL = `Post-Stim_NE` - `Pre-Stim_NE`) %>%
    select(idx, deltaDA_RL, deltaSE_RL, deltaNE_RL)

  nt.cl <- merge(rl.nt.cl, ug.nt.cl, all.y = TRUE)
  cl.Oz <- merge(cl.df, nt.cl)

  # HDRS trajectories data
  cl.Oz.fig <- cl.Oz %>%
    filter(!session %in% c("fmri")) %>%
    mutate(idx     = factor(idx),
           sess_fig = factor(session,
                             levels = c("pre stim", "post stim", "week 1", "month 1", "month 2",
                                        "month 3", "month 4", "month 5", "month 6"),
                             labels = c("Baseline", "DBS", "Week 1", "Month 1", "Month 2",
                                        "Month 3", "Month 4", "Month 5", "Month 6"))) %>%
    group_by(idx) %>%
    mutate(HDRS_CP      = (HDRS - HDRS[session == "pre stim"]) / HDRS[session == "pre stim"],
           remission    = ifelse(HDRS < 8, 1, 0),
           responder    = factor(ifelse(abs(HDRS_CP) > .5, 1, 0), levels = c(1, 0),
                                 labels = c("Responder", "Nonresponder")),
           M6_responder = responder[sess_fig == "Month 6"],
           M6_remission = remission[sess_fig == "Month 6"]) %>%
    ungroup() %>%
    as.data.frame()

  hdrs_baseline <- cl.Oz.fig %>%
    filter(sess_fig == "Baseline") %>%
    select(idx, HDRS_base = HDRS)


  # HDRS trajectories
  panel_E <- cl.Oz.fig %>%
    filter(sess_fig %in% c("Baseline", "Month 6")) %>%
    left_join(hdrs_baseline, by = "idx") %>%
    transmute(panel            = "E",
              level            = "subject",
              idx              = as.character(idx),
              session          = as.character(sess_fig),
              HDRS,
              HDRS_change      = HDRS - HDRS_base,
              outcome_category = as.character(M6_responder))

  source_data <- panel_E

  write.csv(source_data,
            here::here("data/figures/figure1_panel_E_source_data.csv"),
            row.names = FALSE)

  rm(rl.EST, ug.EST, rl.EST.Reward, ug.EST.Offer, cl.df,
     rl.nt.cl, ug.nt.cl, nt.cl, cl.Oz, cl.Oz.fig, hdrs_baseline,
     panel_E, source_data)
}

# ===== FIGURE 2 =====

cat("  Figure 2...\n")
{
  nt_dir <- here::here("data/nt/processed")

  load(file.path(nt_dir, "UG_RL_NT-Continuous_9-23-25.RData"))

  rl.EST.Reward <- rl.EST %>%
    filter(event == "Reward") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim", "Post-Stim"))) %>%
    select(!c(O, R, P, M, O10, Oz10, R10, Rz10, P10, Pz10, M10, Mz10))

  ug.EST.Offer <- ug.EST %>%
    filter(event == "Offer") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim", "Post-Stim"))) %>%
    select(!c(O, R, P, M, O10, Oz10, R10, Rz10, P10, Pz10, M10, Mz10))

  # Trial-level summaries
  rl.trial <- rl.EST.Reward %>%
    filter(nt %in% c("DA", "SE")) %>%
    pivot_longer(cols = c("Oz", "Rz", "Pz", "Mz"), names_to = "nt_metric", values_to = "nt_val") %>%
    mutate(block_type = factor(cond, levels = c("Mixed", "Negative", "Positive"))) %>%
    group_by(idx, stim, nt, nt_metric, trial, block_type, trial_within_block, prev_rew_raw) %>%
    summarise(mTrial = mean(nt_val),
              prev_rew = ifelse(prev_rew_raw == 0,
                                ifelse(block_type == "Negative" | block_type == "Mixed", -10, prev_rew_raw),
                                ifelse(prev_rew_raw == 1,
                                       ifelse(block_type == "Positive" | block_type == "Mixed", 10, 0),
                                       prev_rew_raw))) %>%
    ungroup()

  ug.trial <- ug.EST.Offer %>%
    filter(nt %in% c("DA", "SE")) %>%
    pivot_longer(cols = c("Oz", "Rz", "Pz", "Mz"), names_to = "nt_metric", values_to = "nt_val") %>%
    mutate(offer_bin = ifelse(offer < 4, "$1-3",
                              ifelse(offer < 7, "$4-6", "$7-9")),
           offer_bin = factor(offer_bin, levels = c("$1-3", "$4-6", "$7-9"))) %>%
    group_by(idx, stim, nt, nt_metric, trial, offer, offer_bin) %>%
    summarise(mTrial = mean(nt_val)) %>%
    ungroup() %>%
    mutate(offerz = scale(offer)[, 1])

  # Panel A: trial-level traces
  panel_A_RL <- rl.trial %>%
    filter(nt_metric == "Oz") %>%
    group_by(idx, stim, nt, trial) %>%
    summarise(Oz = mean(mTrial)) %>%
    mutate(nt = factor(nt, levels = c("DA", "SE"), labels = c("DA", "5-HT")),
           stim = factor(stim, levels = c("Pre-Stim", "Post-Stim"))) %>%
    transmute(panel = "A", level = "trial", task = "RL", nt, stim, idx,
              trial, condition = NA_character_, Oz)

  panel_A_UG <- ug.trial %>%
    filter(nt_metric == "Oz") %>%
    group_by(idx, stim, nt, trial) %>%
    summarise(Oz = mean(mTrial)) %>%
    mutate(nt = factor(nt, levels = c("DA", "SE"), labels = c("DA", "5-HT")),
           stim = factor(stim, levels = c("Pre-Stim", "Post-Stim"))) %>%
    transmute(panel = "A", level = "trial", task = "UG", nt, stim, idx,
              trial, condition = NA_character_, Oz)

  # Panel B: subject-level averages by stim
  rl.avg <- rl.trial %>%
    filter(nt_metric == "Oz") %>%
    group_by(idx, stim, nt) %>%
    summarise(Oz = mean(mTrial)) %>%
    mutate(nt = factor(nt, levels = c("DA", "SE"), labels = c("DA", "5-HT")),
           stim = factor(stim, levels = c("Pre-Stim", "Post-Stim"))) %>%
    transmute(panel = "B", level = "subject", task = "RL", nt, stim, idx,
              trial = NA_integer_, condition = NA_character_, Oz)

  ug.avg <- ug.trial %>%
    filter(nt_metric == "Oz") %>%
    group_by(idx, stim, nt) %>%
    summarise(Oz = mean(mTrial)) %>%
    mutate(nt = factor(nt, levels = c("DA", "SE"), labels = c("DA", "5-HT")),
           stim = factor(stim, levels = c("Pre-Stim", "Post-Stim"))) %>%
    transmute(panel = "B", level = "subject", task = "UG", nt, stim, idx,
              trial = NA_integer_, condition = NA_character_, Oz)

  # Panel C: subject-level by condition/offer bin
  rl.task <- rl.trial %>%
    filter(nt_metric == "Oz") %>%
    group_by(idx, stim, nt, block_type) %>%
    summarise(Oz = mean(mTrial)) %>%
    mutate(nt = factor(nt, levels = c("DA", "SE"), labels = c("DA", "5-HT")),
           stim = factor(stim, levels = c("Pre-Stim", "Post-Stim")),
           condition = as.character(block_type)) %>%
    transmute(panel = "C", level = "subject", task = "RL", nt, stim, idx,
              trial = NA_integer_, condition, Oz)

  ug.task <- ug.trial %>%
    filter(nt_metric == "Oz") %>%
    group_by(idx, stim, nt, offer_bin) %>%
    summarise(Oz = mean(mTrial)) %>%
    mutate(nt = factor(nt, levels = c("DA", "SE"), labels = c("DA", "5-HT")),
           stim = factor(stim, levels = c("Pre-Stim", "Post-Stim")),
           condition = as.character(offer_bin)) %>%
    transmute(panel = "C", level = "subject", task = "UG", nt, stim, idx,
              trial = NA_integer_, condition, Oz)

  source_data <- bind_rows(panel_A_RL, panel_A_UG,
                           rl.avg, ug.avg,
                           rl.task, ug.task) %>%
    mutate(across(where(is.factor), as.character))

  write.csv(source_data,
            here::here("data/figures/figure2_source_data.csv"),
            row.names = FALSE)

  rm(rl.EST, ug.EST, rl.EST.Reward, ug.EST.Offer,
     rl.trial, ug.trial, panel_A_RL, panel_A_UG,
     rl.avg, ug.avg, rl.task, ug.task, source_data)
}


# ===== FIGURE 3 =====

cat("  Figure 3...\n")
{
  nt_dir  <- here::here("data/nt/processed")
  beh_dir <- here::here("data/behavior")

  load(file.path(nt_dir, "UG_RL_NT-Continuous_9-23-25.RData"))

  rl.EST.Reward <- rl.EST %>%
    filter(event == "Reward") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim", "Post-Stim")))

  ug.EST.Offer <- ug.EST %>%
    filter(event == "Offer") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim", "Post-Stim")))

  rl.beh <- read.csv(file.path(beh_dir, "rl-data_deid_07-10-25.csv"))
  ug.beh <- read.csv(file.path(beh_dir, "ug-data_deid_07-10-25.csv"))

  ug.beh.means <- ug.beh %>%
    group_by(idx, sess) %>%
    summarise(mChoice = mean(rej == 0), mRT = mean(rt), mLogRT = mean(log(rt))) %>%
    mutate(idx = factor(idx), sess = recode(sess, "fMRI" = "Baseline")) %>%
    as.data.frame()

  ug.mood.means <- ug.beh %>%
    group_by(idx, sess) %>%
    filter(!is.na(rt_mood)) %>%
    summarise(mMood = mean(mood), mRT = mean(rt_mood), mLogRT = mean(log(rt))) %>%
    mutate(idx = factor(idx), sess = recode(sess, "fMRI" = "Baseline")) %>%
    as.data.frame()

  rl.beh.means <- rl.beh %>%
    group_by(idx, sess) %>%
    summarise(mChoice = mean(opt), mRT = mean(rt), mLogRT = mean(log(rt))) %>%
    mutate(idx = factor(idx), sess = recode(sess, "fMRI" = "Baseline")) %>%
    as.data.frame()

  add_dbs_rows <- function(df, value_cols) {
    unique_ids <- unique(df$idx)
    dbs_rows   <- data.frame(idx = unique_ids, sess = "DBS")
    for (col in value_cols) dbs_rows[[col]] <- NA
    bind_rows(df, dbs_rows) %>%
      mutate(sess = factor(sess,
                           levels = c("Baseline", "Pre-Stim", "DBS", "Post-Stim", "Week 1",
                                      "Month 1", "Month 2", "Month 3", "Month 4", "Month 5", "Month 6"))) %>%
      arrange(idx, sess)
  }

  rl.beh.means  <- add_dbs_rows(rl.beh.means,  c("mChoice", "mRT", "mLogRT"))
  ug.beh.means  <- add_dbs_rows(ug.beh.means,  c("mChoice", "mRT", "mLogRT"))
  ug.mood.means <- add_dbs_rows(ug.mood.means, c("mMood", "mRT", "mLogRT"))

  # Build top panel rows (behavioral timelines)
  panel_RL_choice <- rl.beh.means %>%
    transmute(panel = "RL_choice", level = "subject", task = "RL",
              idx = as.character(idx), sess = as.character(sess),
              measure = "mChoice", value = mChoice,
              se = NA_real_, n = NA_integer_,
              nt_change = NA_real_, beh_change = NA_real_, month = NA_character_)

  panel_RL_RT <- rl.beh.means %>%
    transmute(panel = "RL_RT", level = "subject", task = "RL",
              idx = as.character(idx), sess = as.character(sess),
              measure = "mLogRT", value = mLogRT,
              se = NA_real_, n = NA_integer_,
              nt_change = NA_real_, beh_change = NA_real_, month = NA_character_)

  panel_UG_choice <- ug.beh.means %>%
    transmute(panel = "UG_choice", level = "subject", task = "UG",
              idx = as.character(idx), sess = as.character(sess),
              measure = "mChoice", value = mChoice,
              se = NA_real_, n = NA_integer_,
              nt_change = NA_real_, beh_change = NA_real_, month = NA_character_)

  panel_UG_RT <- ug.beh.means %>%
    transmute(panel = "UG_RT", level = "subject", task = "UG",
              idx = as.character(idx), sess = as.character(sess),
              measure = "mLogRT", value = mLogRT,
              se = NA_real_, n = NA_integer_,
              nt_change = NA_real_, beh_change = NA_real_, month = NA_character_)

  panel_UG_mood <- ug.mood.means %>%
    transmute(panel = "UG_mood", level = "subject", task = "UG",
              idx = as.character(idx), sess = as.character(sess),
              measure = "mMood", value = mMood,
              se = NA_real_, n = NA_integer_,
              nt_change = NA_real_, beh_change = NA_real_, month = NA_character_)

  # NT change scores for correlation panels
  rl.nt <- rl.EST.Reward %>%
    filter(nt %in% c("DA", "SE")) %>%
    group_by(idx, nt, stim) %>%
    summarise(est = mean(Oz)) %>%
    group_by(idx, nt) %>%
    summarise(dff = est[stim == "Post-Stim"] - est[stim == "Pre-Stim"]) %>%
    group_by(idx) %>%
    reframe(DA = dff[nt == "DA"], SE = dff[nt == "SE"])

  ug.nt <- ug.EST.Offer %>%
    filter(nt %in% c("DA", "SE")) %>%
    group_by(idx, nt, stim) %>%
    summarise(est = mean(Oz)) %>%
    group_by(idx, nt) %>%
    summarise(dff = est[stim == "Post-Stim"] - est[stim == "Pre-Stim"]) %>%
    group_by(idx) %>%
    reframe(DA = dff[nt == "DA"], SE = dff[nt == "SE"])

  rl.beh.change <- rl.beh.means %>%
    mutate(mRT = mLogRT) %>%
    select(!mLogRT) %>%
    pivot_wider(names_from = sess, values_from = c(mChoice, mRT)) %>%
    mutate(mRT_PostStim_Change = `mRT_Post-Stim` - mRT_Baseline,
           mRT_W1_Change  = `mRT_Week 1`  - mRT_Baseline,
           mRT_M1_Change  = `mRT_Month 1` - mRT_Baseline,
           mRT_M2_Change  = `mRT_Month 2` - mRT_Baseline,
           mRT_M3_Change  = `mRT_Month 3` - mRT_Baseline,
           mRT_M4_Change  = `mRT_Month 4` - mRT_Baseline,
           mRT_M5_Change  = `mRT_Month 5` - mRT_Baseline,
           mRT_M6_Change  = `mRT_Month 6` - mRT_Baseline) %>%
    select(idx, mRT_PostStim_Change, mRT_W1_Change, mRT_M1_Change, mRT_M2_Change,
           mRT_M3_Change, mRT_M4_Change, mRT_M5_Change, mRT_M6_Change)

  rl.beh.nt <- merge(rl.beh.change, rl.nt)

  rl.beh.nt.long <- rl.beh.nt %>%
    pivot_longer(cols = c(mRT_PostStim_Change, mRT_W1_Change, mRT_M1_Change, mRT_M2_Change,
                          mRT_M3_Change, mRT_M4_Change, mRT_M5_Change, mRT_M6_Change),
                 names_to = "change", values_to = "beh_change") %>%
    mutate(month     = str_extract(change, "(?<=_)[^_]+(?=_Change)"),
           nt_change = SE)

  ug.mood.change <- ug.mood.means %>%
    mutate(mRT = mLogRT) %>%
    select(!mLogRT) %>%
    pivot_wider(names_from = sess, values_from = c(mMood, mRT)) %>%
    mutate(mMood_PostStim_Change = `mMood_Post-Stim` - mMood_Baseline,
           mMood_W1_Change  = `mMood_Week 1`  - mMood_Baseline,
           mMood_M1_Change  = `mMood_Month 1` - mMood_Baseline,
           mMood_M2_Change  = `mMood_Month 2` - mMood_Baseline,
           mMood_M3_Change  = `mMood_Month 3` - mMood_Baseline,
           mMood_M4_Change  = `mMood_Month 4` - mMood_Baseline,
           mMood_M5_Change  = `mMood_Month 5` - mMood_Baseline,
           mMood_M6_Change  = `mMood_Month 6` - mMood_Baseline) %>%
    select(idx, mMood_PostStim_Change, mMood_W1_Change, mMood_M1_Change, mMood_M2_Change,
           mMood_M3_Change, mMood_M4_Change, mMood_M5_Change, mMood_M6_Change)

  ug.mood.nt <- merge(ug.mood.change, ug.nt)

  ug.mood.nt.long <- ug.mood.nt %>%
    pivot_longer(cols = c(mMood_PostStim_Change, mMood_W1_Change, mMood_M1_Change, mMood_M2_Change,
                          mMood_M3_Change, mMood_M4_Change, mMood_M5_Change, mMood_M6_Change),
                 names_to = "change", values_to = "beh_change") %>%
    mutate(month     = str_extract(change, "(?<=_)[^_]+(?=_Change)"),
           nt_change = DA)

  panel_corr_5HT_RL <- rl.beh.nt.long %>%
    filter(month %in% c("M1", "M3", "M6")) %>%
    transmute(panel = "corr_5HT_RL", level = "subject", task = "RL",
              idx = as.character(idx), sess = NA_character_,
              measure = "rt", value = NA_real_,
              se = NA_real_, n = NA_integer_,
              nt_change, beh_change, month)

  panel_corr_DA_UG <- ug.mood.nt.long %>%
    filter(month %in% c("M1", "M3", "M6")) %>%
    transmute(panel = "corr_DA_UG", level = "subject", task = "UG",
              idx = as.character(idx), sess = NA_character_,
              measure = "mood", value = NA_real_,
              se = NA_real_, n = NA_integer_,
              nt_change, beh_change, month)

  source_data <- bind_rows(panel_RL_choice, panel_RL_RT,
                           panel_UG_choice, panel_UG_RT, panel_UG_mood,
                           panel_corr_5HT_RL, panel_corr_DA_UG)

  write.csv(source_data,
            here::here("data/figures/figure3_source_data.csv"),
            row.names = FALSE)

  rm(rl.EST, ug.EST, rl.EST.Reward, ug.EST.Offer,
     rl.beh, ug.beh, rl.beh.means, ug.beh.means, ug.mood.means,
     rl.nt, ug.nt, rl.beh.change, rl.beh.nt, rl.beh.nt.long,
     ug.mood.change, ug.mood.nt, ug.mood.nt.long,
     panel_RL_choice, panel_RL_RT, panel_UG_choice, panel_UG_RT, panel_UG_mood,
     panel_corr_5HT_RL, panel_corr_DA_UG, source_data)
}


# ===== FIGURE 4 =====

cat("  Figure 4...\n")
{
  nt_dir   <- here::here("data/nt/processed")
  clin_dir <- here::here("data/clinical")

  load(file.path(nt_dir, "UG_RL_NT-Continuous_9-23-25.RData"))

  rl.EST.Reward <- rl.EST %>%
    filter(event == "Reward") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim", "Post-Stim"))) %>%
    select(!c(O, R, P, M, O10, Oz10, R10, Rz10, P10, Pz10, M10, Mz10))

  ug.EST.Offer <- ug.EST %>%
    filter(event == "Offer") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim", "Post-Stim"))) %>%
    select(!c(O, R, P, M, O10, Oz10, R10, Rz10, P10, Pz10, M10, Mz10))

  cl.df <- read.csv(file.path(clin_dir, "clinical-data_deid_07-10-25.csv"))

  ug.nt.cl <- ug.EST.Offer %>%
    pivot_longer(cols = c("Oz", "Rz", "Pz", "Mz"), names_to = "nt_metric", values_to = "nt_val") %>%
    group_by(idx, stim, nt, nt_metric) %>%
    summarise(mTrial = mean(nt_val)) %>%
    filter(nt_metric == "Oz") %>%
    select(idx, stim, nt, mTrial) %>%
    pivot_wider(values_from = mTrial, names_from = c("stim", "nt")) %>%
    mutate(deltaDA_UG = `Post-Stim_DA` - `Pre-Stim_DA`,
           deltaSE_UG = `Post-Stim_SE` - `Pre-Stim_SE`,
           deltaNE_UG = `Post-Stim_NE` - `Pre-Stim_NE`) %>%
    select(idx, deltaDA_UG, deltaSE_UG, deltaNE_UG)

  rl.nt.cl <- rl.EST.Reward %>%
    pivot_longer(cols = c("Oz", "Rz", "Pz", "Mz"), names_to = "nt_metric", values_to = "nt_val") %>%
    group_by(idx, stim, nt, nt_metric) %>%
    summarise(mTrial = mean(nt_val)) %>%
    filter(nt_metric == "Oz") %>%
    select(idx, stim, nt, mTrial) %>%
    pivot_wider(values_from = mTrial, names_from = c("stim", "nt")) %>%
    mutate(deltaDA_RL = `Post-Stim_DA` - `Pre-Stim_DA`,
           deltaSE_RL = `Post-Stim_SE` - `Pre-Stim_SE`,
           deltaNE_RL = `Post-Stim_NE` - `Pre-Stim_NE`) %>%
    select(idx, deltaDA_RL, deltaSE_RL, deltaNE_RL)

  nt.cl <- merge(rl.nt.cl, ug.nt.cl, all.y = TRUE)
  cl.Oz <- merge(cl.df, nt.cl)

  # HDRS trajectories data
  cl.Oz.fig <- cl.Oz %>%
    filter(!session %in% c("fmri")) %>%
    mutate(idx     = factor(idx),
           sess_fig = factor(session,
                             levels = c("pre stim", "post stim", "week 1", "month 1", "month 2",
                                        "month 3", "month 4", "month 5", "month 6"),
                             labels = c("Baseline", "DBS", "Week 1", "Month 1", "Month 2",
                                        "Month 3", "Month 4", "Month 5", "Month 6"))) %>%
    group_by(idx) %>%
    mutate(HDRS_CP      = (HDRS - HDRS[session == "pre stim"]) / HDRS[session == "pre stim"],
           remission    = ifelse(HDRS < 8, 1, 0),
           responder    = factor(ifelse(abs(HDRS_CP) > .5, 1, 0), levels = c(1, 0),
                                 labels = c("Responder", "Nonresponder")),
           M6_responder = responder[sess_fig == "Month 6"],
           M6_remission = remission[sess_fig == "Month 6"]) %>%
    ungroup() %>%
    as.data.frame()

  hdrs_baseline <- cl.Oz.fig %>%
    filter(sess_fig == "Baseline") %>%
    select(idx, HDRS_base = HDRS)

  # NT change metrics
  data_m6 <- cl.Oz %>%
    filter(session == "month 6") %>%
    select(idx, HDRS, deltaDA_UG, deltaSE_UG, deltaDA_RL, deltaSE_RL) %>%
    left_join(cl.Oz %>% filter(session == "pre stim") %>% select(idx, HDRS_base = HDRS), by = "idx") %>%
    mutate(deltaPerHDRS     = (HDRS_base - HDRS) / HDRS_base,
           change_magnitude = HDRS - HDRS_base,
           synergy_score    = case_when(
             deltaDA_UG > 0 & deltaSE_UG > 0 ~ deltaDA_UG * deltaSE_UG,
             deltaDA_UG < 0 & deltaSE_UG < 0 ~ abs(deltaDA_UG * deltaSE_UG),
             TRUE                             ~ -abs(deltaDA_UG * deltaSE_UG)
           ),
           response_category = case_when(
             HDRS <= 8          ~ "Remission",
             deltaPerHDRS > 0.5 ~ "Responder",
             TRUE               ~ "Non-Responder"
           ),
           response_category = factor(response_category,
                                      levels = c("Non-Responder", "Responder", "Remission")),
           change_pattern = case_when(
             deltaDA_UG > 0 & deltaSE_UG > 0 ~ "Both Increase",
             deltaDA_UG < 0 & deltaSE_UG < 0 ~ "Both Decrease",
             deltaDA_UG > 0 & deltaSE_UG < 0 ~ "DA\u2191/5-HT\u2193",
             TRUE                             ~ "DA\u2193/5-HT\u2191"
           ))

  # Panel A: DA vs 5-HT scatter
  panel_A <- data_m6 %>%
    transmute(panel            = "A",
              level            = "subject",
              idx              = as.character(idx),
              session          = NA_character_,
              HDRS,
              HDRS_change      = change_magnitude,
              outcome_category = as.character(response_category),
              deltaDA          = deltaDA_UG,
              deltaSE          = deltaSE_UG,
              task             = "UG",
              nt_pattern       = change_pattern,
              synergy_score    = NA_real_,
              HDRS_m6          = NA_real_,
              HDRS_change_m6   = NA_real_)

  # Panel B: synergy vs HDRS
  panel_B <- data_m6 %>%
    transmute(panel            = "B",
              level            = "subject",
              idx              = as.character(idx),
              session          = NA_character_,
              HDRS,
              HDRS_change      = HDRS - HDRS_base,
              outcome_category = as.character(response_category),
              deltaDA          = NA_real_,
              deltaSE          = NA_real_,
              task             = NA_character_,
              nt_pattern       = change_pattern,
              synergy_score,
              HDRS_m6          = HDRS,
              HDRS_change_m6   = change_magnitude)

  source_data <- bind_rows(panel_A, panel_B)

  write.csv(source_data,
            here::here("data/figures/figure4_source_data.csv"),
            row.names = FALSE)

  rm(rl.EST, ug.EST, rl.EST.Reward, ug.EST.Offer, cl.df,
     rl.nt.cl, ug.nt.cl, nt.cl, cl.Oz, cl.Oz.fig, hdrs_baseline, data_m6,
     panel_A, panel_B, source_data)
}


# ===== FIGURE EXTENDED 1 =====
cat(" Extended Figure 1...\n") 
{
  clin_dir <- here::here("data/clinical")

  cl.df <- read.csv(file.path(clin_dir, "clinical-data_deid_07-10-25.csv"))

  # HDRS trajectories data
  cl.fig <- cl.df %>%
    filter(!session %in% c("fmri")) %>%
    mutate(idx     = factor(idx),
           sess_fig = factor(session,
                             levels = c("pre stim", "post stim", "week 1", "month 1", "month 2",
                                        "month 3", "month 4", "month 5", "month 6"),
                             labels = c("Baseline", "DBS", "Week 1", "Month 1", "Month 2",
                                        "Month 3", "Month 4", "Month 5", "Month 6"))) %>%
    group_by(idx) %>%
    mutate(HDRS_CP      = (HDRS - HDRS[session == "pre stim"]) / HDRS[session == "pre stim"],
           remission    = ifelse(HDRS < 8, 1, 0),
           responder    = factor(ifelse(abs(HDRS_CP) > .5, 1, 0), levels = c(1, 0),
                                 labels = c("Responder", "Nonresponder")),
           M6_responder = responder[sess_fig == "Month 6"],
           M6_remission = remission[sess_fig == "Month 6"]) %>%
    ungroup() %>%
    as.data.frame()

  hdrs_baseline <- cl.fig %>%
    filter(sess_fig == "Baseline") %>%
    select(idx, HDRS_base = HDRS)

  # HDRS trajectories
  source_data <- cl.fig %>%
    filter(sess_fig %in% c("Baseline", "Week 1", "Month 1", "Month 2", "Month 3", "Month 4", "Month 5", "Month 6")) %>%
    left_join(hdrs_baseline, by = "idx") %>%
    transmute(level            = "subject",
              idx              = as.character(idx),
              session          = as.character(sess_fig),
              HDRS,
              HDRS_change      = HDRS - HDRS_base,
              outcome_category = as.character(responder))

  write.csv(source_data,
            here::here("data/figures/figureExtended1_source_data.csv"),
            row.names = FALSE)

  rm(cl.df,cl.fig, hdrs_baseline, source_data)
}

# ===== Extended FIGURE 2 =====

cat(" Extended Figure 2...\n")
{
  library(caTools)

  nt_dir <- here::here("data/nt/processed")

  load(file.path(nt_dir, "UG_RL_NT-Raw_10-08-25.RData"))

  rl.EST.Reward <- rl_proc %>%
    filter(event == "Reward")

  sm.rl <- rl.EST.Reward %>%
    filter(!is.na(nM)) %>%
    mutate(stim   = factor(stim, levels = c("Pre-Stim", "Post-Stim")),
           sample = ifelse(sample == 51, 0, (sample - 51) * 100),
           sm     = runmean(nM, 3, align = "right")) %>%
    group_by(idx, stim, nt, trial) %>%
    mutate(szm = runmean(nMz, 3, align = "right")) %>%
    filter(sample %in% c(-400:800))

  fig.data <- sm.rl %>%
    filter(event        == "Reward",
           stim         == "Post-Stim",
           nt           == "SE",
           trial        == 28,
           idx          == 5,
           !is.na(prev_rew))

  panel_C <- fig.data %>%
    transmute(
      panel         = "C",
      level         = "trial",
      time_sample   = sample,
      nM_raw        = nM,
      nM_smoothed   = sm,
      z_score       = szm,
      in_auc_window = (sample >= 0 & sample <= 500)
    )

  source_data <- panel_C

  write.csv(source_data,
            here::here("data/figures/figureExtended2_source_data.csv"),
            row.names = FALSE)

  rm(rl_proc, rl.EST.Reward, sm.rl, fig.data, panel_C, source_data)
}

# ===== Extended FIGURE 2 =====
cat("  Figure Extended 3...\n")
{
  nt_dir <- here::here("data/nt/processed")

  load(file.path(nt_dir, "UG_RL_NT-Continuous_9-23-25.RData"))

  rl.EST.Reward <- rl.EST %>%
    filter(event == "Reward") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim", "Post-Stim")))

  ug.EST.Offer <- ug.EST %>%
    filter(event == "Offer") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim", "Post-Stim")))

  ug.trial <- ug.EST.Offer %>%
    pivot_longer(cols = c("Oz", "Rz", "Pz", "Mz"), names_to = "nt_metric", values_to = "nt_val") %>%
    mutate(offer_bin = ifelse(offer < 4, "$1-3", ifelse(offer < 7, "$4-6", "$7-9")),
           offer_bin = factor(offer_bin, levels = c("$1-3", "$4-6", "$7-9"))) %>%
    group_by(idx, stim, nt, nt_metric, trial, offer, offer_bin) %>%
    summarise(mTrial = mean(nt_val)) %>%
    ungroup()

  rl.trial <- rl.EST.Reward %>%
    pivot_longer(cols = c("Oz", "Rz", "Pz", "Mz"), names_to = "nt_metric", values_to = "nt_val") %>%
    mutate(block_type = factor(cond, levels = c("Mixed", "Negative", "Positive"))) %>%
    group_by(idx, stim, nt, nt_metric, trial, block_type, trial_within_block, outcome) %>%
    summarise(mTrial = mean(nt_val)) %>%
    ungroup()

  # Panel A: RL subject-level averages (NE only)
  panel_A <- rl.trial %>%
    filter(nt_metric == "Oz", nt == "NE") %>%
    group_by(idx, stim, nt) %>%
    summarise(Oz = mean(mTrial)) %>%
    mutate(panel = "A", level = "subject", task = "RL", trial = NA_real_)

  # Panel B: UG subject-level averages (NE only)
  panel_B <- ug.trial %>%
    filter(nt_metric == "Oz", nt == "NE") %>%
    group_by(idx, stim, nt) %>%
    summarise(Oz = mean(mTrial)) %>%
    mutate(panel = "B", level = "subject", task = "UG", trial = NA_real_)

  # Panel C: RL trial-level (NE only)
  panel_C <- rl.trial %>%
    filter(nt_metric == "Oz", nt == "NE") %>%
    group_by(idx, stim, nt, trial) %>%
    summarise(Oz = mean(mTrial)) %>%
    mutate(panel = "C", level = "trial", task = "RL")

  # Panel D: UG trial-level (NE only)
  panel_D <- ug.trial %>%
    filter(nt_metric == "Oz", nt == "NE") %>%
    group_by(idx, stim, nt, trial) %>%
    summarise(Oz = mean(mTrial)) %>%
    mutate(panel = "D", level = "trial", task = "UG")

  source_data <- bind_rows(panel_A, panel_B, panel_C, panel_D) %>%
    select(panel, level, task, nt, stim, idx, trial, Oz)

  write.csv(source_data,
            here::here("data/figures/figureExtended3_source_data.csv"),
            row.names = FALSE)

  rm(rl.EST, ug.EST, rl.EST.Reward, ug.EST.Offer,
     rl.trial, ug.trial, panel_A, panel_B, panel_C, panel_D, source_data)
}


# ===== FIGURE EXTENDED 4 =====

cat("  Figure Extended 4...\n")
{
  clin_dir <- here::here("data/clinical")

  cl.df <- read.csv(file.path(clin_dir, "clinical-data_deid_07-10-25.csv"))

  cl.O.fig <- cl.df %>%
    filter(!session %in% c("fmri")) %>%
    mutate(id       = factor(idx),
           sess_fig = factor(session,
                             levels = c("pre stim", "post stim", "week 1", "month 1",
                                        "month 2", "month 3", "month 4", "month 5", "month 6"),
                             labels = c("Baseline", "DBS", "Week 1", "Month 1",
                                        "Month 2", "Month 3", "Month 4", "Month 5", "Month 6"))) %>%
    group_by(id) %>%
    mutate(HDRS_CP     = (HDRS - HDRS[session == "pre stim"]) / HDRS[session == "pre stim"],
           remission   = ifelse(HDRS < 8, 1, 0),
           responder   = factor(ifelse(abs(HDRS_CP) > .5, 1, 0), levels = c(1, 0),
                                labels = c("Responder", "Nonresponder")),
           M6_responder = responder[sess_fig == "Month 6"],
           M6_remission = remission[sess_fig == "Month 6"],
           cohort       = ifelse(id %in% c(1:5), "RC+S", "PC")) %>%
    ungroup() %>%
    filter(!sess_fig %in% c("Month 2", "Month 4", "Month 5")) %>%
    as.data.frame()

  source_data <- cl.O.fig %>%
    mutate(panel           = "A",
           level           = "subject",
           pct_change_HDRS = HDRS_CP) %>%
    select(panel, level, idx, session, HDRS, pct_change_HDRS, cohort, responder)

  write.csv(source_data,
            here::here("data/figures/figureExtended4_source_data.csv"),
            row.names = FALSE)

  rm(cl.df, cl.O.fig, source_data)
}


# ===== FIGURE SUPPLEMENT 1 =====

cat("  Figure Supplement 1...\n")
{
  se_fn <- function(x) sd(x) / sqrt(length(x))

  nt_dir   <- here::here("data/nt/processed")
  clin_dir <- here::here("data/clinical")

  load(file.path(nt_dir, "UG_RL_NT-Continuous_Raw-Model-Output_11-25-25.RData"))

  rl.EST.Reward <- rl.EST %>%
    filter(event == "Reward") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim", "Post-Stim"))) %>%
    select(!c(R, P, M, O10, Oz10, R10, Rz10, P10, Pz10, M10, Mz10))

  ug.EST.Offer <- ug.EST %>%
    filter(event == "Offer") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim", "Post-Stim"))) %>%
    select(!c(R, P, M, O10, Oz10, R10, Rz10, P10, Pz10, M10, Mz10))

  rl.subj_summary <- rl.EST.Reward %>%
    filter(nt %in% c("DA", "SE")) %>%
    pivot_longer(cols = c("Oz", "O"), names_to = "nt_metric", values_to = "nt_val") %>%
    filter(nt_metric == "O") %>%
    group_by(idx, stim, nt) %>%
    summarise(Oz = mean(nt_val, na.rm = TRUE),
              se = se_fn(nt_val),
              .groups = "drop") %>%
    mutate(nt   = factor(nt, levels = c("DA", "SE"), labels = c("DA", "5-HT")),
           stim = factor(stim, levels = c("Pre-Stim", "Post-Stim")))

  ug.subj_summary <- ug.EST.Offer %>%
    filter(nt %in% c("DA", "SE")) %>%
    pivot_longer(cols = c("Oz", "O"), names_to = "nt_metric", values_to = "nt_val") %>%
    filter(nt_metric == "O") %>%
    group_by(idx, stim, nt) %>%
    summarise(Oz = mean(nt_val, na.rm = TRUE),
              se = se_fn(nt_val),
              .groups = "drop") %>%
    mutate(nt   = factor(nt, levels = c("DA", "SE"), labels = c("DA", "5-HT")),
           stim = factor(stim, levels = c("Pre-Stim", "Post-Stim")))

  cl.df <- read.csv(file.path(clin_dir, "clinical-data_deid_07-10-25.csv")) %>%
    filter(session == "pre stim") %>%
    select(idx, HDRS)

  rl.cl.df <- merge(rl.subj_summary, cl.df, by = "idx")
  ug.cl.df <- merge(ug.subj_summary, cl.df, by = "idx")

  panel_A <- rl.subj_summary %>% mutate(panel = "A", level = "subject", task = "RL", HDRS = NA_real_)
  panel_B <- ug.subj_summary %>% mutate(panel = "B", level = "subject", task = "UG", HDRS = NA_real_)
  panel_C <- rl.cl.df %>% filter(stim == "Pre-Stim") %>% mutate(panel = "C", level = "subject", task = "RL")
  panel_D <- ug.cl.df %>% filter(stim == "Pre-Stim") %>% mutate(panel = "D", level = "subject", task = "UG")

  source_data <- bind_rows(panel_A, panel_B, panel_C, panel_D) %>%
    select(panel, level, task, nt, stim, idx, Oz, se, HDRS)

  write.csv(source_data,
            here::here("data/figures/figureSupplement1_source_data.csv"),
            row.names = FALSE)

  rm(rl.EST, ug.EST, rl.EST.Reward, ug.EST.Offer,
     rl.subj_summary, ug.subj_summary, cl.df, rl.cl.df, ug.cl.df,
     panel_A, panel_B, panel_C, panel_D, source_data)
}


# ===== FIGURE SUPPLEMENT 2 =====

cat("  Figure Supplement 2...\n")
{
  library(lmerTest)
  library(emmeans)

  nt_dir <- here::here("data/nt/processed")

  load(file.path(nt_dir, "UG_RL_NT-Continuous_9-23-25.RData"))

  rl.EST.Reward <- rl.EST %>%
    filter(event == "Reward") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim", "Post-Stim")))

  ug.EST.Offer <- ug.EST %>%
    filter(event == "Offer") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim", "Post-Stim")))

  rl.premerge <- rl.EST.Reward %>%
    mutate(trial_z = scale(trial)[, 1]) %>%
    select(idx, stim, nt, Oz, trial_z) %>%
    mutate(task = "RL")

  ug.premerge <- ug.EST.Offer %>%
    mutate(trial_z = scale(trial)[, 1]) %>%
    select(idx, stim, nt, Oz, trial_z) %>%
    mutate(task = "UG")

  task.merge <- rbind(rl.premerge, ug.premerge) %>% mutate(task = factor(task))

  sess.lm.DA <- lmer(Oz ~ stim * task * trial_z + (1 + stim + task | idx),
                     data    = task.merge[task.merge$nt == "DA", ],
                     REML    = FALSE,
                     control = lmerControl(optimizer = "nmkbw"))
  sess.lm.SE <- lmer(Oz ~ stim * task * trial_z + (1 + stim + task | idx),
                     data    = task.merge[task.merge$nt == "SE", ],
                     REML    = FALSE,
                     control = lmerControl(optimizer = "nmkbw"))
  sess.lm.NE <- lmer(Oz ~ stim * task * trial_z + (1 + stim + task | idx),
                     data    = task.merge[task.merge$nt == "NE", ],
                     REML    = FALSE,
                     control = lmerControl(optimizer = "nmkbw"))

  da_slopes_df <- emtrends(sess.lm.DA, ~ stim * task, var = "trial_z") %>% as.data.frame()
  se_slopes_df <- emtrends(sess.lm.SE, ~ stim * task, var = "trial_z") %>% as.data.frame()
  ne_slopes_df <- emtrends(sess.lm.NE, ~ stim * task, var = "trial_z") %>% as.data.frame()

  build_slope_rows <- function(df, nt_label) {
    df %>%
      transmute(nt          = nt_label,
                task        = as.character(task),
                stim        = as.character(stim),
                slope       = trial_z.trend,
                slope_lower = lower.CL,
                slope_upper = upper.CL)
  }

  slope_data <- bind_rows(
    build_slope_rows(da_slopes_df, "DA"),
    build_slope_rows(se_slopes_df, "5-HT"),
    build_slope_rows(ne_slopes_df, "NE")
  )

  trial_q <- as.numeric(quantile(task.merge$trial_z, probs = c(.25, .75)))

  early_late_data <- task.merge %>%
    group_by(idx, nt, task, stim) %>%
    mutate(period = ifelse(trial_z <= min(trial_q), "Early",
                           ifelse(trial_z >= max(trial_q), "Late", "Middle"))) %>%
    filter(period != "Middle") %>%
    group_by(idx, nt, task, stim, period) %>%
    summarise(Oz = mean(Oz), .groups = "drop") %>%
    mutate(nt = recode(nt, "SE" = "5-HT"),
           nt = factor(nt, levels = c("DA", "5-HT", "NE")))

  panel_A <- slope_data %>%
    mutate(panel        = "A",
           level        = "group",
           idx          = NA_character_,
           trial_period = NA_character_,
           Oz           = NA_real_) %>%
    select(panel, level, task, nt, stim, idx, trial_period, Oz, slope, slope_lower, slope_upper)

  panel_B <- early_late_data %>%
    filter(task == "UG") %>%
    mutate(panel        = "B",
           level        = "subject",
           trial_period = period,
           slope        = NA_real_,
           slope_lower  = NA_real_,
           slope_upper  = NA_real_) %>%
    select(panel, level, task, nt, stim, idx, trial_period, Oz, slope, slope_lower, slope_upper)

  panel_C <- early_late_data %>%
    filter(task == "RL") %>%
    mutate(panel        = "C",
           level        = "subject",
           trial_period = period,
           slope        = NA_real_,
           slope_lower  = NA_real_,
           slope_upper  = NA_real_) %>%
    select(panel, level, task, nt, stim, idx, trial_period, Oz, slope, slope_lower, slope_upper)

  source_data <- bind_rows(panel_A, panel_B, panel_C)

  write.csv(source_data,
            here::here("data/figures/figureSupplement2_source_data.csv"),
            row.names = FALSE)

  rm(rl.EST, ug.EST, rl.EST.Reward, ug.EST.Offer,
     rl.premerge, ug.premerge, task.merge,
     sess.lm.DA, sess.lm.SE, sess.lm.NE,
     da_slopes_df, se_slopes_df, ne_slopes_df, slope_data,
     early_late_data, panel_A, panel_B, panel_C, source_data)
}


# ===== FIGURE SUPPLEMENT 3 =====

cat("  Figure Supplement 3...\n")
{
  beh_dir <- here::here("data/behavior")

  rl.beh <- read.csv(file.path(beh_dir, "rl-data_deid_07-10-25.csv"))
  ug.beh <- read.csv(file.path(beh_dir, "ug-data_deid_07-10-25.csv"))

  sessions <- c("Baseline", "Pre-Stim", "Post-Stim", "DBS", "Month 1", "Month 3", "Month 6")

  add_dbs_rows <- function(df, value_cols) {
    unique_ids <- unique(df$idx)
    dbs_rows   <- data.frame(idx = unique_ids, sess = "DBS")
    for (col in value_cols) dbs_rows[[col]] <- NA
    bind_rows(df, dbs_rows) %>%
      mutate(sess = factor(sess,
                           levels = c("Baseline", "Pre-Stim", "DBS", "Post-Stim", "Week 1",
                                      "Month 1", "Month 2", "Month 3", "Month 4", "Month 5", "Month 6"))) %>%
      arrange(idx, sess)
  }

  ug.beh.means <- ug.beh %>%
    group_by(idx, sess) %>%
    summarise(mChoice = mean(rej == 0), mRT = mean(rt), mLogRT = mean(log(rt))) %>%
    mutate(idx = factor(idx), sess = recode(sess, "fMRI" = "Baseline")) %>%
    as.data.frame()

  ug.mood.means <- ug.beh %>%
    group_by(idx, sess) %>%
    filter(!is.na(rt_mood)) %>%
    summarise(mMood = mean(mood), mRT = mean(rt_mood), mLogRT = mean(log(rt))) %>%
    mutate(idx = factor(idx), sess = recode(sess, "fMRI" = "Baseline")) %>%
    as.data.frame()

  rl.beh.means <- rl.beh %>%
    group_by(idx, sess) %>%
    summarise(mChoice = mean(opt), mRT = mean(rt), mLogRT = mean(log(rt))) %>%
    mutate(idx = factor(idx), sess = recode(sess, "fMRI" = "Baseline")) %>%
    as.data.frame()

  rl.beh.means  <- add_dbs_rows(rl.beh.means,  c("mChoice", "mRT", "mLogRT"))
  ug.beh.means  <- add_dbs_rows(ug.beh.means,  c("mChoice", "mRT", "mLogRT"))
  ug.mood.means <- add_dbs_rows(ug.mood.means, c("mMood", "mRT", "mLogRT"))

  rl.beh.means.plot  <- rl.beh.means  %>% filter(sess %in% sessions) %>% mutate(sess = factor(sess, levels = sessions))
  ug.beh.means.plot  <- ug.beh.means  %>% filter(sess %in% sessions) %>% mutate(sess = factor(sess, levels = sessions))
  ug.mood.means.plot <- ug.mood.means %>% filter(sess %in% sessions) %>% mutate(sess = factor(sess, levels = sessions))

  panel_RL_choice <- rl.beh.means.plot %>%
    transmute(panel = "RL_choice", level = "subject", task = "RL",
              idx, sess, mChoice, mLogRT = NA, mMood = NA)

  panel_RL_RT <- rl.beh.means.plot %>%
    transmute(panel = "RL_RT", level = "subject", task = "RL",
              idx, sess, mChoice = NA, mLogRT, mMood = NA)

  panel_UG_choice <- ug.beh.means.plot %>%
    transmute(panel = "UG_choice", level = "subject", task = "UG",
              idx, sess, mChoice, mLogRT = NA, mMood = NA)

  panel_UG_RT <- ug.beh.means.plot %>%
    transmute(panel = "UG_RT", level = "subject", task = "UG",
              idx, sess, mChoice = NA, mLogRT, mMood = NA)

  panel_UG_mood <- ug.mood.means.plot %>%
    transmute(panel = "UG_mood", level = "subject", task = "UG",
              idx, sess, mChoice = NA, mLogRT = NA, mMood)

  source_data <- bind_rows(panel_RL_choice, panel_RL_RT,
                           panel_UG_choice, panel_UG_RT, panel_UG_mood)

  write.csv(source_data,
            here::here("data/figures/figureSupplement3_source_data.csv"),
            row.names = FALSE)

  rm(rl.beh, ug.beh, rl.beh.means, ug.beh.means, ug.mood.means,
     rl.beh.means.plot, ug.beh.means.plot, ug.mood.means.plot,
     panel_RL_choice, panel_RL_RT, panel_UG_choice, panel_UG_RT, panel_UG_mood, source_data)
}


# ===== FIGURE SUPPLEMENT 4 =====

cat("  Figure Supplement 4...\n")
{
  beh_dir  <- here::here("data/behavior")
  clin_dir <- here::here("data/clinical")

  cl.df  <- read.csv(file.path(clin_dir, "clinical-data_deid_07-10-25.csv"))
  rl.beh <- read.csv(file.path(beh_dir, "rl-data_deid_07-10-25.csv"))
  ug.beh <- read.csv(file.path(beh_dir, "ug-data_deid_07-10-25.csv"))

  cl.change <- cl.df %>%
    filter(!session %in% c("fmri", "post stim")) %>%
    mutate(sess = recode(session,
                         "pre stim" = "Baseline",
                         "week 1"   = "Week 1",
                         "month 1"  = "Month 1",
                         "month 2"  = "Month 2",
                         "month 3"  = "Month 3",
                         "month 4"  = "Month 4",
                         "month 5"  = "Month 5",
                         "month 6"  = "Month 6")) %>%
    select(!c(MADRS, session)) %>%
    group_by(idx) %>%
    mutate(HDRS_d = HDRS - HDRS[sess == "Baseline"])

  rl.beh.change <- rl.beh %>%
    filter(rt < 10) %>%
    group_by(idx, sess) %>%
    summarise(mChoice = mean(opt), mRT = mean(rt), mLogRT = mean(log(rt))) %>%
    mutate(idx = factor(idx), sess = recode(sess, "fMRI" = "Baseline"))

  ug.beh.change <- ug.beh %>%
    filter(rt < 10) %>%
    group_by(idx, sess) %>%
    summarise(mChoice = mean(rej == 0), mRT = mean(rt),
              mLogRT = mean(log(rt)), mMood = mean(mood, na.rm = TRUE)) %>%
    mutate(idx = factor(idx), sess = recode(sess, "fMRI" = "Baseline")) %>%
    as.data.frame()

  rl.beh.cl   <- merge(rl.beh.change, cl.change)
  ug.beh.cl   <- merge(ug.beh.change, cl.change)

  rl.beh.cl.l <- rl.beh.cl %>%
    pivot_longer(cols = c("mChoice", "mRT", "mLogRT"), names_to = "beh")
  ug.beh.cl.l <- ug.beh.cl %>%
    pivot_longer(cols = c("mChoice", "mRT", "mLogRT", "mMood"), names_to = "beh")

  target_sessions <- c("Month 1", "Month 3", "Month 6")

  panel_A <- rl.beh.cl.l %>%
    filter(sess %in% target_sessions, beh == "mChoice") %>%
    transmute(panel = "A", level = "subject", task = "RL", measure = "mChoice",
              idx, session = sess, beh_change = value, HDRS_change = HDRS_d)

  panel_B <- rl.beh.cl.l %>%
    filter(sess %in% target_sessions, beh == "mLogRT") %>%
    transmute(panel = "B", level = "subject", task = "RL", measure = "mLogRT",
              idx, session = sess, beh_change = value, HDRS_change = HDRS_d)

  panel_C <- ug.beh.cl.l %>%
    filter(sess %in% target_sessions, beh == "mMood") %>%
    transmute(panel = "C", level = "subject", task = "UG", measure = "mMood",
              idx, session = sess, beh_change = value, HDRS_change = HDRS_d)

  panel_D <- ug.beh.cl.l %>%
    filter(sess %in% target_sessions, beh == "mChoice") %>%
    transmute(panel = "D", level = "subject", task = "UG", measure = "mChoice",
              idx, session = sess, beh_change = value, HDRS_change = HDRS_d)

  panel_E <- ug.beh.cl.l %>%
    filter(sess %in% target_sessions, beh == "mLogRT") %>%
    transmute(panel = "E", level = "subject", task = "UG", measure = "mLogRT",
              idx, session = sess, beh_change = value, HDRS_change = HDRS_d)

  source_data <- bind_rows(panel_A, panel_B, panel_C, panel_D, panel_E)

  write.csv(source_data,
            here::here("data/figures/figureSupplement4_source_data.csv"),
            row.names = FALSE)

  rm(cl.df, rl.beh, ug.beh, cl.change, rl.beh.change, ug.beh.change,
     rl.beh.cl, ug.beh.cl, rl.beh.cl.l, ug.beh.cl.l,
     panel_A, panel_B, panel_C, panel_D, panel_E, source_data)
}


# ===== FIGURE SUPPLEMENT 5 =====

cat("  Figure Supplement 5...\n")
{
  nt_dir  <- here::here("data/nt/processed")
  beh_dir <- here::here("data/behavior")

  rl.beh <- read.csv(file.path(beh_dir, "rl-data_deid_07-10-25.csv"))
  ug.beh <- read.csv(file.path(beh_dir, "ug-data_deid_07-10-25.csv"))

  load(file.path(nt_dir, "UG_RL_NT-Continuous_9-23-25.RData"))

  # tasks.EST.alt uses cross-task z-scoring
  tasks.EST.alt <- tasks.EST.alt %>%
    group_by(task, nt, event, idx) %>%
    mutate(trial_z = scale(trial)[, 1]) %>%
    ungroup()

  rl.EST.Shared <- tasks.EST.alt %>%
    filter(task == "RL", event == "Shared")

  ug.EST.Shared <- tasks.EST.alt %>%
    filter(task == "UG", event == "Shared")

  rl.task.OR <- rl.beh %>%
    filter(sess %in% c("Pre-Stim", "Post-Stim")) %>%
    mutate(stim = factor(sess, levels = c("Pre-Stim", "Post-Stim"))) %>%
    select(!c(session, sess, score))

  ug.task.OR <- ug.beh %>%
    filter(sess %in% c("Pre-Stim", "Post-Stim")) %>%
    mutate(stim = factor(sess, levels = c("Pre-Stim", "Post-Stim"))) %>%
    select(!c(session, sess))

  rl.EST.Shared <- merge(rl.EST.Shared, rl.task.OR, by = c("idx", "stim", "trial"))
  ug.EST.Shared <- merge(ug.EST.Shared, ug.task.OR)

  ug.trial <- ug.EST.Shared %>%
    filter(nt %in% c("DA", "SE")) %>%
    pivot_longer(cols = c("Oz", "Rz", "Pz", "Mz"), names_to = "nt_metric", values_to = "nt_val") %>%
    mutate(offer_bin = ifelse(offer < 4, "$1-3", ifelse(offer < 7, "$4-6", "$7-9")),
           offer_bin = factor(offer_bin, levels = c("$1-3", "$4-6", "$7-9"))) %>%
    group_by(idx, stim, nt, nt_metric, trial, offer, offer_bin) %>%
    summarise(mTrial = mean(nt_val)) %>%
    ungroup()

  rl.trial <- rl.EST.Shared %>%
    filter(nt %in% c("DA", "SE")) %>%
    pivot_longer(cols = c("Oz", "Rz", "Pz", "Mz"), names_to = "nt_metric", values_to = "nt_val") %>%
    mutate(block_type = factor(cond, levels = c("Mixed", "Negative", "Positive"))) %>%
    group_by(idx, stim, nt, nt_metric, trial, block_type, trial_within_block) %>%
    summarise(mTrial = mean(nt_val)) %>%
    ungroup()

  recode_nt <- function(x) factor(recode(x, "SE" = "5-HT"), levels = c("DA", "5-HT"))

  panel_A <- rl.trial %>%
    filter(nt_metric == "Oz") %>%
    group_by(idx, stim, nt, trial) %>%
    summarise(Oz = mean(mTrial)) %>%
    ungroup() %>%
    mutate(nt = recode_nt(nt), panel = "A", level = "trial", task = "RL", condition = NA_character_)

  panel_B <- ug.trial %>%
    filter(nt_metric == "Oz") %>%
    group_by(idx, stim, nt, trial) %>%
    summarise(Oz = mean(mTrial)) %>%
    ungroup() %>%
    mutate(nt = recode_nt(nt), panel = "B", level = "trial", task = "UG", condition = NA_character_)

  panel_C <- rl.trial %>%
    filter(nt_metric == "Oz") %>%
    group_by(idx, stim, nt) %>%
    summarise(Oz = mean(mTrial)) %>%
    ungroup() %>%
    mutate(nt = recode_nt(nt), panel = "C", level = "subject", task = "RL",
           trial = NA_integer_, condition = NA_character_)

  panel_D <- ug.trial %>%
    filter(nt_metric == "Oz") %>%
    group_by(idx, stim, nt) %>%
    summarise(Oz = mean(mTrial)) %>%
    ungroup() %>%
    mutate(nt = recode_nt(nt), panel = "D", level = "subject", task = "UG",
           trial = NA_integer_, condition = NA_character_)

  panel_E <- rl.trial %>%
    filter(nt_metric == "Oz") %>%
    group_by(idx, stim, nt, block_type) %>%
    summarise(Oz = mean(mTrial)) %>%
    ungroup() %>%
    mutate(nt        = recode_nt(nt),
           block_type = factor(block_type,
                               levels = c("Positive", "Mixed", "Negative"),
                               labels = c("Reward", "Mixed", "Punishment")),
           panel     = "E", level = "subject", task = "RL",
           trial     = NA_integer_, condition = as.character(block_type))

  panel_F <- ug.trial %>%
    filter(nt_metric == "Oz") %>%
    group_by(idx, stim, nt, offer_bin) %>%
    summarise(Oz = mean(mTrial)) %>%
    ungroup() %>%
    mutate(nt        = recode_nt(nt),
           panel     = "F", level = "subject", task = "UG",
           trial     = NA_integer_, condition = as.character(offer_bin))

  source_data <- bind_rows(
    panel_A, panel_B, panel_C, panel_D, panel_E, panel_F
  ) %>%
    select(panel, level, task, nt, stim, idx, trial, condition, Oz) %>%
    mutate(across(where(is.factor), as.character))

  write.csv(source_data,
            here::here("data/figures/figureSupplement5_source_data.csv"),
            row.names = FALSE)

  rm(rl.beh, ug.beh, rl.EST, ug.EST, tasks.EST.alt,
     rl.EST.Shared, ug.EST.Shared, rl.task.OR, ug.task.OR,
     rl.trial, ug.trial, panel_A, panel_B, panel_C, panel_D, panel_E, panel_F, source_data)
}


# ===== FIGURE SUPPLEMENT 6 =====

cat("  Figure Supplement 6...\n")
{
  nt_dir <- here::here("data/nt/processed")

  load(file.path(nt_dir, "UG_RL_NT-Continuous_9-23-25.RData"))

  rl.EST.Reward <- rl.EST %>%
    filter(event == "Reward") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim", "Post-Stim"))) %>%
    select(!c(O, R, P, M, O10, Oz10, R10, Rz10, P10, Pz10, M10, Mz10))

  ug.EST.Offer <- ug.EST %>%
    filter(event == "Offer") %>%
    mutate(stim = factor(stim, levels = c("Pre-Stim", "Post-Stim"))) %>%
    select(!c(O, R, P, M, O10, Oz10, R10, Rz10, P10, Pz10, M10, Mz10))

  ug.trial <- ug.EST.Offer %>%
    filter(nt %in% c("DA", "SE")) %>%
    pivot_longer(cols = c("Oz", "Rz", "Pz", "Mz"), names_to = "nt_metric", values_to = "nt_val") %>%
    mutate(offer_bin = ifelse(offer < 4, "$1-3", ifelse(offer < 7, "$4-6", "$7-9")),
           offer_bin = factor(offer_bin, levels = c("$1-3", "$4-6", "$7-9"))) %>%
    group_by(idx, stim, nt, nt_metric, trial, offer, offer_bin) %>%
    summarise(mTrial = mean(nt_val)) %>%
    ungroup()

  rl.trial <- rl.EST.Reward %>%
    filter(nt %in% c("DA", "SE")) %>%
    pivot_longer(cols = c("Oz", "Rz", "Pz", "Mz"), names_to = "nt_metric", values_to = "nt_val") %>%
    mutate(block_type = factor(cond, levels = c("Mixed", "Negative", "Positive"))) %>%
    group_by(idx, stim, nt, nt_metric, trial, block_type, trial_within_block, prev_rew_raw) %>%
    summarise(mTrial = mean(nt_val)) %>%
    ungroup()

  panel_A <- rl.trial %>%
    filter(nt_metric == "Oz") %>%
    group_by(idx, stim, nt, trial) %>%
    summarise(Oz = mean(mTrial)) %>%
    ungroup() %>%
    mutate(panel = "A", level = "trial", task = "RL")

  panel_B <- ug.trial %>%
    filter(nt_metric == "Oz") %>%
    group_by(idx, stim, nt, trial) %>%
    summarise(Oz = mean(mTrial)) %>%
    ungroup() %>%
    mutate(panel = "B", level = "trial", task = "UG")

  source_data <- bind_rows(
    panel_A %>% select(panel, level, task, nt, stim, idx, trial, Oz),
    panel_B %>% select(panel, level, task, nt, stim, idx, trial, Oz)
  )

  write.csv(source_data,
            here::here("data/figures/figureSupplement6_source_data.csv"),
            row.names = FALSE)

  rm(rl.EST, ug.EST, rl.EST.Reward, ug.EST.Offer,
     rl.trial, ug.trial, panel_A, panel_B, source_data)
}


cat("Done! All figure source data CSVs written to data/figures/\n")
