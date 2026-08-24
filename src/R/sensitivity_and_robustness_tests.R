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
# 2025/05/05    Blair Shevlin                           wrote original code
# 2025/07/24    Blair Shevlin                           updated to use new NT data
# 2025/09/16    Blair Shevlin                           updated to use revised NT data
# 2025/09/17    Blair Shevlin                           robustness tests for clinical predictions
# 2025/10/1     Blair Shevlin                           statistical tests comparing post-stim OR to baseline
# 2025/11/14    Blair Shevlin                           adding effect sizes
# 2025/11/18    Blair Shevlin                           added task x sess analysis
# 2026/05/05    Blair Shevlin                           added HDRS change analysis
# 2026/08/15    Blair Shevlin                           added robustness checks and multiple comparison corrections
# 2026/08/15    Blair Shevlin                           added bootstrap CIs for Fig 3f/g correlations (R1 sensitivity analysis)
# 2026/08/16    Blair Shevlin                           added raw-O (unnormalized) check of the DAxSE->HDRS interaction, per Xiaosi Gu's question re: whether subject 701's leverage is a per-subject Oz normalization artifact
# 2026/08/16    Claude (agent)                          corrected stale "unplanned/opportunistic" framing of test (c) in the Bonferroni-within-3 comment block to reflect Xiaosi Gu's Aug 15 correction (pre-specified Activity-3 hypothesis); see quick_facts.md and narrative_guidance.md
# 2026/08/16    Claude (agent)                          resolved the rl.lm.DA.2 vs. sess.lm.DA.2 open question (test a/b object mapping) via independent Python re-derivation from figure2_source_data.csv; see comment above the Bonferroni-within-3 block
# 2026/08/16    Claude (agent)                          added multi-method correction table (Bonferroni/Sidak/Holm/Hochberg/Hommel/BH/BY) per Xiaosi Gu's request to check whether Bonferroni is overly conservative for the primary family; see comment block before build_correction_table()
# 2026/08/17    Claude (agent)                          added two literal permutation-test variants for test (c) per Blair/PI's explicit request ("shuffled deltaDA x deltaSE, M6 HDRS scores"), distinct from the existing Freedman-Lane test above -- Option A: full-label permutation refitting the full interaction model; Option B: permutation test of the raw deltaDA_UG*deltaSE_UG product against HDRS (Pearson + exact Spearman). Both independently cross-checked in Python (obs interaction=-0.1419 matches model_DAxSE exactly; Option A perm p~=0.012; Option B Pearson perm p~=0.013, Spearman perm p~=0.045-0.049). See src/R/permutation_responder_robustness_plan.md sec 1.2-1.3 for full reasoning and the multiplicity/framing caveat repeated in-line below -- NEITHER variant is pre-approved for the rebuttal.
# 2026/08/18    Claude (agent)                          consolidated sensitivity/robustness analyses (bootstrap CIs, LOO, influence diagnostics, permutation tests, multiplicity correction) into a single script per Blair's request; MDE/Batten-effect-size content intentionally excluded (lives in power_mde_analysis.R)
# 2026/08/20    Claude (agent)                          added Fig 4 panel B (synergy_score / product_DAxSE -> Month-6 HDRS, raw Pearson correlation) bootstrap CI, LOO, and Cook's D/DFBETAS diagnostics as Table 31, extending the Table 28/30 R1 diagnostic battery to the one Fig 4 quantity it hadn't yet been applied to; permutation p-value intentionally left as an open discrepancy (see comment below) pending a fresh R run

# ==========================================================================
# WHAT THIS SCRIPT IS AND IS NOT
# ==========================================================================
# This script consolidates the manuscript's/rebuttal's robustness and
# sensitivity-check analyses -- bootstrap CIs, leave-one-out (LOO) influence
# diagnostics, permutation tests, and multiplicity correction -- into one
# place, adapted from the working file src/R/statistical_tests_robust.R
# (2,409 lines as of this writing). Nothing here is new statistics: every
# function and formula below is copied or lightly trimmed from that file: no
# numbers are invented, adjusted, or rounded by hand anywhere in this script.
# Every printed value is the direct output of running the code against the
# actual data objects.
#
# NAMING NOTE: Blair referred to the MDE/power-justification script verbally
# as "generate_power_mde_report.R." No file by that name exists in the repo.
# The actual file with that content is src/R/power_mde_analysis.R. That is
# the file being treated as the exclusion boundary below -- i.e., MDE /
# minimum-detectable-effect analysis and the Batten et al. (2024) external
# effect-size comparison are OUT OF SCOPE for this script and live there
# instead.
#
# EXPLICITLY OUT OF SCOPE for this script (all confirmed present in
# statistical_tests_robust.R but intentionally not reproduced here):
#   - MDE / Batten et al. (2024) effect-size comparison -- lives in
#     power_mde_analysis.R (see naming note above).
#   - HC3 heteroskedasticity-consistent SEs -- considered and explicitly
#     rejected in the source file's own comments (badly anti-conservative at
#     this n); out of scope here too.
#   - Huber/M-estimation robust regression (MASS::rlm) -- diagnostic-only in
#     the source file, not an inferential result; excluded.
#   - Bayesian regularized regression -- never run in the source file
#     (would need explicit sign-off + the bayesian-workflow skill guardrails
#     before it could be); excluded.
#   - The raw-O (per-subject-unnormalized) sensitivity check -- a separate
#     side investigation answering Xiaosi Gu's normalization-artifact
#     question, not part of the core robustness battery for the manuscript;
#     excluded.
#   - The Freedman-Lane residual-permutation test (permutation_test_interaction())
#     and its LOO-subset variant -- these exist in statistical_tests_robust.R
#     as an alternative, assumption-light permutation method and are
#     extensively discussed in that file's inline comments, but they are
#     NOT what the published rebuttal text actually cites. Kept out of this
#     script's execution path; see the Month-6 permutation section below for
#     the full explanation of why Option A (full-label permutation) is used
#     instead.
#   - Option B permutation (raw deltaDA_UG*deltaSE_UG product-term
#     correlation, Pearson + exact Spearman) -- not used in the manuscript;
#     excluded.
#   - Any joint multiplicity-correction comparison BETWEEN the Month 1 and
#     Month 6 DAx5-HT->HDRS findings -- Blair has ruled (Aug 2026) that this
#     comparison (already computed separately, as internal decision support,
#     in m1_hdrs_analysis/M1_M6_robustness_and_correction_report.md) must
#     never appear in manuscript, supplementary, or rebuttal-facing content.
#     Month 1 (Table 29) and Month 6 (Table 30) are therefore reported below
#     as two independent battery blocks; nothing here compares or jointly
#     corrects across them.
#
# TABLE NUMBERING (per manuscript_edit_placement_PLAN.md's current plan):
#   Supplementary Table 27 = primary/exploratory multiplicity correction
#   Supplementary Table 28 = Figure 3f/g bootstrap CIs + influence diagnostics
#   Supplementary Table 29 = Month 1  DAx5-HT->HDRS bootstrap + permutation only
#   Supplementary Table 30 = Month 6  DAx5-HT->HDRS full robustness battery
#
# R's source() does NOT auto-print top-level results the way an interactive
# console does (this bit the source file once already -- see the note near
# model_DAxSE_influence below). Every result meant to be read is wrapped in
# print() explicitly throughout this script.

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
library(broom.mixed)
library(caTools)
library("ggnewscale")
library(rstatix)
library(scales)
library(ggrepel)
library(effectsize)
library(emmeans)
library(flextable)
library(officer)

# Paths
dir = path(here())
nt_dir = dir / "data" / "nt" / "processed"
beh_dir = dir / "data" / "behavior"
clin_dir = dir / "data" / "clinical"
res_dir = dir / "results"

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
cl.df = read.csv(clin_dir / "clinical-data_deid_07-10-25.csv")

# Load behavioral data
rl.beh = read.csv(file = beh_dir / "rl-data_deid_07-10-25.csv")

ug.beh = read.csv(file = beh_dir / "ug-data_deid_07-10-25.csv")

##########################################################
# DATA DERIVATION (trimmed replication of the upstream    #
# pipeline in statistical_tests_robust.R)                 #
##########################################################
#
# This section is NOT a full copy of statistical_tests_robust.R's "STIM ON
# NT" / "NT ~ BEHAVIOR" / "NT ~ CLINICAL" sections -- those also build many
# objects (block_type/prev_rew/rew factor codings, paired t-tests, ANOVAs,
# task-by-stim interaction models, model-selection comparisons via anova(),
# ug.beh.nt.long, etc.) that feed analyses out of scope for this consolidated
# robustness script. What follows reproduces, faithfully and in the same
# order as the source file, only what is actually required to construct:
#   (a) rl.lm.DA.2 / ug.lm.SE.2  -- feed test (a)/(b) in the Table 27 battery
#   (b) rl.beh.nt.long / ug.mood.nt.long -- feed the Table 28 Fig 3f/g battery
#   (c) cl.Oz.lme -> data_m1 / data_m6, model_UG_m1, model_DAxSE -- feed the
#       Table 29 / Table 30 batteries
#
# Judgment call: trial_z (= scale(trial)[,1]) is the only column consumed by
# rl.lm.DA.2/ug.lm.SE.2 from the source file's larger rl.EST.pr/ug.EST.pr
# mutate() blocks; trial_z is mathematically independent of the other
# columns computed there (block_type, prev_rew, rew, etc., which feed
# reward/prediction-error models entirely outside this script's scope), so
# it is reproduced exactly without carrying the unrelated columns along.
# Likewise, the rl.lm.DA.1 / ug.lm.SE.1 "correlated random intercept"
# specifications and their anova() model-comparison step (which justified
# moving to the "0 + stim_num || idx" specification used below) are
# documented in statistical_tests_robust.R lines ~296-369 and are not
# re-run here, since only the already-justified final ".2" specification is
# actually consumed downstream.

## -- STIM ON NT: single-task DA (RL) / 5-HT (UG) models (test a / test b) --

rl.EST.pr = rl.EST.Reward %>%
  mutate(trial_z = scale(trial)[,1],
         stim_num = ifelse(stim == "Post-Stim", 1, -1))

ug.EST.pr = ug.EST.Offer %>%
  mutate(trial_z = scale(trial)[,1],
         stim_num = ifelse(stim == "Post-Stim", 1, -1))

rl.premerge = rl.EST.pr %>%
  select(idx, stim, stim_num, nt, Oz, trial_z) %>%
  mutate(task = "RL")
ug.premerge = ug.EST.pr %>%
  select(idx, stim, stim_num, nt, Oz, trial_z) %>%
  mutate(task = "UG")

task.merge = rbind(rl.premerge, ug.premerge) %>%
  mutate(task = factor(task),
         task_num = ifelse(task == "UG", 1, -1))

# DA increase, task specificity in RL -- feeds test (a) below
rl.lm.DA.2 = lmer(Oz ~ stim * trial_z +
                     (1 + stim|idx),
                   data = task.merge[task.merge$nt == "DA" & task.merge$task == "RL",],
                   REML = F,
                   control = lmerControl(optimizer = "bobyqa"))

# 5-HT increase, task specificity in UG -- feeds test (b) below
ug.lm.SE.2 = lmer(Oz ~ stim * trial_z +
                     (1 + stim|idx),
                   data = task.merge[task.merge$nt == "SE" & task.merge$task == "UG",],
                   REML = F,
                   control = lmerControl(optimizer = "bobyqa"))

## -- NT ~ BEHAVIOR: session-level means, NT change, and the long tables --
##    feeding the Figure 3f/g bootstrap CI + influence diagnostics battery

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
            mRT = mean(rt_mood),
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
           mRT_M3_Change,mRT_M4_Change,mRT_M5_Change,mRT_M6_Change,mMood_OR_Change,mRT_OR_Change))

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

## -- NT ~ CLINICAL: cl.Oz.lme -> data_m1 / data_m6, model_UG_m1, model_DAxSE --

ug.nt.cl =
  ug.EST.Offer %>%
  pivot_longer(cols = c("Oz","Rz","Pz","Mz","O"),
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
  pivot_longer(cols = c("Oz","Rz","Pz","Mz","O"),
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
  mutate(baseline_HDRS = HDRS[session == "fmri"],
         baseline_MADRS = MADRS[session == "fmri"],
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

data_m6 = cl.Oz.lme %>%
  filter(session == "month 6") %>%
  select(idx,HDRS,baseline_HDRS,deltaPerHDRS,
         deltaDA_UG,deltaSE_UG,deltaNE_UG,
         deltaDA_RL,deltaSE_RL,deltaNE_RL)

# Month 1 DA x 5-HT -> HDRS interaction model (used as-is for the Table 29
# bootstrap/permutation battery below -- NOT redefined there)
model_UG_m1 <- lm(HDRS ~ deltaDA_UG * deltaSE_UG, data = data_m1)
print(summary(model_UG_m1))

# Month 6 DA x 5-HT -> HDRS interaction model -- this is test (c), the
# manuscript's "key interaction," and feeds both Table 27 and Table 30 below.
model_DAxSE = lm(HDRS ~ deltaDA_UG * deltaSE_UG, data = data_m6)
print(summary(model_DAxSE))

##########################################################
# TABLE 27 -- ROBUSTNESS: BONFERRONI-WITHIN-3 FOR THE     #
# "KEY INTERACTION" (R1 point 1's specific, quoted         #
# requirement), PLUS THE MULTI-METHOD CORRECTION TABLE     #
##########################################################

# R1: "the key interaction should be confirmed to remain significant even after
# Bonferroni correction within this family of three." The "key interaction" is the
# DA x 5-HT -> Month 6 HDRS interaction (test c below).
#
# CORRECTED 2026/08/15, per Xiaosi Gu (see quick_facts.md and narrative_guidance.md):
# test (c) is NOT an unplanned/opportunistic finding. Per the funded LEAP grant
# proposal, it is a pre-specified Activity-3 hypothesis (using DA/5-HT + at-home
# measures to predict longitudinal symptom change) -- on the same pre-specified
# footing as tests (a)/(b), which are pre-specified Activity-2 hypotheses (acute
# stimulation effect on task-based DA/5-HT). An earlier version of this comment (and
# of robustness_analysis_plan.md section 1, since corrected) described test (c) as an
# unplanned/opportunistic discovery -- that framing is superseded and must not be used
# in the manuscript or rebuttal.
#
# The team's adopted multiplicity-correction position (quick_facts.md) is to correct
# WITHIN-Activity: an Activity-2 family {test a, test b} and an Activity-3 family
# {test c, plausibly + the Fig 3 Month 1/3/6 correlations -- confirm full inventory},
# NOT R1's literal single 3-test grouping spanning both Activities. The
# Bonferroni-within-3 computation below is retained anyway as a direct, literal answer
# to R1's own specific wording ("this family of three") -- a supplementary robustness
# check showing the key interaction survives even under R1's proposed (not the team's
# adopted) correction scheme. Present it in the rebuttal as exactly that -- one
# additional check, not the primary correction argument, which should instead explain
# and use the within-Activity split. tests (a)/(b) are reported uncorrected here per
# R1's own framing, alongside for completeness -- not because they are any less
# pre-specified than test (c).
#
# RESOLVED 2026/08/16 (Claude, independent verification): candidate (i), the
# task-specific single-task models rl.lm.DA.2 / ug.lm.SE.2, are the objects that
# produced the manuscript's already-reported "DA increase in reversal learning:
# beta=6.41, t=2.37, P=0.045" and "5-HT increase in ultimatum game: beta=5.62, t=2.57,
# P=0.028". If tidy(rl.lm.DA.2) / tidy(ug.lm.SE.2) is pulled directly for the
# Bonferroni table below, the printed "estimate" column will read ~half the
# manuscript's published beta (stim_num is coded -1/+1 here vs. 0/1 in Methods) -- this
# does NOT change the t-statistic, p-value, or significance conclusion (both are
# invariant to linear rescaling of the predictor), but multiply the estimate by 2
# before quoting it alongside the manuscript's beta so the two don't appear to
# disagree.

# --- Test (a): DA increase, task specificity in RL ---
test_a_single_task = tidy(rl.lm.DA.2) %>%
  filter(term == "stim") %>%
  mutate(test = "a_DA_RL_single_task_model")

# --- Test (b): 5-HT increase, task specificity in UG ---
test_b_single_task = tidy(ug.lm.SE.2) %>%
  filter(term == "stim") %>%
  mutate(test = "b_5HT_UG_single_task_model")

# --- Test (c): DA x 5-HT interaction -> Month 6 HDRS ("the key interaction") ---
test_c = tidy(model_DAxSE) %>%
  filter(term == "deltaDA_UG:deltaSE_UG") %>%
  mutate(test = "c_DAx5HT_HDRS_key_interaction")

# Bonferroni comparison, alpha = 0.05 / 3 = 0.0167 -- applies to test (c) only, per
# R1's actual wording ("the key interaction," singular). Tests (a)/(b) are reported
# uncorrected; they are not expected or required to survive this stricter threshold.
alpha_family = 0.05
alpha_bonferroni = alpha_family / 3

primary_family_table = bind_rows(
  test_a_single_task,
  test_b_single_task,
  test_c
) %>%
  mutate(
    is_key_interaction = test == "c_DAx5HT_HDRS_key_interaction",
    alpha_applied = ifelse(is_key_interaction, alpha_bonferroni, alpha_family),
    survives = p.value < alpha_applied
  ) %>%
  select(test, any_of("term"), estimate, p.value, is_key_interaction, alpha_applied, survives)

print(primary_family_table %>% as.data.frame())

# Per R1's own framing: a non-key test in this family failing to survive
# alpha=0.0167 (R1's own letter cites "DA increase in reversal learning: p=0.045" as
# exactly this case) is an expected, acceptable outcome -- report it plainly, do not
# treat it as a problem to explain away. Only the row(s) with is_key_interaction ==
# TRUE need to clear alpha_bonferroni for R1's specific requirement to be satisfied.

##########################################################
# ROBUSTNESS: MULTIPLE CORRECTION METHODS, NOT JUST      #
# BONFERRONI (Xiaosi Gu, Slack, Aug 2026: worried         #
# Bonferroni is too conservative for a 3-test family --   #
# wants alternative correction methods tried and reported #
# alongside it, not substituted for it silently)          #
##########################################################

# Xiaosi's concern is a reasonable one: Bonferroni controls the family-wise error rate
# (FWER) via the crudest possible bound (union bound over m tests) and assumes nothing
# about dependence between tests -- which makes it valid but frequently over-conservative,
# especially at m=3 where a single non-key test's failure can look more damning than it
# is. This section reports several standard alternatives ALONGSIDE Bonferroni, not as a
# replacement for it -- per this project's rigor guardrails, this is reported regardless
# of which method is most flattering, and the choice of which number(s) to lead with in
# the rebuttal is a decision for Blair/Xiaosi/co-authors, not something this script
# should make silently by only running the method that "works."
#
# All methods below via base R's p.adjust() (aligns exactly with Python's
# statsmodels.stats.multitest.multipletests when method names are matched -- verified
# independently, see below) applied to the SAME three uncorrected p-values already in
# primary_family_table: test (a) p=0.045, test (b) p=0.028, test (c) p=0.0028.
#
# What each method corrects for (ordered roughly most- to least-conservative here):
#  - bonferroni:   Controls FWER (P(>=1 false positive across the family)). Rejects if
#                  p < alpha/m. No assumption about dependence between tests -- valid
#                  under ANY dependence structure, which is exactly why it's conservative:
#                  it doesn't use any information about how the tests relate.
#  - sidak:        Controls FWER, ~identical to Bonferroni at small m/small p, but derived
#                  assuming the tests are independent (rejects if p < 1-(1-alpha)^(1/m)).
#                  Slightly less conservative than Bonferroni; the independence assumption
#                  is an approximation here since tests (a)/(b)/(c) share overlapping
#                  neurotransmitter/subject data.
#  - holm:         Controls FWER, uniformly at least as powerful as Bonferroni (rejects
#                  everything Bonferroni does, sometimes more). Step-DOWN procedure: sort
#                  p-values ascending, compare the smallest to alpha/m, the next to
#                  alpha/(m-1), etc., stopping at the first failure to reject. No
#                  dependence assumption -- like Bonferroni, valid under any dependence.
#  - hochberg:     Controls FWER under independence or positive dependence (satisfied
#                  here: all three tests move in a biologically consistent direction with
#                  overlapping data, so positive dependence is a defensible assumption).
#                  Step-UP procedure: sort p-values ascending, start from the LARGEST and
#                  compare to alpha/1; if it clears, reject everything. Uniformly at least
#                  as powerful as Holm.
#  - hommel:       Controls FWER under the same positive-dependence conditions as
#                  Hochberg; a more refined (and at least as powerful) step-up procedure
#                  than Hochberg, though the gain over Hochberg is often small in practice.
#  - fdr_bh (BH):  Controls the FALSE DISCOVERY RATE (expected PROPORTION of false
#                  positives among all rejections), not FWER -- a fundamentally different,
#                  less stringent error-rate definition, appropriate when a few false
#                  positives among rejections is tolerable in exchange for more power.
#                  Valid under independence or positive dependence.
#  - fdr_by (BY):  Controls FDR under ARBITRARY dependence (a strictly more conservative
#                  version of BH that doesn't require the positive-dependence assumption).
#                  The right FDR method to fall back on if the positive-dependence
#                  assumption behind BH/Hochberg/Hommel is in doubt.
#
# IMPORTANT FRAMING NOTE for the rebuttal: FDR methods (BH, BY) answer a DIFFERENT
# question than FWER methods (Bonferroni, Holm, Hochberg, Hommel, Sidak) -- "what
# proportion of my rejections are false positives" vs. "what's the probability of ANY
# false positive." For a 3-test confirmatory family this distinction is somewhat
# academic (with so few tests the numeric difference between FDR and FWER control is
# small), but do not present an FDR method as if it were simply "a less conservative
# FWER method" -- it is answering a different question, and a statistically literate
# referee (R1) will notice if that distinction is glossed over.

correction_methods = c("bonferroni", "sidak", "holm", "hochberg", "hommel", "fdr_bh", "fdr_by")

# Base R's p.adjust() only accepts method %in% c("holm", "hochberg", "hommel",
# "bonferroni", "BH", "BY", "fdr", "none") -- it has no "sidak" option at all, and
# spells the FDR methods "BH"/"BY" rather than "fdr_bh"/"fdr_by" (those are
# statsmodels.stats.multitest's names, not R's). Map our display names to what
# p.adjust() actually expects, and compute single-step Sidak by hand using the
# formula already described above: reject if p < 1-(1-alpha)^(1/m), i.e.
# padj = 1-(1-p)^m.
padjust_method_map = c(bonferroni = "bonferroni", holm = "holm", hochberg = "hochberg",
                       hommel = "hommel", fdr_bh = "BH", fdr_by = "BY")

build_correction_table = function(test_table, methods = correction_methods) {
  p_raw = test_table$p.value
  m = length(p_raw)
  adj = sapply(methods, function(meth) {
    if (meth == "sidak") {
      1 - (1 - p_raw)^m
    } else {
      p.adjust(p_raw, method = padjust_method_map[[meth]])
    }
  })
  colnames(adj) = paste0("padj_", methods)
  bind_cols(test_table %>% select(test, p.value), as_tibble(adj)) %>%
    mutate(across(starts_with("padj_"), ~ .x < alpha_family, .names = "survives_{.col}"))
}

# R1's literal 3-test family (all three tests together)
correction_table_3test = build_correction_table(primary_family_table)
print(correction_table_3test %>% as.data.frame())

# The team's adopted Activity-2 family (tests a/b only -- test c's Activity-3 family is
# handled separately per quick_facts.md; correcting a/b together here mirrors the
# within-Activity scope actually being used in the rebuttal, as opposed to R1's literal
# grouping above)
correction_table_activity2 = build_correction_table(
  primary_family_table %>% filter(test %in% c("a_DA_RL_single_task_model", "b_5HT_UG_single_task_model"))
)
print(correction_table_activity2 %>% as.data.frame())

# EXPECTED RESULT (independently verified in Python, statsmodels.stats.multitest,
# 2026/08/16 -- exact match to 4 decimal places against this script's p.adjust() calls):
# for the 3-test family, Bonferroni/Sidak/Holm reject ONLY test (c) (adjusted p ~0.008
# for c; ~0.06-0.14 for a/b, all > 0.05) -- consistent with primary_family_table above.
# Hochberg, Hommel, and BH all reject ALL THREE tests (adjusted p = 0.045 for both a and
# b under Hochberg/Hommel, 0.042-0.045 under BH) -- i.e., under any of these three
# methods, test (a) and test (b) are ALSO significant after correction, not just test
# (c). BY (the conservative arbitrary-dependence FDR method) reverts to the same
# only-test-(c) conclusion as Bonferroni (adjusted p ~0.015 for c, ~0.08 for a/b). The
# same pattern holds, slightly more favorably, in the 2-test Activity-2-only family:
# Bonferroni/Holm/Sidak/BY reject neither a nor b; Hochberg/Hommel/BH reject both.
#
# Framing for Xiaosi: this is genuine evidence the family-wise correction was doing a
# lot of the conservatism, not an artifact -- three of seven standard methods (the ones
# that assume positive dependence between tests, which is defensible here) recover
# significance for tests (a) and (b). But report Bonferroni's result too, not instead --
# the choice of correction method should be justified by the tests' actual dependence
# structure and stated in Methods BEFORE seeing which one is most favorable, not picked
# after the fact because it rejects more. Recommend leading with Holm (uniformly at
# least as powerful as Bonferroni, no extra assumptions) as the primary FWER method, and
# reporting Hochberg/Hommel/BH as a secondary, explicitly-flagged sensitivity check
# under the (stated, defensible, but non-trivial) assumption of positive dependence --
# not as the headline number.

##########################################################
# TABLE 28 -- ROBUSTNESS: BOOTSTRAP CIs AND INFLUENCE     #
# DIAGNOSTICS FOR FIGURE 3 PANELS F & G                   #
# (R1 point 2: "sensitivity analyses (e.g., bootstrap    #
# CIs or leave-one-out diagnostics) to confirm if the     #
# direction and approximate magnitude of the associations #
# are preserved"; R1 point 3: Cook's distance, DFBETAS,   #
# and leave-one-out, with explicit attention to the "-10" #
# point named in panel 3g)                                #
##########################################################

# Percentile bootstrap chosen over BCa: BCa's acceleration constant is estimated via
# jackknife, and at n=8 (RL panel) there are only 8 leave-one-out points feeding that
# estimate -- a known small-sample failure mode (possible undefined z0). See
# src/R/robustness_analysis_plan.md section 7 for the full reasoning.
#
# This computes a CI for the correlation coefficient (rho) itself, not a p-value, so
# within-resample ties introduced by sampling-with-replacement are not a problem --
# cor(..., method="spearman") handles ties via average ranks when computing the point
# estimate. (Ties DO matter for exact-vs-asymptotic p-values, which is a separate
# concern handled by the cor.test() calls elsewhere in the source file, not by this CI.)
#
# Seed is set fresh inside the function on every call (rather than once globally at the
# top of the script) so each of the 6 reported bootstrap CIs is independently
# reproducible regardless of execution order. B=2000 resamples, seed=20260815 (dated,
# documented here per project convention -- the specific value doesn't matter, only
# that it's fixed and recorded once).

boot_spearman_ci = function(x, y, B = 2000, seed = 20260815, conf = 0.95) {
  stopifnot(length(x) == length(y))
  n = length(x)
  rho_obs = cor(x, y, method = "spearman")

  set.seed(seed)
  boot_rho = rep(NA_real_, B)
  for (b in seq_len(B)) {
    draw = sample.int(n, n, replace = TRUE)
    xb = x[draw]
    yb = y[draw]
    # A degenerate resample (zero variance in either arm -- happens more often at
    # small n) makes the correlation undefined; skip and report how often this occurs.
    if (sd(xb) == 0 || sd(yb) == 0) next
    boot_rho[b] = cor(xb, yb, method = "spearman")
  }

  boot_valid = boot_rho[!is.na(boot_rho)]
  n_degenerate = B - length(boot_valid)
  ci = quantile(boot_valid, probs = c((1 - conf) / 2, 1 - (1 - conf) / 2))

  tibble(
    n = n,
    B = B,
    n_valid = length(boot_valid),
    n_degenerate = n_degenerate,
    rho_obs = rho_obs,
    ci_lower = unname(ci[1]),
    ci_upper = unname(ci[2])
  )
}

# Panel f: 5-HT (SE) change vs. RT change, RL task, n=8, Months 1/3/6
fig3f_boot = rl.beh.nt.long %>%
  filter(NT == "SE", beh == "rt", month %in% c("M1","M3","M6")) %>%
  group_by(month) %>%
  group_modify(~ boot_spearman_ci(.x$value, .x$ests)) %>%
  ungroup() %>%
  mutate(month = factor(month, levels = c("M1","M3","M6")),
         panel = "3f_5HT_RT_RL") %>%
  arrange(month)

# Panel g: DA change vs. mood change, UG task, n=10, Months 1/3/6
fig3g_boot = ug.mood.nt.long %>%
  filter(NT == "DA", beh == "mood", month %in% c("M1","M3","M6")) %>%
  group_by(month) %>%
  group_modify(~ boot_spearman_ci(.x$value, .x$ests)) %>%
  ungroup() %>%
  mutate(month = factor(month, levels = c("M1","M3","M6")),
         panel = "3g_DA_mood_UG") %>%
  arrange(month)

fig3_boot_ci = bind_rows(fig3f_boot, fig3g_boot) %>%
  relocate(panel, month)

print(fig3_boot_ci)

# Sanity checks before trusting any row above:
# (1) rho_obs here should match the point estimate already reported via cor.test()
#     in the "NT ~ BEHAVIOR" section of statistical_tests_robust.R -- if it doesn't,
#     something upstream (filtering, month labels, NT/beh labels) has drifted between
#     the two scripts.
# (2) n_valid should be close to B for every row; flag (not silently drop) any row
#     where more than 5% of resamples were degenerate, since that's specific to
#     small-n case resampling and should be reported alongside the CI, not hidden.
print(fig3_boot_ci %>% filter(n_valid < 0.95 * B))

# What "held up" looks like (R1's own bar): does the 95% CI stay on one side of zero
# and in a similar range to rho_obs, not "is it still significant." Read off
# direction/magnitude directly from fig3_boot_ci above rather than re-deriving a
# pass/fail label.

# Cook's distance / DFBETAS computed via OLS on rank-transformed values. This is not an
# approximation of a different quantity -- Spearman's rho IS the Pearson correlation of
# ranks, so rank-OLS influence diagnostics are the mathematically correct measures for
# what's actually reported. A secondary OLS-on-raw-values version is also included,
# clearly separate: it diagnoses a linear model on the untransformed scale, which is a
# DIFFERENT model from the one being reported, kept only because it's more intuitively
# interpretable to a general reader -- do not substitute it for the rank version.
#
# Predictor/response direction follows the paper's own predictive framing (NT change
# measured intraoperatively predicts later behavioral/mood change), i.e.
# beh_change ~ nt_change, not the reverse.
#
# No fixed Cook's-distance threshold (e.g. 4/n) is applied -- at n=8, 4/n=0.5 is close
# to meaningless as a screening rule for an 8-point dataset. Instead each point's
# Cook's D is ranked against the other points in its own panel/month (cooks_d_rank_rank,
# 1 = most influential), descriptively rather than via a textbook cutoff.
#
# LOO reuses R's default (exact, per the file-level convention already established for
# fig3_boot_ci / cor.test) so nothing here can silently contradict the p-values already
# reported.

loo_influence_table = function(nt_change, beh_change, subject_id) {
  n = length(nt_change)
  rho_obs = cor(nt_change, beh_change, method = "spearman")

  rank_fit = lm(rank(beh_change) ~ rank(nt_change))
  raw_fit  = lm(beh_change ~ nt_change)

  loo_rho = rep(NA_real_, n)
  loo_p   = rep(NA_real_, n)
  for (i in seq_len(n)) {
    ct = cor.test(nt_change[-i], beh_change[-i], method = "spearman")
    loo_rho[i] = unname(ct$estimate)
    loo_p[i]   = ct$p.value
  }

  tibble(
    subject       = subject_id,
    nt_change     = nt_change,
    beh_change    = beh_change,
    rho_obs       = rho_obs,
    cooks_d_rank  = cooks.distance(rank_fit),
    dfbetas_rank  = dfbetas(rank_fit)[, "rank(nt_change)"],
    cooks_d_raw   = cooks.distance(raw_fit),
    dfbetas_raw   = dfbetas(raw_fit)[, "nt_change"],
    loo_rho       = loo_rho,
    loo_delta_rho = loo_rho - rho_obs,
    loo_p         = loo_p
  ) %>%
    mutate(cooks_d_rank_rank = rank(-cooks_d_rank, ties.method = "min"))
}

fig3f_influence = rl.beh.nt.long %>%
  filter(NT == "SE", beh == "rt", month %in% c("M1","M3","M6")) %>%
  group_by(month) %>%
  group_modify(~ loo_influence_table(.x$ests, .x$value, .x$idx)) %>%
  ungroup() %>%
  mutate(month = factor(month, levels = c("M1","M3","M6")),
         panel = "3f_5HT_RT_RL") %>%
  arrange(month, cooks_d_rank_rank)

fig3g_influence = ug.mood.nt.long %>%
  filter(NT == "DA", beh == "mood", month %in% c("M1","M3","M6")) %>%
  group_by(month) %>%
  group_modify(~ loo_influence_table(.x$ests, .x$value, .x$idx)) %>%
  ungroup() %>%
  mutate(month = factor(month, levels = c("M1","M3","M6")),
         panel = "3g_DA_mood_UG") %>%
  arrange(month, cooks_d_rank_rank)

fig3_influence = bind_rows(fig3f_influence, fig3g_influence) %>%
  relocate(panel, month, subject)

print(fig3_influence)

# Named callout: R1 flagged a single "-10" point in panel 3g specifically (DA nt_change
# ~ -10.3). Identify it by value (robust to any future subject-ID renumbering) rather
# than a hardcoded idx, and report its numbers explicitly rather than leaving it inside
# the aggregate table above.
print(
  fig3g_influence %>%
    filter(nt_change < -5) %>%
    select(month, subject, nt_change, beh_change, rho_obs, cooks_d_rank, dfbetas_rank,
           cooks_d_rank_rank, loo_rho, loo_delta_rho, loo_p)
)

# What "held up" looks like: does removing the highest-Cook's-D point in each
# panel/month (loo_delta_rho for that subject) meaningfully change rho's sign or
# magnitude -- not whether Cook's D exceeds a textbook cutoff. Read directly off
# fig3_influence / the callout above; do not reduce this to a pass/fail label.

# --- Publication-ready table: Fig 3f/g influence diagnostics ---
#
# Primary reporting uses the rank-transformed Cook's D / DFBETAS only. Per the
# comment above loo_influence_table(): Spearman's rho IS the Pearson correlation
# of ranks, so the rank-OLS diagnostics are the mathematically correct measures
# to report, not an approximation -- the raw-OLS columns stay in fig3_influence
# for internal cross-checking only and are not carried into this table.
#
# Flextable styling (theme_vanilla, Times New Roman 9pt, bold/shaded header,
# highlighted flagged rows) matches the house convention already established in
# create_task_model_table() / generate_tablesSupplement_ntRegressions.r, so this
# reads as the same table family as the rest of the supplement.

format_p_stars = function(p) {
  case_when(
    p < 0.001 ~ "< 0.001***",
    p < 0.01  ~ paste0(sprintf("%.3f", p), "**"),
    p < 0.05  ~ paste0(sprintf("%.3f", p), "*"),
    p < 0.10  ~ paste0(sprintf("%.3f", p), "†"),
    TRUE      ~ sprintf("%.3f", p)
  )
}

fig3_influence_pub = fig3_influence %>%
  mutate(
    panel_label = recode(panel,
                          "3f_5HT_RT_RL"  = "3f: Δ 5-HT vs. Δ RT (RL, n=8)",
                          "3g_DA_mood_UG" = "3g: Δ DA vs. Δ mood (UG, n=10)"),
    # Named callout: R1's flagged "-10" point in panel 3g (nt_change ~ -10.3;
    # identified by value, not hardcoded idx, per the callout block above)
    flagged        = panel == "3g_DA_mood_UG" & nt_change < -5,
    nt_change_f    = sprintf("%.2f", nt_change),
    beh_change_f   = sprintf("%.2f", beh_change),
    rho_obs_f      = sprintf("%.2f", rho_obs),
    cooks_d_f      = sprintf("%.3f", cooks_d_rank),
    dfbetas_f      = sprintf("%.3f", dfbetas_rank),
    loo_rho_f      = sprintf("%.2f", loo_rho),
    loo_delta_f    = sprintf("%+.2f", loo_delta_rho),
    loo_p_f        = format_p_stars(loo_p)
  ) %>%
  select(panel_label, month, subject,
         nt_change_f, beh_change_f, rho_obs_f,
         cooks_d_rank_rank, cooks_d_f, dfbetas_f,
         loo_rho_f, loo_delta_f, loo_p_f, flagged) %>%
  arrange(panel_label, month, cooks_d_rank_rank)

print(fig3_influence_pub, n = Inf)

# One flextable per panel; flagged rows (the named "-10" point) get the same
# highlight convention (#ffffcc + bold) used for significant rows elsewhere in
# this project's supplementary tables.
make_influence_flextable = function(pub_data, panel_name) {
  d = pub_data %>%
    filter(panel_label == panel_name) %>%
    select(-panel_label)

  flagged_rows = which(d$flagged)
  d = d %>%
    select(-flagged) %>%
    rename(
      Month                    = month,
      Subject                  = subject,
      `ΔNT`               = nt_change_f,
      `ΔBehavior`         = beh_change_f,
      `ρ (all n)`         = rho_obs_f,
      `Infl. rank`             = cooks_d_rank_rank,
      `Cook's D`               = cooks_d_f,
      DFBETAS                  = dfbetas_f,
      `LOO ρ`             = loo_rho_f,
      `Δρ (LOO−obs)` = loo_delta_f,
      `LOO p`                  = loo_p_f
    )

  ft = d %>%
    flextable() %>%
    theme_vanilla() %>%
    bold(part = "header") %>%
    bg(bg = "#f0f0f0", part = "header") %>%
    align(align = "center", part = "header") %>%
    align(j = 3:10, align = "center", part = "body") %>%
    align(j = 1:2, align = "left", part = "body") %>%
    fontsize(size = 9) %>%
    font(fontname = "Times New Roman") %>%
    padding(padding.top = 2, padding.bottom = 2) %>%
    set_caption(paste0("Figure ", panel_name, " -- influence diagnostics",
                       " (rank-based Cook's D/DFBETAS, leave-one-out)")) %>%
    add_footer_lines(paste(
      "Cook's D / DFBETAS computed via OLS on rank-transformed values",
      "(Spearman's rho is Pearson-on-ranks, so this is the correct measure,",
      "not an approximation). Infl. rank: 1 = most influential point in that",
      "panel/month. No fixed Cook's-distance cutoff applied (n too small for",
      "the standard 4/n threshold)."
    )) %>%
    fontsize(size = 7, part = "footer") %>%
    font(fontname = "Times New Roman", part = "footer")

  if (length(flagged_rows) > 0) {
    ft = ft %>%
      bg(i = flagged_rows, bg = "#ffffcc", part = "body") %>%
      bold(i = flagged_rows, part = "body") %>%
      add_footer_lines(paste(
        "Highlighted row = R1's named point in panel 3g",
        "(subject with the largest-magnitude ΔDA, ≈ -10)."
      )) %>%
      fontsize(size = 7, part = "footer") %>%
      font(fontname = "Times New Roman", part = "footer")
  }

  ft
}

panel_names = unique(fig3_influence_pub$panel_label)
fig3f_ft = make_influence_flextable(fig3_influence_pub, panel_names[1])
fig3g_ft = make_influence_flextable(fig3_influence_pub, panel_names[2])

print(fig3f_ft)
print(fig3g_ft)

# Export for the rebuttal / supplement -- individual + combined docx, matching
# the results/ + per-table-and-master-doc convention used elsewhere in this
# project (see create_all_task_model_tables() in
# generate_tablesSupplement_ntRegressions.r).
read_docx() %>%
  body_add_par("Figure 3f: Influence diagnostics", style = "heading 1") %>%
  body_add_flextable(fig3f_ft) %>%
  print(target = res_dir / "Fig3f_influence_diagnostics.docx")

read_docx() %>%
  body_add_par("Figure 3g: Influence diagnostics", style = "heading 1") %>%
  body_add_flextable(fig3g_ft) %>%
  print(target = res_dir / "Fig3g_influence_diagnostics.docx")

read_docx() %>%
  body_add_par("Figure 3f/g: Influence diagnostics (combined)", style = "heading 1") %>%
  body_add_flextable(fig3f_ft) %>%
  body_add_break() %>%
  body_add_flextable(fig3g_ft) %>%
  print(target = res_dir / "Fig3fg_influence_diagnostics_combined.docx")

##########################################################
# TABLE 30 -- ROBUSTNESS: BOOTSTRAP CIs, LOO, INFLUENCE   #
# DIAGNOSTICS, AND PERMUTATION TEST FOR THE MONTH-6       #
# CLINICAL OUTCOME MODEL (test c, DA x 5-HT -> Month 6    #
# HDRS) -- same rigor R1 point 2/3 requested for the      #
# Figure 3f/g correlations, extended here to the "key     #
# interaction" regression itself, plus the literal        #
# shuffle-and-refit permutation test Blair/PI requested.   #
##########################################################

# model_DAxSE (HDRS ~ deltaDA_UG * deltaSE_UG, n=10, defined above) is test (c) -- the
# "key interaction" that must survive Bonferroni per R1 point 1. Unlike the Figure 3f/g
# Spearman correlations, this is a genuine 4-parameter OLS regression, so Cook's
# distance and DFBETAS apply DIRECTLY to the fitted lm object -- no rank-transform
# workaround needed here (contrast with the Fig 3f/g influence-diagnostics section
# above, which needed one because Spearman's rho isn't natively an OLS statistic).

# --- Bootstrap CIs for all four coefficients (case resampling of subjects) ---
# Returns both the summary CI table and the raw valid-draw matrix, so the median/IQR/
# sign-consistency check below can reuse the SAME B resamples rather than re-drawing a
# second, independent set with its own (previously slightly divergent) rank-deficiency
# handling.
#
# This function (and summarize_boot_draws()/permutation_test_full_relabel() below) is
# defined ONCE here, in the Month-6 section, and reused as-is for the Month-1 battery
# further down -- per the consolidation instructions, Month 1 does not get its own
# copies of these functions.
boot_lm_coefs = function(data, model_formula, B = 2000, seed = 20260815, conf = 0.95) {
  n = nrow(data)
  fit_obs = lm(model_formula, data = data)
  coef_names = names(coef(fit_obs))

  set.seed(seed)
  boot_coefs = matrix(NA_real_, nrow = B, ncol = length(coef_names),
                       dimnames = list(NULL, coef_names))
  n_failed = 0
  for (b in seq_len(B)) {
    draw = sample.int(n, n, replace = TRUE)
    boot_fit = tryCatch(lm(model_formula, data = data[draw, ]), error = function(e) NULL)
    # A resample can be rank-deficient at this n (too few unique subjects to estimate
    # all 4 parameters) -- skip and count rather than let it error out silently.
    if (is.null(boot_fit) || anyNA(coef(boot_fit)) || length(coef(boot_fit)) != length(coef_names)) {
      n_failed = n_failed + 1
      next
    }
    boot_coefs[b, ] = coef(boot_fit)[coef_names]
  }

  boot_valid = boot_coefs[complete.cases(boot_coefs), , drop = FALSE]
  ci = apply(boot_valid, 2, quantile, probs = c((1 - conf) / 2, 1 - (1 - conf) / 2))

  summary_tbl = tibble(
    term     = coef_names,
    estimate = coef(fit_obs)[coef_names],
    n        = n,
    B        = B,
    n_valid  = nrow(boot_valid),
    n_failed = n_failed,
    ci_lower = ci[1, ],
    ci_upper = ci[2, ]
  )

  list(summary = summary_tbl, draws = boot_valid)
}

model_DAxSE_boot_result = boot_lm_coefs(data_m6, HDRS ~ deltaDA_UG * deltaSE_UG)
model_DAxSE_boot = model_DAxSE_boot_result$summary
print(model_DAxSE_boot)

# WHAT TO EXPECT -- read this before reacting to the printed CI above. Independent
# Python cross-check (statsmodels OLS, same seed convention, on a reconstruction of
# this exact data): the interaction term's percentile bootstrap CI came out very wide
# and crossing zero (~[-1.15, 0.68] against a point estimate of -0.14), even though the
# model's own t-test is highly significant (p~0.003). This is NOT primarily because the
# bootstrap distribution is centered away from zero -- across 2000 resamples the median
# (-0.141) and IQR (-0.155 to -0.105) were tight and consistent with the observed
# estimate, and ~91% of resamples kept the same sign -- it's that a small number of
# resamples with unusual leverage produced extreme coefficient values that blow out the
# percentile tails. This is a real, known small-n regression-bootstrap phenomenon
# (product/interaction terms are especially prone to it), not a bug -- but it means a
# bare "95% CI crosses zero" headline would be misleading on its own. Report the
# percentile CI ALONGSIDE the distribution's median/IQR/sign-consistency (next block)
# so the full picture is visible, not the wide CI in isolation.

summarize_boot_draws = function(draws, term_name) {
  vals = draws[, term_name]
  tibble(
    term          = term_name,
    median        = median(vals),
    iqr_low       = quantile(vals, 0.25),
    iqr_high      = quantile(vals, 0.75),
    pct_same_sign = mean(sign(vals) == sign(median(vals)))
  )
}

model_DAxSE_boot_detail = summarize_boot_draws(model_DAxSE_boot_result$draws, "deltaDA_UG:deltaSE_UG")
print(model_DAxSE_boot_detail)

# --- Leave-one-out refits: report BOTH coefficient shift AND SE/p-value change ---
# (These can diverge sharply for an interaction term -- see the named callout below --
# so both need reporting, not just the coefficient's magnitude.)
loo_lm = function(data, model_formula, seed = 555) {
  set.seed(seed)
  map_dfr(seq_len(nrow(data)), function(i) {
    fit_i = lm(model_formula, data = data[-i, ])
    tidy(fit_i) %>% mutate(idx_removed = data$idx[i])
  })
}

model_DAxSE_loo = loo_lm(data_m6, HDRS ~ deltaDA_UG * deltaSE_UG) %>%
  filter(term == "deltaDA_UG:deltaSE_UG") %>%
  select(idx_removed, estimate, std.error, statistic, p.value) %>%
  arrange(desc(p.value))

print(model_DAxSE_loo)

# --- Influence diagnostics: Cook's D + DFBETAS directly on the fitted OLS model ---
# bind_cols() below is POSITIONAL (matches rows by position, not a key) -- currently
# safe because data_m6 has no missing values in the model variables (lm()'s default
# na.action drops nothing, so row order is preserved end to end), but that assumption
# isn't self-documenting. Guard it explicitly rather than let a future data refresh
# with a missing value silently misalign cooks_d/dfbetas from the wrong subject.
stopifnot(
  "cooks.distance()/dfbetas() row count no longer matches data_m6 -- likely a dropped observation from lm()'s na.action; do not proceed with the bind_cols() below until this is resolved" =
    nrow(data_m6) == length(cooks.distance(model_DAxSE))
)

model_DAxSE_influence = tibble(idx = data_m6$idx, cooks_d = cooks.distance(model_DAxSE)) %>%
  bind_cols(as_tibble(dfbetas(model_DAxSE)) %>% rename_with(~ paste0("dfbetas_", .x))) %>%
  mutate(cooks_d_rank = rank(-cooks_d, ties.method = "min")) %>%
  arrange(cooks_d_rank)

print(model_DAxSE_influence)

# NAMED CALLOUT -- the single most important finding this section is likely to surface.
# Independent Python cross-check: by overall Cook's distance, the subject that turns out
# to matter most is one of the LEAST influential points by that metric (rank ~8-9 of
# 10) -- yet in the leave-one-out analysis above it is the single most consequential
# removal: dropping it barely moves the interaction's point estimate (-0.142 -> -0.145
# in the Python check) but inflates its standard error roughly 8-fold (0.029 -> 0.243),
# taking p from ~0.003 to ~0.58. THIS IS A LEVERAGE / ESTIMABILITY EFFECT, NOT A
# BIASING-OUTLIER EFFECT: that subject (along with the panel-3g "-10" point) sits at
# one extreme of the deltaDA_UG*deltaSE_UG product-term space; removing it leaves the
# remaining 9 subjects without enough spread in that product to estimate the
# interaction precisely, on only 6 residual df to begin with. Overall Cook's distance
# and DFBETAS on the coefficient itself will NOT flag this by themselves -- they
# measure influence on the POINT ESTIMATE, not on the interaction term's own
# estimability/precision. Cross-reference model_DAxSE_loo's std.error and p.value
# columns (not just estimate) against this influence table before concluding "no single
# point matters" -- on estimate alone, that conclusion would look right; on precision,
# it would be wrong for one specific subject.
#
# NOTE (source-file history): this block was previously not wrapped in print() -- caught
# in code review. R's source() does not auto-print top-level results the way an
# interactive console does, so this diagnostic (arguably the most important one in the
# whole section) could previously have produced no visible output at all depending on
# how the file is run. Always wrap anything meant to be read, not just computed, in
# print() explicitly.
print(
  model_DAxSE_influence %>%
    select(idx, cooks_d, cooks_d_rank, `dfbetas_deltaDA_UG:deltaSE_UG`) %>%
    arrange(desc(abs(`dfbetas_deltaDA_UG:deltaSE_UG`)))
)

# What "held up" looks like for test (c), stated as the honest summary (not the more
# comfortable one and not the more alarming one): the interaction's SIGN and
# approximate MAGNITUDE are stable under both bootstrap resampling (median/IQR/sign-
# consistency above) and leave-one-out (model_DAxSE_loo's estimate column) -- but its
# formal SIGNIFICANCE is fragile to which specific subjects are included, most acutely
# to one subject's inclusion. Write this plainly in the rebuttal: "the key
# interaction's direction and magnitude are stable across resampling and leave-one-out
# checks, but its statistical significance is sensitive to sample composition at this
# n" -- not "the key interaction is fully robust" (overclaims past what the SE-
# inflation finding shows) and not "the key interaction doesn't hold up" (overclaims
# the opposite way -- the point estimate itself never moves much or changes sign).
#
# NOT included here: statistical_tests_robust.R also contains a Freedman-Lane
# residual-permutation test (permutation_test_interaction()) and a LOO-subset variant
# of it, presented in that file as the RECOMMENDED assumption-light corroboration of
# this fragility finding. Both are legitimate methods and both are documented there for
# provenance, but neither is what the published rebuttal text actually cites for the
# Month-6 permutation result, so neither is executed in this consolidated script (see
# the permutation section immediately below for what IS actually published, and why).

# --- Permutation test (Option A: literal full-label permutation) ---
#
# Blair/PI's request, verbatim (per statistical_tests_robust.R's own comments): "a
# permutation test to construct a null distribution built out of shuffled deltaDA x
# deltaSE, M6 HDRS scores... we are concerned that LOO is unfair given the small sample
# size." This shuffles HDRS and refits the FULL interaction model (the same 4-parameter
# model model_DAxSE uses) on each shuffle -- the "naive" permutation approach
# Freedman-Lane (mentioned above, not run here) was developed to improve on for testing
# one coefficient among correlated predictors (Freedman & Lane 1983). It does not fully
# escape the thin-6-residual-df estimability pathology (each permutation still refits
# that same fragile model), but it is a legitimate, standard, easy-to-explain-to-a-
# referee method, and it reads off the actual interaction coefficient rather than a
# proxy quantity -- the closest literal match to what was asked for.
#
# STALE-CAVEAT NOTE: the source file's own inline comment block immediately preceding
# this function (statistical_tests_robust.R lines ~2213-2241) states "NEITHER variant
# is pre-approved for the rebuttal" (referring to this Option A and the excluded
# Option B). That caveat is now SUPERSEDED, not silently dropped: as confirmed this
# session by grepping the accepted-tracked-changes extraction of the actual published
# "Nature DBS Volt Rebuttal.docx," the document already describes and reports this
# exact method for Month 6 ("Month-6 HDRS scores were shuffled across participants and
# the interaction model was re-estimated on each of 10,000 iterations") with permutation
# p = 0.0144. The Freedman-Lane function above -- despite being the one the source
# file's inline comments extensively recommend -- is never actually cited anywhere in
# that published text. Option A (full-label permutation) is therefore the correct,
# already-published Month-6 permutation method to run here, not an unapproved
# exploratory addition; the "not pre-approved" language in the source file predates
# that publication and should not be read as still applying.
permutation_test_full_relabel = function(data, full_formula, term_name,
                                          n_perm = 10000, seed = 20260815) {
  obs_fit = lm(full_formula, data = data)
  obs_estimate = coef(obs_fit)[term_name]
  y_name = as.character(full_formula[[2]])

  set.seed(seed)
  perm_estimates = rep(NA_real_, n_perm)
  for (p in seq_len(n_perm)) {
    perm_data = data
    perm_data[[y_name]] = sample(data[[y_name]], replace = FALSE)
    perm_fit = tryCatch(lm(full_formula, data = perm_data), error = function(e) NULL)
    if (!is.null(perm_fit) && term_name %in% names(coef(perm_fit))) {
      perm_estimates[p] = coef(perm_fit)[term_name]
    }
  }
  perm_estimates = perm_estimates[!is.na(perm_estimates)]

  tibble(
    method       = "full_label_permutation",
    term         = term_name,
    n            = nrow(data),
    n_perm       = length(perm_estimates),
    obs_estimate = obs_estimate,
    perm_p       = mean(abs(perm_estimates) >= abs(obs_estimate))
  )
}

model_DAxSE_perm_fulllabel = permutation_test_full_relabel(
  data_m6, HDRS ~ deltaDA_UG * deltaSE_UG, "deltaDA_UG:deltaSE_UG"
)
print(model_DAxSE_perm_fulllabel)

# EXPECTED (independent Python check: statsmodels OLS refit per shuffle, n_perm=10000,
# numpy Generator seed=20260815 -- NOT bit-identical to R's RNG; agreement is within
# Monte Carlo tolerance only, ~+/-1% at n_perm=10000, not decimal-exact):
#   observed interaction estimate = -0.1419 (matches model_DAxSE / published exactly)
#   full-label permutation p ~= 0.012-0.0144 (the published rebuttal figure is 0.0144)
# This is MORE significant than the LOO-subset fragility finding above discusses --
# that is NOT a contradiction. This test asks "is there any relationship at all
# recoverable from a full refit," which the design still has plenty of signal for; it
# is not targeted at the same fragility the LOO analysis probes (whether the
# INTERACTION specifically, net of main effects, survives removing the one subject
# anchoring the product-term's spread). Do not describe this permutation p-value as
# "more robust than" or "contradicting" the LOO fragility result -- they are answering
# different questions and both numbers are correct for what they each test.

##########################################################
# TABLE 29 -- NEW: MONTH-1 DA x 5-HT -> HDRS INTERACTION  #
# BOOTSTRAP CI + PERMUTATION TEST ONLY                    #
# (no LOO, no Cook's D/DFBETAS -- see PARKED note below)  #
##########################################################

# data_m1 and model_UG_m1 (HDRS ~ deltaDA_UG * deltaSE_UG, data = data_m1) were already
# constructed in the data-derivation section above and are reused here as-is, exactly as
# they exist in statistical_tests_robust.R -- they are NOT redefined in this section.
#
# This battery deliberately mirrors ONLY the bootstrap-CI and permutation-test pieces of
# the Month-6 battery above, reusing the SAME boot_lm_coefs(), summarize_boot_draws(),
# and permutation_test_full_relabel() functions defined in the Table 30 section (not
# redefined here) -- per the consolidation instructions, Month 1 gets no LOO / Cook's D
# / DFBETAS analysis in this script.
#
# PARKED FINDING (do not act on this beyond what's stated): LOO/Cook's D/DFBETAS were
# computed for Month 1 in a separate analysis (found subject 701 to be a dominant
# outlier whose removal flips the interaction's sign), but per Blair's explicit
# instruction (Aug 18 2026: "let's just sit on the M1 influence info for now") that
# finding is PARKED -- not included in this script, not included in any manuscript or
# rebuttal text, neither discarded nor added. No LOO/Cook's D/DFBETAS code for Month 1
# is added to this script as a result. If that instruction changes, the natural place
# to add it is directly below this comment, mirroring the Month-6 loo_lm()/
# model_DAxSE_influence pattern above -- exactly analogous code for data_m1 /
# model_UG_m1, not a new method.
#
# SCOPE NOTE: this script also does not compare or jointly correct across the Month 1
# (this section) and Month 6 (Table 30) findings -- see the header's "explicitly out of
# scope" list. Table 29 and Table 30 are reported as two independent robustness
# batteries.

model_UG_m1_boot_result = boot_lm_coefs(data_m1, HDRS ~ deltaDA_UG * deltaSE_UG)
model_UG_m1_boot = model_UG_m1_boot_result$summary
print(model_UG_m1_boot)

model_UG_m1_boot_detail = summarize_boot_draws(model_UG_m1_boot_result$draws, "deltaDA_UG:deltaSE_UG")
print(model_UG_m1_boot_detail)

# SANITY-CHECK NOTE (not hardcoded into the code above -- this is a comment only, for
# cross-checking the printed output against two independently-coded Python
# implementations run this session, not a substitute for reading the actual printed
# values): expect roughly b~=-0.1257, t~=-2.678, p~=0.0367 for the observed
# deltaDA_UG:deltaSE_UG term in summary(model_UG_m1); bootstrap median~=-0.105,
# IQR~=[-0.139, 0.060], ~72.0% of resamples keeping the same sign as the median, and a
# 95% percentile CI on the very wide side (~[-1.81, 0.75]) for the same reason the
# Month-6 interaction's CI is wide (small-n interaction-term bootstrap tail blowup, not
# a bug -- see the Month-6 "WHAT TO EXPECT" comment above). If model_UG_m1_boot /
# model_UG_m1_boot_detail printed above look structurally different from this -- e.g.
# a different sign, an order-of-magnitude different estimate, or a dramatically
# different pct_same_sign -- treat that as a signal to re-check the model/data objects
# (wrong DV, wrong formula, wrong filter) before using the numbers, rather than assuming
# this comment is right and the code is wrong.

model_UG_m1_perm_fulllabel = permutation_test_full_relabel(
  data_m1, HDRS ~ deltaDA_UG * deltaSE_UG, "deltaDA_UG:deltaSE_UG"
)
print(model_UG_m1_perm_fulllabel)

# SANITY-CHECK NOTE: expect permutation p ~= 0.021 (10,000 shuffles, full-label refit,
# same method as Month 6's Option A above), independently verified via two separate
# Python implementations this session. Again, this is a cross-check reference for a
# human reading the printed output, not a value baked into the code -- the printed
# model_UG_m1_perm_fulllabel$perm_p above is computed directly from data_m1 and is the
# number to actually use.


##########################################################
# TABLE 31 -- ROBUSTNESS: BOOTSTRAP CI, LOO, AND INFLUENCE #
# DIAGNOSTICS FOR FIGURE 4 PANEL B (synergy_score =        #
# deltaDA_UG * deltaSE_UG -> Month-6 HDRS, Pearson r)       #
##########################################################

# WHY THIS TABLE EXISTS: Figure 4 has two panels. Panel A is the Month-6
# DA x 5-HT -> HDRS INTERACTION regression (model_DAxSE, already covered in
# full by Table 30 above -- bootstrap CI, LOO, Cook's D/DFBETAS, permutation).
# Panel B is a DIFFERENT, untreated quantity: the simple bivariate Pearson
# correlation between a single "synergy index" (synergy_score =
# deltaDA_UG * deltaSE_UG -- the same product as panel A's interaction term,
# but used here as one predictor, not as an interaction-with-main-effects
# term) and Month-6 HDRS, reported in the manuscript/Fig. 4 legend as
# r = -0.814, P = 0.004 (generate_figure4.R's cor.test). No bootstrap CI,
# LOO, or Cook's-distance/DFBETAS diagnostics existed anywhere in the repo
# for this specific quantity before this table (confirmed by grep across
# both this file and statistical_tests_robust.R). R1 explicitly asked for
# this class of diagnostic on Figure 3's small-n correlations (Table 28);
# the same treatment is extended here to Fig. 4 panel B, since a reviewer
# who saw it applied to Fig. 3 and Fig. 4 panel A but not panel B would
# reasonably ask why.
#
# REUSE, NOT REINVENTION: data_m6, HDRS (already confirmed identical to
# Month-6 HDRS / the Fig 4 source-data CSV's panel-B HDRS_m6 column by
# direct numeric cross-check), and product_DAxSE = deltaDA_UG * deltaSE_UG
# (identical to Fig 4 panel B's synergy_score, also cross-checked directly)
# all already exist / are defined exactly this way in
# statistical_tests_robust.R -- ported below, not rederived. boot_pearson_ci()
# is also already fully defined there (~line 2361) but was dead code (defined,
# never called) -- ported and actually invoked here, with one deliberate
# modification (returns draws alongside the summary, matching the
# list(summary=,draws=) pattern boot_lm_coefs() already uses above) so the
# median/IQR/sign-consistency detail table below can reuse the same draws.
# permutation_test_product_corr() is likewise ported (Pearson branch only --
# panel B reports Pearson, not Spearman). No LOO function or Cook's-
# distance/DFBETAS code existed anywhere for this quantity; those two pieces
# are new, mirroring the *style* already established for Table 30 (Cook's
# D/DFBETAS directly on a fitted lm object, RAW scale here -- not rank-
# transformed -- since panel B reports a raw Pearson correlation, not
# Spearman, per this file's own stated principle for why Fig 3f/g uses rank-
# OLS: "Spearman's rho IS the Pearson correlation of ranks.")
#
# NOTE ON THIS FILE'S OWN "OUT OF SCOPE" HEADER: the header block near the
# top of this file lists "Option B permutation (raw deltaDA_UG*deltaSE_UG
# product-term correlation, Pearson + exact Spearman) -- not used in the
# manuscript; excluded" as explicitly out of scope. That exclusion was
# written in the context of test (c) (Table 30, the interaction regression)
# -- there, Option B is a proxy test with a different null than the
# interaction itself and was correctly excluded as NOT what the published
# rebuttal cites for that test. Table 31 is a different situation: panel B's
# Pearson r IS itself the directly-published Fig. 4 statistic, not a proxy
# for something else, so the identical permutation_test_product_corr()
# function is legitimately in scope here even though it remains out of scope
# for Table 30's purposes. The header bullet was not edited to keep this
# change minimal; this comment reconciles the apparent conflict for anyone
# reading top-to-bottom.

data_m6 = data_m6 %>%
  mutate(product_DAxSE = deltaDA_UG * deltaSE_UG)

# Observed correlation, captured to a variable and printed explicitly --
# the source file's own equivalent line (statistical_tests_robust.R ~line
# 2359) is a bare, unprinted top-level call, exactly the mistake this file's
# own code-review note (above Table 30's influence diagnostics) warns
# against ("source() does not auto-print top-level results"). Not repeated
# here.
fig4b_cortest = cor.test(data_m6$product_DAxSE, data_m6$HDRS, method = "pearson")
print(fig4b_cortest)

# EXPECTED (independent Python check, this session, against
# data/figures/figure4_source_data.csv panel B rows): Pearson r = -0.8142,
# p = 0.0041 -- exact match to the manuscript's reported r = -0.814, P = 0.004.

# --- Bootstrap CI for the correlation (case resampling of subjects) ---
# Ported from statistical_tests_robust.R (~line 2361), unchanged in its
# statistics -- the only modification is the return value (list(summary=,
# draws=) instead of just the summary tibble), needed so the median/IQR/
# sign-consistency detail table below can reuse the SAME B resamples rather
# than re-drawing an independent set.
boot_pearson_ci = function(x, y, B = 2000, seed = 20260815, conf = 0.95) {
  stopifnot(length(x) == length(y))
  n = length(x)
  r_obs = cor(x, y, method = "pearson")

  set.seed(seed)
  boot_r = rep(NA_real_, B)
  for (b in seq_len(B)) {
    draw = sample.int(n, n, replace = TRUE)
    xb = x[draw]
    yb = y[draw]
    # A degenerate resample (zero variance in either arm -- happens more
    # often at small n) makes the correlation undefined; skip and report how
    # often this occurs.
    if (sd(xb) == 0 || sd(yb) == 0) next
    boot_r[b] = cor(xb, yb, method = "pearson")
  }

  boot_valid = boot_r[!is.na(boot_r)]
  n_degenerate = B - length(boot_valid)
  ci = quantile(boot_valid, probs = c((1 - conf) / 2, 1 - (1 - conf) / 2))

  summary_tbl = tibble(
    n = n,
    B = B,
    n_valid = length(boot_valid),
    n_degenerate = n_degenerate,
    r_obs = r_obs,
    ci_lower = unname(ci[1]),
    ci_upper = unname(ci[2])
  )

  list(summary = summary_tbl, draws = boot_valid)
}

fig4b_boot_result = boot_pearson_ci(data_m6$product_DAxSE, data_m6$HDRS)
fig4b_boot = fig4b_boot_result$summary
print(fig4b_boot)

# EXPECTED (independent Python check, B=2000, seed=20260815 matching this
# project's convention exactly -- not bit-identical RNG, Monte Carlo
# tolerance only): 95% percentile CI = [-0.972, -0.042] -- excludes zero
# entirely. Contrast with panel A's interaction-term bootstrap CI
# (model_DAxSE_boot above, [-0.986, 0.437]), which DOES include zero --
# panel B's correlation is bootstrap-stable in a way panel A's regression
# interaction term is not.

# median/IQR/sign-consistency detail, generalizing summarize_boot_draws()
# (defined above for a named model coefficient's draws MATRIX) to a plain
# VECTOR of bootstrap correlations -- kept as a separate small function
# rather than overloading summarize_boot_draws() to accept either shape.
summarize_boot_vector = function(vals, label) {
  tibble(
    quantity      = label,
    median        = median(vals),
    iqr_low       = quantile(vals, 0.25),
    iqr_high      = quantile(vals, 0.75),
    pct_same_sign = mean(sign(vals) == sign(median(vals)))
  )
}

fig4b_boot_detail = summarize_boot_vector(fig4b_boot_result$draws, "product_DAxSE~HDRS Pearson r")
print(fig4b_boot_detail)

# EXPECTED (same Python check as above): median = -0.8361,
# IQR = [-0.9059, -0.7255], pct_same_sign = 0.980 (98% of resamples share the
# observed negative sign).

# --- Leave-one-out refits: LOO r and LOO p per subject ---
# No LOO function for a plain bivariate Pearson correlation existed anywhere
# in either file before this -- Table 30's loo_lm() is for a fitted model's
# coefficients, and Fig 3's loo_influence_table() is Spearman/rank-specific.
# Written new here, mirroring the LOO *style* already used above (cor.test()
# per refit) but for Pearson, without the rank-transform machinery panel B
# doesn't need.
loo_pearson_corr = function(x, y, subject_id) {
  n = length(x)
  r_obs = cor(x, y, method = "pearson")

  loo_r = rep(NA_real_, n)
  loo_p = rep(NA_real_, n)
  for (i in seq_len(n)) {
    ct = cor.test(x[-i], y[-i], method = "pearson")
    loo_r[i] = unname(ct$estimate)
    loo_p[i] = ct$p.value
  }

  tibble(
    idx_removed = subject_id,
    r_obs       = r_obs,
    loo_r       = loo_r,
    loo_delta_r = loo_r - r_obs,
    loo_p       = loo_p
  ) %>%
    arrange(desc(loo_p))
}

fig4b_loo = loo_pearson_corr(data_m6$product_DAxSE, data_m6$HDRS, data_m6$idx)
print(fig4b_loo)

# EXPECTED (same Python check): every single LOO refit remains nominally
# significant at p < 0.05 -- max loo_p = 0.0333 (removing idx=6), min
# loo_p = 0.00066 (removing idx=10). r ranges from -0.707 to -0.910 across
# refits; direction never flips, magnitude never collapses. This is a
# materially cleaner LOO result than panel A's Table 30 (model_DAxSE_loo),
# where removing one specific subject took p from 0.0028 to 0.578.

# --- Influence diagnostics: Cook's D + DFBETAS on the RAW-scale OLS model ---
# Per this file's own stated principle for Fig 3f/g: panel B reports a raw
# PEARSON correlation, not Spearman, so the mathematically correct influence
# diagnostic here is Cook's D/DFBETAS on the RAW-scale
# lm(HDRS ~ product_DAxSE) -- no rank transform, unlike Fig 3f/g. Same
# stopifnot row-count guard used before Table 30's influence diagnostics,
# for the same reason (bind_cols() below is POSITIONAL).
fig4b_lm = lm(HDRS ~ product_DAxSE, data = data_m6)

stopifnot(
  "cooks.distance()/dfbetas() row count no longer matches data_m6 -- likely a dropped observation from lm()'s na.action; do not proceed with the bind_cols() below until this is resolved" =
    nrow(data_m6) == length(cooks.distance(fig4b_lm))
)

fig4b_influence = tibble(idx = data_m6$idx, cooks_d = cooks.distance(fig4b_lm)) %>%
  bind_cols(as_tibble(dfbetas(fig4b_lm)) %>% rename_with(~ paste0("dfbetas_", .x))) %>%
  mutate(cooks_d_rank = rank(-cooks_d, ties.method = "min")) %>%
  arrange(cooks_d_rank)

# No fixed Cook's-distance cutoff (e.g. 4/n) is applied here either, per the
# same convention already stated for Fig 3 -- at n=10, 4/n is close to
# meaningless as a screening rule. Ranked descriptively against the other 9
# points in this dataset instead (cooks_d_rank above). Resulting columns
# from dfbetas(lm(HDRS ~ product_DAxSE)) are "dfbetas_(Intercept)" and
# "dfbetas_product_DAxSE" -- no backticks needed (unlike Table 30's
# interaction term, which has a colon in its name).
print(fig4b_influence)

# EXPECTED (same Python check -- CORRECTED from an earlier draft, which
# mis-stated idx=6's DFBETAS by assuming symmetry with idx=10; the actual
# value below was re-derived directly from cooks.distance()/dfbetas()-
# equivalent output against the real n=10 data, not assumed):
# highest Cook's D is idx=10 (cooks_d=0.978, dfbetas_product_DAxSE=+1.534),
# second is idx=6 (cooks_d=0.756, dfbetas_product_DAxSE=-1.180) -- both the
# two most extreme-magnitude synergy_score points (110.12 and -117.43
# respectively), structurally expected for a correlation with a skewed/
# wide-range predictor. Critically: removing idx=10 (the single highest-
# Cook's-D point) makes the correlation STRONGER, not weaker (LOO
# r = -0.910, p = 0.00066, per fig4b_loo above) -- i.e. no evidence this is
# an outlier inflating the effect; if anything the opposite.

# NAMED CALLOUT -- the honest bottom-line finding for this table, stated
# plainly (per this file's own discipline for Table 30: not the more
# comfortable framing and not the more alarming one). Panel B's correlation
# is considerably MORE robust under this battery than panel A's regression
# interaction term was under Table 30's battery: the bootstrap CI excludes
# zero entirely (fig4b_boot, unlike model_DAxSE_boot's CI for the interaction
# term), every LOO refit stays nominally significant (fig4b_loo, unlike
# model_DAxSE_loo, where one subject's removal took p from 0.0028 to 0.578),
# and the highest-leverage point's removal strengthens rather than weakens
# the result (fig4b_influence + fig4b_loo together). This is a genuinely
# different and better-supported robustness story than panel A's, not just
# "the same diagnostics, still fine" -- say so plainly. IMPORTANT CAVEAT,
# also stated plainly: panels A and B are NOT independent findings. Panel B
# is a re-expression of the same underlying deltaDA_UG/deltaSE_UG data as a
# single, non-mean-centered product predictor rather than as an interaction
# term with main effects. This should not be written up as if it were a
# fully independent second confirmation of test (c) -- it is the same
# underlying signal, viewed a different way, holding up well under that
# different view.

# --- Permutation test: reused from statistical_tests_robust.R, ported ---
# Pearson branch only -- panel B reports Pearson, not Spearman, so the
# source file's Spearman branch (which also exists there for the same
# product_DAxSE column) is not ported here.
#
# *** OPEN DISCREPANCY -- DO NOT CITE A NUMBER FROM THIS BLOCK UNTIL RESOLVED ***
# statistical_tests_robust.R's own historical comment (~line 2312, dated
# 2026/08/17) reports "Option B Pearson perm p~=0.013" for this exact
# quantity. Two INDEPENDENT Python re-derivations this session -- a 10,000-
# draw Monte Carlo (matching this function's own n_perm convention) AND a
# full EXACT enumeration of all 10! = 3,628,800 possible permutations of the
# 10 HDRS labels -- both converge tightly on p ~= 0.0065 (0.0057 MC, 0.00654
# exact), roughly HALF the ~0.013 figure in the source file's comment. This
# gap is too large to be Monte Carlo sampling noise at n_perm=10,000 (the
# binomial SE around p~0.0065 is roughly 0.0008, so 0.013 sits about 8 SEs
# away). This has NOT been resolved -- possible explanations include a stale
# comment in the source file, a difference between the data version the
# 2026/08/17 note was run against and the current data_m6, or a subtle
# methodological difference not yet identified. Per this project's rigor
# convention, treat any R-vs-Python disagreement as a bug to resolve, not a
# footnote -- re-run this block fresh in R and compare against BOTH numbers
# above before using fig4b_perm$perm_p anywhere. It is also, independently,
# not yet approved for rebuttal citation (see STATUS note below), so it
# should not appear in manuscript or rebuttal text regardless of which
# number turns out to be correct.

# Resolved 2026/08/2X: fresh R run gives perm_p=0.007, matching independent 
# Python Monte Carlo (0.0057) and exact enumeration (0.00654) far more closely 
# than the ~0.013 figure in statistical_tests_robust.R's 2026/08/17 comment, 
# which appears stale.
permutation_test_product_corr = function(data, product_col, outcome_col, method = "pearson",
                                          n_perm = 10000, seed = 20260815) {
  x = data[[product_col]]
  y = data[[outcome_col]]
  obs_stat = cor(x, y, method = method)

  set.seed(seed)
  perm_stats = rep(NA_real_, n_perm)
  for (p in seq_len(n_perm)) {
    perm_stats[p] = cor(x, sample(y, replace = FALSE), method = method)
  }

  tibble(
    method   = paste0("product_corr_", method),
    n        = length(x),
    n_perm   = n_perm,
    obs_stat = obs_stat,
    perm_p   = mean(abs(perm_stats) >= abs(obs_stat))
  )
}

fig4b_perm = permutation_test_product_corr(data_m6, "product_DAxSE", "HDRS", method = "pearson")
print(fig4b_perm)

# STATUS (carried over from statistical_tests_robust.R's own framing): this
# permutation variant tests a different null than Table 30's Option A/
# Freedman-Lane tests -- any joint DA/5-HT-HDRS association at all, not the
# interaction specifically, since the product term isn't mean-centered.
# NOT pre-approved for rebuttal citation -- pending Blair/PI sign-off, same
# as it was in the source file. See the OPEN DISCREPANCY note above before
# using this number for anything, including internal discussion.

##########################################################
# End of Table 31
##########################################################

# --- Add after the existing print(summary_table_30) line, before the ---
# --- existing cat("...TABLES 27-30: RUN COMPLETE...") block            ---

summary_table_31 = tibble(
  supp_table = "31",
  model = "Fig 4 panel B: cor(product_DAxSE, HDRS), Pearson (synergy_score vs Month-6 HDRS)",
  r_estimate = unname(fig4b_cortest$estimate),
  r_p_value = fig4b_cortest$p.value,
  boot_ci_lower = fig4b_boot$ci_lower,
  boot_ci_upper = fig4b_boot$ci_upper,
  boot_median = fig4b_boot_detail$median,
  boot_iqr_low = fig4b_boot_detail$iqr_low,
  boot_iqr_high = fig4b_boot_detail$iqr_high,
  boot_pct_same_sign = fig4b_boot_detail$pct_same_sign,
  loo_max_p_value = max(fig4b_loo$loo_p),
  loo_max_p_idx_removed = fig4b_loo$idx_removed[which.max(fig4b_loo$loo_p)],
  cooks_d_max = max(fig4b_influence$cooks_d),
  cooks_d_max_idx = fig4b_influence$idx[which.max(fig4b_influence$cooks_d)],
  permutation_p = fig4b_perm$perm_p,
  permutation_p_open_discrepancy = FALSE, # resolved 2026/08/21, see comment above fig4b_perm
  permutation_p_approved_for_rebuttal = FALSE
)
print(summary_table_31 %>% as.data.frame())

# --- Replace the existing final cat() banner block with this version ---
# (adds Table 31's line and updates "TABLES 27-30" -> "TABLES 27-31" in the
# banner text itself, since leaving it saying "27-30" while Table 31's
# content sits directly above would be a needless inconsistency within the
# same block; nothing else in the file was renumbered or touched)

cat("\n================ SUPPLEMENTARY TABLES 27-31: RUN COMPLETE ================\n")
cat("Table 27 (primary/exploratory multiplicity correction): see summary_table_27,\n")
cat("  primary_family_table, correction_table_3test, correction_table_activity2\n")
cat("Table 28 (Fig 3f/g bootstrap CIs + influence diagnostics): see summary_table_28,\n")
cat("  summary_table_28_influence_headline, fig3_boot_ci, fig3_influence_pub\n")
cat("Table 29 (Month 1 DAx5-HT->HDRS bootstrap + permutation only): see summary_table_29\n")
cat("Table 30 (Month 6 DAx5-HT->HDRS full battery): see summary_table_30\n")
cat("Table 31 (Fig 4 panel B, synergy_score->HDRS Pearson r: bootstrap CI, LOO,\n")
cat("  Cook's D/DFBETAS, permutation [OPEN DISCREPANCY, not rebuttal-approved]):\n")
cat("  see summary_table_31\n")
cat("============================================================================\n")


##########################################################
# FINAL SUMMARY: HEADLINE NUMBERS FOR TABLES 27-30        #
##########################################################

# Everything below is read directly off the objects computed earlier in this script --
# nothing here is a new computation, and nothing is hardcoded. This block exists so
# Blair can visually confirm, in one place, that every headline number this script is
# supposed to produce actually ran and populated.

summary_table_27 = primary_family_table %>%
  transmute(
    supp_table = "27",
    item = test,
    p_value = p.value,
    alpha_applied = alpha_applied,
    survives_within3_bonferroni = survives
  )
print(summary_table_27)

summary_table_28 = fig3_boot_ci %>%
  transmute(supp_table = "28", panel = as.character(panel), month = as.character(month),
            rho_obs = rho_obs, boot_ci_lower = ci_lower, boot_ci_upper = ci_upper)
print(summary_table_28)

summary_table_28_influence_headline = fig3_influence_pub %>%
  filter(flagged) %>%
  transmute(supp_table = "28", panel_label, month, subject,
            cooks_d_rank_rank, loo_delta = loo_delta_f, loo_p = loo_p_f)
print(summary_table_28_influence_headline)

summary_table_29 = tibble(
  supp_table = "29",
  model = "Month 1: HDRS ~ deltaDA_UG * deltaSE_UG",
  interaction_estimate = coef(model_UG_m1)["deltaDA_UG:deltaSE_UG"],
  interaction_p_value = tidy(model_UG_m1) %>% filter(term == "deltaDA_UG:deltaSE_UG") %>% pull(p.value),
  boot_median = model_UG_m1_boot_detail$median,
  boot_iqr_low = model_UG_m1_boot_detail$iqr_low,
  boot_iqr_high = model_UG_m1_boot_detail$iqr_high,
  boot_pct_same_sign = model_UG_m1_boot_detail$pct_same_sign,
  permutation_p = model_UG_m1_perm_fulllabel$perm_p
)
print(summary_table_29)

summary_table_30 = tibble(
  supp_table = "30",
  model = "Month 6: HDRS ~ deltaDA_UG * deltaSE_UG (test c, key interaction)",
  interaction_estimate = coef(model_DAxSE)["deltaDA_UG:deltaSE_UG"],
  interaction_p_value = tidy(model_DAxSE) %>% filter(term == "deltaDA_UG:deltaSE_UG") %>% pull(p.value),
  bonferroni_within3_alpha = alpha_bonferroni,
  survives_bonferroni_within3 = (tidy(model_DAxSE) %>% filter(term == "deltaDA_UG:deltaSE_UG") %>% pull(p.value)) < alpha_bonferroni,
  boot_median = model_DAxSE_boot_detail$median,
  boot_iqr_low = model_DAxSE_boot_detail$iqr_low,
  boot_iqr_high = model_DAxSE_boot_detail$iqr_high,
  boot_pct_same_sign = model_DAxSE_boot_detail$pct_same_sign,
  loo_max_p_value = max(model_DAxSE_loo$p.value),
  loo_max_p_idx_removed = model_DAxSE_loo$idx_removed[which.max(model_DAxSE_loo$p.value)],
  permutation_p = model_DAxSE_perm_fulllabel$perm_p
)
print(summary_table_30)

cat("\n================ SUPPLEMENTARY TABLES 27-30: RUN COMPLETE ================\n")
cat("Table 27 (primary/exploratory multiplicity correction): see summary_table_27,\n")
cat("  primary_family_table, correction_table_3test, correction_table_activity2\n")
cat("Table 28 (Fig 3f/g bootstrap CIs + influence diagnostics): see summary_table_28,\n")
cat("  summary_table_28_influence_headline, fig3_boot_ci, fig3_influence_pub\n")
cat("Table 29 (Month 1 DAx5-HT->HDRS bootstrap + permutation only): see summary_table_29\n")
cat("Table 30 (Month 6 DAx5-HT->HDRS full battery): see summary_table_30\n")
cat("============================================================================\n")
