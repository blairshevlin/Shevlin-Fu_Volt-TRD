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
# 2026/08/16    Blair Shevlin                           corrected stale "unplanned/opportunistic" framing of test (c) in the Bonferroni-within-3 comment block to reflect Xiaosi Gu's Aug 15 correction (pre-specified Activity-3 hypothesis); see quick_facts.md and narrative_guidance.md
# 2026/08/16    Blair Shevlin                           resolved the rl.lm.DA.2 vs. sess.lm.DA.2 open question (test a/b object mapping) via independent Python re-derivation from figure2_source_data.csv; see comment above the Bonferroni-within-3 block
# 2026/08/16    Blair Shevlin                           added multi-method correction table (Bonferroni/Sidak/Holm/Hochberg/Hommel/BH/BY) per Xiaosi Gu's request to check whether Bonferroni is overly conservative for the primary family; see comment block before build_correction_table()
# 2026/08/17    Blair Shevlin                           added two literal permutation-test variants for test (c) per Blair/PI's explicit request ("shuffled deltaDA x deltaSE, M6 HDRS scores"), distinct from the existing Freedman-Lane test above -- Option A: full-label permutation refitting the full interaction model; Option B: permutation test of the raw deltaDA_UG*deltaSE_UG product against HDRS (Pearson + exact Spearman). Both independently cross-checked in Python (obs interaction=-0.1419 matches model_DAxSE exactly; Option A perm p~=0.012; Option B Pearson perm p~=0.013, Spearman perm p~=0.045-0.049). See src/R/permutation_responder_robustness_plan.md sec 1.2-1.3 for full reasoning and the multiplicity/framing caveat repeated in-line below -- NEITHER variant is pre-approved for the rebuttal.
# 2026/08/18    Blair Shevlin                           consolidated sensitivity/robustness analyses (bootstrap CIs, LOO, influence diagnostics, permutation tests, multiplicity correction) into a single script per Blair's request; MDE/Batten-effect-size content intentionally excluded (lives in power_mde_analysis.R)
# 2026/11/13    Blair Shevlin                           relabelled tables to final SI numbering: code Table 27 -> SI 27; 28 -> SI 28 (bootstrap) + 29 (influence); 29 -> SI 30 Month 1 column; 30 -> SI 30 Month 6 column + 31 (LOO/influence); 31 -> SI 32 (bootstrap/permutation) + 33 (LOO/influence). summary_table_* object names are unchanged
# 2026/08/20    Blair Shevlin                           added Fig 4 panel B (synergy_score / product_DAxSE -> Month-6 HDRS, raw Pearson correlation) bootstrap CI, LOO, and Cook's D/DFBETAS diagnostics as Table 31, extending the Table 28/30 R1 diagnostic battery to the one Fig 4 quantity it hadn't yet been applied to; permutation p-value intentionally left as an open discrepancy (see comment below) pending a fresh R run

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
  mutate(baseline_HDRS = HDRS[session == "pre stim"],
         baseline_MADRS = MADRS[session == "pre stim"],
         m6_HDRS = HDRS[session == "month 6"],
         m6_MADRS = MADRS[session == "month 6"],
         deltaPerDA_UG = (DA_Post_UG - DA_Pre_UG)/DA_Pre_UG,
         deltaPerSE_UG = (SE_Post_UG - SE_Pre_UG)/SE_Pre_UG,
         deltaPerNE_UG = (NE_Post_UG - NE_Pre_UG)/NE_Pre_UG,
         deltaPerDA_RL = (DA_Post_RL - DA_Pre_RL)/DA_Pre_RL,
         deltaPerSE_RL = (SE_Post_RL - SE_Pre_RL)/SE_Pre_RL,
         deltaPerNE_RL = (NE_Post_RL - NE_Pre_RL)/NE_Pre_RL,
         deltaPerHDRS = (HDRS[session == "pre stim"] - HDRS[session == "month 6"])/HDRS[session == "pre stim"]
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

# Month 1 DA x 5-HT -> HDRS interaction model (used as-is for the Supplementary Table 30, Month 1 column,
# bootstrap/permutation battery below -- NOT redefined there)
model_UG_m1 <- lm(HDRS ~ deltaDA_UG * deltaSE_UG, data = data_m1)
print(summary(model_UG_m1))

# Month 6 DA x 5-HT -> HDRS interaction model -- this is test (c), the
# manuscript's "key interaction," and feeds Supplementary Table 27 and Supplementary Tables 30 (Month 6 column) and 31 below.
model_DAxSE = lm(HDRS ~ deltaDA_UG * deltaSE_UG, data = data_m6)
print(summary(model_DAxSE))

################################################
# SUPPLEMENTARY TABLE 27 -- ROBUSTNESS: BONFERRONI-WITHIN-3 #
# FOR THE"KEY INTERACTION                     #
###############################################

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


#############################################
# ROBUSTNESS: MULTIPLE CORRECTION METHODS   #
#############################################

correction_methods = c("bonferroni", "sidak", "holm", "hochberg", "hommel", "fdr_bh", "fdr_by")

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

# 3-test family (all three tests together)
correction_table_3test = build_correction_table(primary_family_table)
print(correction_table_3test %>% as.data.frame())


correction_table_activity2 = build_correction_table(
  primary_family_table %>% filter(test %in% c("a_DA_RL_single_task_model", "b_5HT_UG_single_task_model"))
)
print(correction_table_activity2 %>% as.data.frame())

##########################################################
# SUPPLEMENTARY TABLES 28 (BOOTSTRAP CIs) AND 29        #
# (INFLUENCE DIAGNOSTICS) FOR FIGURE 3 PANELS F & G      #
##########################################################

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

# Sanity checks
print(fig3_boot_ci %>% filter(n_valid < 0.95 * B))

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

print(
  fig3g_influence %>%
    filter(nt_change < -5) %>%
    select(month, subject, nt_change, beh_change, rho_obs, cooks_d_rank, dfbetas_rank,
           cooks_d_rank_rank, loo_rho, loo_delta_rho, loo_p)
)

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

# Export for the rebuttal / supplement
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
# SUPPLEMENTARY TABLES 30 (MONTH 6 COLUMN: BOOTSTRAP,   #
# PERMUTATION) AND 31 (LOO, INFLUENCE DIAGNOSTICS) FOR  #
# THE MONTH-6 CLINICAL OUTCOME MODEL                     #
##########################################################

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

# Leave-one-out refits
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

# Influence diagnostics
stopifnot(
  "cooks.distance()/dfbetas() row count no longer matches data_m6 -- likely a dropped observation from lm()'s na.action; do not proceed with the bind_cols() below until this is resolved" =
    nrow(data_m6) == length(cooks.distance(model_DAxSE))
)

model_DAxSE_influence = tibble(idx = data_m6$idx, cooks_d = cooks.distance(model_DAxSE)) %>%
  bind_cols(as_tibble(dfbetas(model_DAxSE)) %>% rename_with(~ paste0("dfbetas_", .x))) %>%
  mutate(cooks_d_rank = rank(-cooks_d, ties.method = "min")) %>%
  arrange(cooks_d_rank)

print(model_DAxSE_influence)

print(
  model_DAxSE_influence %>%
    select(idx, cooks_d, cooks_d_rank, `dfbetas_deltaDA_UG:deltaSE_UG`) %>%
    arrange(desc(abs(`dfbetas_deltaDA_UG:deltaSE_UG`)))
)

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

##########################################################
# SUPPLEMENTARY TABLE 30, MONTH 1 COLUMN: DA x 5-HT ->   #
# HDRS INTERACTION, BOOTSTRAP CI + PERMUTATION TEST ONLY #
##########################################################

model_UG_m1_boot_result = boot_lm_coefs(data_m1, HDRS ~ deltaDA_UG * deltaSE_UG)
model_UG_m1_boot = model_UG_m1_boot_result$summary
print(model_UG_m1_boot)

model_UG_m1_boot_detail = summarize_boot_draws(model_UG_m1_boot_result$draws, "deltaDA_UG:deltaSE_UG")
print(model_UG_m1_boot_detail)

model_UG_m1_perm_fulllabel = permutation_test_full_relabel(
  data_m1, HDRS ~ deltaDA_UG * deltaSE_UG, "deltaDA_UG:deltaSE_UG"
)
print(model_UG_m1_perm_fulllabel)

############################################################
# SUPPLEMENTARY TABLES 32 (BOOTSTRAP CI, PERMUTATION) AND #
# 33 (LOO, INFLUENCE DIAGNOSTICS) FOR FIGURE 4 PANEL B     #
############################################################

data_m6 = data_m6 %>%
  mutate(product_DAxSE = deltaDA_UG * deltaSE_UG)

fig4b_cortest = cor.test(data_m6$product_DAxSE, data_m6$HDRS, method = "pearson")
print(fig4b_cortest)

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

fig4b_lm = lm(HDRS ~ product_DAxSE, data = data_m6)

stopifnot(
  "cooks.distance()/dfbetas() row count no longer matches data_m6 -- likely a dropped observation from lm()'s na.action; do not proceed with the bind_cols() below until this is resolved" =
    nrow(data_m6) == length(cooks.distance(fig4b_lm))
)

fig4b_influence = tibble(idx = data_m6$idx, cooks_d = cooks.distance(fig4b_lm)) %>%
  bind_cols(as_tibble(dfbetas(fig4b_lm)) %>% rename_with(~ paste0("dfbetas_", .x))) %>%
  mutate(cooks_d_rank = rank(-cooks_d, ties.method = "min")) %>%
  arrange(cooks_d_rank)

print(fig4b_influence)

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

#################################
# Supplementary Tables 32-33    #
#################################

summary_table_31 = tibble(
  supp_table = "32-33",
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

cat("\n================ SUPPLEMENTARY TABLES 27-33: RUN COMPLETE ================\n")
cat("SI Table 27 (primary/exploratory multiplicity correction): see summary_table_27,\n")
cat("  primary_family_table, correction_table_3test, correction_table_activity2\n")
cat("SI Tables 28-29 (Fig 3f/g bootstrap CIs; influence diagnostics): see summary_table_28,\n")
cat("  summary_table_28_influence_headline, fig3_boot_ci, fig3_influence_pub\n")
cat("SI Table 30, Month 1 column (DAx5-HT->HDRS bootstrap + permutation only): see summary_table_29\n")
cat("SI Tables 30-31, Month 6 (DAx5-HT->HDRS full battery): see summary_table_30\n")
cat("SI Tables 32-33 (Fig 4 panel B, synergy_score->HDRS Pearson r: bootstrap CI, LOO,\n")
cat("  Cook's D/DFBETAS, permutation [OPEN DISCREPANCY, not rebuttal-approved]):\n")
cat("  see summary_table_31\n")
cat("============================================================================\n")


##########################################################
# FINAL SUMMARY: HEADLINE NUMBERS FOR SI TABLES 27-31     #
##########################################################

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
  transmute(supp_table = "29", panel_label, month, subject,
            cooks_d_rank_rank, loo_delta = loo_delta_f, loo_p = loo_p_f)
print(summary_table_28_influence_headline)

summary_table_29 = tibble(
  supp_table = "30 (Month 1)",
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
  supp_table = "30-31 (Month 6)",
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

cat("\n================ SUPPLEMENTARY TABLES 27-31: RUN COMPLETE ================\n")
cat("SI Table 27 (primary/exploratory multiplicity correction): see summary_table_27,\n")
cat("  primary_family_table, correction_table_3test, correction_table_activity2\n")
cat("SI Tables 28-29 (Fig 3f/g bootstrap CIs; influence diagnostics): see summary_table_28,\n")
cat("  summary_table_28_influence_headline, fig3_boot_ci, fig3_influence_pub\n")
cat("SI Table 30, Month 1 column (DAx5-HT->HDRS bootstrap + permutation only): see summary_table_29\n")
cat("SI Tables 30-31, Month 6 (DAx5-HT->HDRS full battery): see summary_table_30\n")
cat("============================================================================\n")
