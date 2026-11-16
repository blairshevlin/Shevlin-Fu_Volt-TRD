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
# 2025/05/05    Blair Shevlin                           wrote original code (statistical_tests.R)
# 2025/07/24    Blair Shevlin                           updated to use new NT data
# 2025/09/16    Blair Shevlin                           updated to use revised NT data
# 2025/11/18    Blair Shevlin                           added task x sess analysis
# 2026/05/05    Blair Shevlin                           added HDRS change analysis
# 2026/08/18    Blair Shevlin                           added sensitivity and robustness tests (sensitivity_and_robustness_tests.R)
# 2026/11/16    Blair Shevlin                           merged statistical_tests.R and sensitivity_and_robustness_tests.R into
#                                                       main_analyses.R; kept only analyses reported in the main text, figure
#                                                       legends and Methods; console output only; fixed the Family 1 test
#                                                       selection (term "stimPost-Stim"), which had silently dropped tests (a)
#                                                       and (b) and stopped the multiplicity table with an error
# 2026/11/16    Blair Shevlin                           session coded pre = -1, post = +1 in every model (was 0/1),
#                                                       so stimulation coefficients are half the post-minus-pre difference
#
# Reproduces every statistic reported in the main text, the main and Extended Data figure legends and the
# Methods of Shevlin, Fu et al., including the full output of the Supplementary Tables those numbers come
# from. Analyses reported only in the Supplementary Information are in supplementary_information.R.
#
#   Section 1  Clinical outcome ................................ Results "Study Procedure and Clinical Outcomes"; Fig. 1e
#   Section 2  Stimulation effects on DA / 5-HT / NE ........... Fig. 2a-b, Supplementary Tables 3-4 (Family 1 tests)
#              Task x stimulation interaction .................. Supplementary Table 5
#              Block type / offer amount ANOVAs ................ Fig. 2 legend (block type, offer amount)
#              Norepinephrine .................................. Extended Data Fig. 3
#   Section 3  Longitudinal behaviour .......................... Fig. 3a-e
#   Section 4  NT change vs. behavioural change ................ Fig. 3f-g
#   Section 5  NT change vs. HDRS-17 ........................... Fig. 4a, Supplementary Table 21
#   Section 6  Multiplicity correction ......................... Supplementary Table 27; Methods "Statistical analysis"
#   Section 7  Robustness: Fig. 3f/g correlations .............. Supplementary Tables 28-29
#   Section 8  Robustness: DA x 5-HT -> HDRS interaction ....... Supplementary Tables 30-31
#   Section 9  Robustness: Fig. 4b product correlation ......... Supplementary Tables 32-33
#
# Run time: under 1 minute.

rm(list = ls())
options(width = 200, pillar.sigfig = 5)

library(tidyverse)
library(fs)
library(here)
library(lmerTest)
library(broom)
library(broom.mixed)
library(effectsize)

# Optional: set export_docx = TRUE to also write Supplementary Tables 3-5 as Word tables in
# results/si_tables/ (Supplementary Information format, ready to copy). Off by default: no files written.
export_docx = FALSE
if (export_docx) source(here::here("src", "R", "write_si_table_docx.R"))

# Paths
dir = path(here())
nt_dir = dir / "data" / "nt" / "processed"
beh_dir = dir / "data" / "behavior"
clin_dir = dir / "data" / "clinical"

#################
# DATA          #
#################

# Subject-level NT Data
load(nt_dir / "UG_RL_NT-Continuous_9-23-25.RData")

# Trial-level estimates at outcome reveal (RL) and offer presentation (UG)
rl.EST.Reward = rl.EST %>%
  filter(event == "Reward") %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))

ug.EST.Offer = ug.EST %>%
  filter(event == "Offer") %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))

# Load clinical data
cl.df = read.csv(clin_dir / "clinical-data_deid_07-10-25.csv")

# Load behavioral data
rl.beh = read.csv(file = beh_dir / "rl-data_deid_07-10-25.csv")
ug.beh = read.csv(file = beh_dir / "ug-data_deid_07-10-25.csv")

rl.EST.pr = rl.EST.Reward %>%
  mutate(block_type = factor(cond, levels = c("Mixed","Negative","Positive")),
         trial_z = scale(trial)[,1])

ug.EST.pr = ug.EST.Offer %>%
  mutate(offer_bin_f = factor(offer_bin, levels = c("Middle","Low","High")),
         trial_z = scale(trial)[,1])

task.merge = rbind(rl.EST.pr %>% select(idx,stim,nt,Oz,trial_z) %>% mutate(task = "RL"),
                   ug.EST.pr %>% select(idx,stim,nt,Oz,trial_z) %>% mutate(task = "UG")) %>%
  mutate(task = factor(task))
contrasts(task.merge$stim) <- c(-1, 1)  # session coding: pre-stimulation = -1, post-stimulation = +1

# Session-level behaviour (trials with RT < 10 s); the preoperative "fMRI" session is the baseline
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

# Stimulation-induced NT change per participant (post minus pre, mean Oz)
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

# Behavioural change from baseline, joined with NT change (long format: one row per NT x measure x session)
rl.beh.nt.long = rl.beh.means %>%
  mutate(mRT = mLogRT) %>%
  select(!mLogRT) %>%
  pivot_wider(names_from = sess, values_from = c(mChoice, mRT)) %>%
  mutate(mChoice_M1_Change = `mChoice_Month 1` - mChoice_Baseline,
         mChoice_M3_Change = `mChoice_Month 3` - mChoice_Baseline,
         mChoice_M6_Change = `mChoice_Month 6` - mChoice_Baseline,
         mRT_M1_Change = `mRT_Month 1` - mRT_Baseline,
         mRT_M3_Change = `mRT_Month 3` - mRT_Baseline,
         mRT_M6_Change = `mRT_Month 6` - mRT_Baseline) %>%
  select(idx, ends_with("_Change")) %>%
  merge(rl.nt) %>%
  pivot_longer(cols = c(DA,SE), names_to = "NT", values_to = "ests") %>%
  pivot_longer(cols = ends_with("_Change"), names_to = "change") %>%
  mutate(month = str_extract(change, "(?<=_)[^_]+(?=_)"),
         beh = ifelse(grepl("mChoice", change), "choice", "rt"))

ug.mood.nt.long = ug.mood.means %>%
  mutate(mRT = mLogRT) %>%
  select(!mLogRT) %>%
  pivot_wider(names_from = sess, values_from = c(mMood, mRT)) %>%
  mutate(mMood_M1_Change = `mMood_Month 1` - mMood_Baseline,
         mMood_M3_Change = `mMood_Month 3` - mMood_Baseline,
         mMood_M6_Change = `mMood_Month 6` - mMood_Baseline) %>%
  select(idx, ends_with("_Change")) %>%
  merge(ug.nt) %>%
  pivot_longer(cols = c(DA,SE), names_to = "NT", values_to = "ests") %>%
  pivot_longer(cols = ends_with("_Change"), names_to = "change") %>%
  mutate(month = str_extract(change, "(?<=_)[^_]+(?=_)"),
         beh = ifelse(grepl("mMood", change), "mood", "rt"))

# NT change (mean Oz, post minus pre) joined with HDRS-17
nt_change_wide = function(est, task_name) {
  est %>%
    group_by(idx,stim,nt) %>%
    summarise(mTrial = mean(Oz), .groups = "drop") %>%
    pivot_wider(values_from = mTrial, names_from = c("stim","nt")) %>%
    transmute(idx,
              deltaDA = `Post-Stim_DA` - `Pre-Stim_DA`,
              deltaSE = `Post-Stim_SE` - `Pre-Stim_SE`,
              deltaNE = `Post-Stim_NE` - `Pre-Stim_NE`) %>%
    rename_with(~ paste0(.x, "_", task_name), -idx)
}

nt.cl = merge(nt_change_wide(ug.EST.Offer, "UG"), nt_change_wide(rl.EST.Reward, "RL"), all.x = TRUE)

cl.Oz.lme = merge(cl.df, nt.cl) %>%
  filter(!session %in% c("post stim")) %>%
  group_by(idx) %>%
  mutate(baseline_HDRS = HDRS[session == "pre stim"],
         deltaPerHDRS = (HDRS[session == "pre stim"] - HDRS[session == "month 6"])/HDRS[session == "pre stim"]) %>%
  filter(!session %in% c("pre stim","fmri")) %>%
  ungroup()

data_m1 = cl.Oz.lme %>% filter(session == "month 1") %>%
  select(idx,HDRS,baseline_HDRS,deltaPerHDRS,deltaDA_UG,deltaSE_UG,deltaNE_UG,deltaDA_RL,deltaSE_RL,deltaNE_RL)
data_m3 = cl.Oz.lme %>% filter(session == "month 3") %>%
  select(idx,HDRS,baseline_HDRS,deltaPerHDRS,deltaDA_UG,deltaSE_UG,deltaNE_UG,deltaDA_RL,deltaSE_RL,deltaNE_RL)
data_m6 = cl.Oz.lme %>% filter(session == "month 6") %>%
  select(idx,HDRS,baseline_HDRS,deltaPerHDRS,deltaDA_UG,deltaSE_UG,deltaNE_UG,deltaDA_RL,deltaSE_RL,deltaNE_RL)

####################################################
# SECTION 1: CLINICAL OUTCOME (Fig. 1e)            #
####################################################
cat("\n===== Section 1: HDRS-17, baseline vs. month 6 (Fig. 1e) =====\n")

print(t.test(data_m6$HDRS, data_m6$baseline_HDRS, paired = TRUE))

# Cohen's d (effectsize), computed on the two sets of scores
data_m6_long = data_m6 %>%
  pivot_longer(cols = c(HDRS, baseline_HDRS), names_to = "time", values_to = "score")
print(effectsize::cohens_d(score ~ time, data = data_m6_long))

cat(sprintf("\nHDRS-17 reduction at month 6: mean %.1f%% (min %.1f%%, max %.1f%%)\n",
            100 * mean(data_m6$deltaPerHDRS), 100 * min(data_m6$deltaPerHDRS), 100 * max(data_m6$deltaPerHDRS)))
cat(sprintf("Responders (>50%% reduction): %d of %d; remission (HDRS-17 <= 7): %d of %d\n",
            sum(data_m6$deltaPerHDRS > 0.5), nrow(data_m6), sum(data_m6$HDRS <= 7), nrow(data_m6)))

####################################################
# SECTION 2: STIMULATION EFFECTS ON NT             #
####################################################

# Fig. 2a-b, Supplementary Tables 3 (RL) and 4 (UG). The DA (RL) and 5-HT (UG) stimulation terms are the
# Family 1 primary tests (Section 6).
cat("\n===== Section 2: Stimulation models, reversal learning (Fig. 2a, Supplementary Table 3) =====\n")
rl.lm.DA = lmer(Oz ~ stim * trial_z + (1+stim|idx),
                data = task.merge[task.merge$nt == "DA" & task.merge$task == "RL",],
                REML=F, control = lmerControl(optimizer = "bobyqa"))
rl.lm.SE = lmer(Oz ~ stim * trial_z + (1+stim|idx),
                data = task.merge[task.merge$nt == "SE" & task.merge$task == "RL",],
                REML=F, control = lmerControl(optimizer = "bobyqa"))
rl.lm.NE = lmer(Oz ~ stim * trial_z + (1+stim|idx),
                data = task.merge[task.merge$nt == "NE" & task.merge$task == "RL",],
                REML=F, control = lmerControl(optimizer = "bobyqa"))
print(summary(rl.lm.DA))
print(summary(rl.lm.SE))
print(summary(rl.lm.NE))

cat("\n===== Section 2: Stimulation models, ultimatum game (Fig. 2b, Supplementary Table 4) =====\n")
ug.lm.DA = lmer(Oz ~ stim * trial_z + (1+stim|idx),
                data = task.merge[task.merge$nt == "DA" & task.merge$task == "UG",],
                REML=F, control = lmerControl(optimizer = "bobyqa"))
ug.lm.SE = lmer(Oz ~ stim * trial_z + (1+stim|idx),
                data = task.merge[task.merge$nt == "SE" & task.merge$task == "UG",],
                REML=F, control = lmerControl(optimizer = "bobyqa"))
ug.lm.NE = lmer(Oz ~ stim * trial_z + (1+stim|idx),
                data = task.merge[task.merge$nt == "NE" & task.merge$task == "UG",],
                REML=F, control = lmerControl(optimizer = "bobyqa"))
print(summary(ug.lm.DA))
print(summary(ug.lm.SE))
print(summary(ug.lm.NE))

cat("\n===== Section 2: Task x stimulation models (Supplementary Table 5) =====\n")
sess.lm.DA = lmer(Oz ~ stim * task * trial_z + (1+stim+task|idx),
                  data = task.merge[task.merge$nt == "DA",],
                  REML=F, control = lmerControl(optimizer = "bobyqa"))
sess.lm.SE = lmer(Oz ~ stim * task * trial_z + (1+stim+task|idx),
                  data = task.merge[task.merge$nt == "SE",],
                  REML=F, control = lmerControl(optimizer = "bobyqa"))
sess.lm.NE = lmer(Oz ~ stim * task * trial_z + (1+stim+task|idx),
                  data = task.merge[task.merge$nt == "NE",],
                  REML=F, control = lmerControl(optimizer = "bobyqa"))
print(summary(sess.lm.DA))
print(summary(sess.lm.SE))
print(summary(sess.lm.NE))

# Optional Word export of Supplementary Tables 3-5 (export_docx, top of script)
if (export_docx) {
  si_model_table = function(models) {
    map_dfr(models, ~ tidy(.x, effects = "fixed"), .id = "Neurotransmitter") %>%
      transmute(Neurotransmitter,
                Effect = recode(term, "(Intercept)" = "Intercept", "stim1" = "Stimulation Session (Pre vs. Post)",
                                "trial_z" = "Trial Number (z-scored)", "stim1:trial_z" = "Stim × Trial",
                                "taskUG" = "Task (RL vs. UG)", "stim1:taskUG" = "Stim × Task",
                                "taskUG:trial_z" = "Task × Trial", "stim1:taskUG:trial_z" = "Stim × Task × Trial"),
                Est. = sprintf("%.2f", estimate), SE = sprintf("%.2f", std.error), T = sprintf("%.2f", statistic),
                P = case_when(p.value < 0.001 ~ "< 0.001***",
                              p.value < 0.01 ~ paste0(sprintf("%.3f", p.value), "**"),
                              p.value < 0.05 ~ paste0(sprintf("%.3f", p.value), "*"),
                              TRUE ~ sprintf("%.3f", p.value)))
  }
  nts = c("Dopamine", "Serotonin", "Norepinephrine")
  write_si_table_docx(si_model_table(setNames(list(rl.lm.DA, rl.lm.SE, rl.lm.NE), nts)),
                      "Supplementary Table 3 (N-RL-M0)", "Supplementary_Table_03.docx")
  write_si_table_docx(si_model_table(setNames(list(ug.lm.DA, ug.lm.SE, ug.lm.NE), nts)),
                      "Supplementary Table 4 (N-UG-M0)", "Supplementary_Table_04.docx")
  write_si_table_docx(si_model_table(setNames(list(sess.lm.DA, sess.lm.SE, sess.lm.NE), nts)),
                      "Supplementary Table 5 (N-B-M0)", "Supplementary_Table_05.docx")
}

cat("\n===== Section 2: Block type (RL) and offer amount (UG) ANOVAs on participant means (Fig. 2 legend) =====\n")
rl.EST.aov = rl.EST.pr %>% group_by(idx,stim,nt,block_type) %>%
  summarise(oz = mean(Oz))
ug.EST.aov = ug.EST.pr %>% group_by(idx,stim,nt,offer_bin_f) %>%
  summarise(oz = mean(Oz))
for (n in c("DA","SE")) {
  cat("\nRL", n, "\n"); print(summary(aov(oz ~ stim * block_type, data = rl.EST.aov[rl.EST.aov$nt == n,])))
  cat("\nUG", n, "\n"); print(summary(aov(oz ~ stim * offer_bin_f, data = ug.EST.aov[ug.EST.aov$nt == n,])))
}

cat("\n===== Section 2: Norepinephrine, participant means pre vs. post (Extended Data Fig. 3a, c) =====\n")
for (task_est in list(list("RL", rl.EST.Reward), list("UG", ug.EST.Offer))) {
  ne = task_est[[2]] %>%
    filter(nt == "NE") %>%
    group_by(idx, stim) %>%
    summarise(m = mean(Oz), .groups = "drop") %>%
    pivot_wider(names_from = stim, values_from = m)
  cat("\n", task_est[[1]], "NE\n"); print(t.test(ne$`Post-Stim`, ne$`Pre-Stim`, paired = TRUE))
}

####################################################
# SECTION 3: LONGITUDINAL BEHAVIOUR (Fig. 3a-e)    #
####################################################
cat("\n===== Section 3: Longitudinal behaviour, paired t-tests vs. preoperative baseline (Fig. 3a-e) =====\n")

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

longitudinal_tests = pmap_dfr(
  list(list(rl.beh.means, rl.beh.means, ug.beh.means, ug.beh.means, ug.mood.means),
       c("mLogRT", "mChoice", "mLogRT", "mChoice", "mMood"),
       c("Fig. 3a RL log RT", "Fig. 3b RL optimal choice", "Fig. 3c UG log RT", "Fig. 3d UG acceptance", "Fig. 3e UG mood")),
  function(d, value_col, label) {
    map_dfr(c("Month 1", "Month 3", "Month 6"), function(m) {
      tt = calculate_t(d, value_col = value_col, time_col = "sess",
                       baseline_level = "Baseline", comparison_level = m)
      tibble(panel = label, comparison = paste("Baseline vs.", m),
             t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value)
    })
  })
print(longitudinal_tests, n = Inf)

####################################################
# SECTION 4: NT CHANGE vs. BEHAVIOURAL CHANGE      #
####################################################
cat("\n===== Section 4: Spearman correlations, NT change vs. behavioural change (Fig. 3f-g) =====\n")

# Fig. 3f: 5-HT change (RL) vs. change in log RT; Fig. 3g: DA change (UG) vs. change in mood
print(rl.beh.nt.long %>%
        filter(NT == "SE", beh == "rt") %>%
        group_by(NT, beh, month) %>%
        do(tidy(cor.test(.$value, .$ests, method = "spearman"))))
print(ug.mood.nt.long %>%
        filter(NT == "DA", beh == "mood") %>%
        group_by(NT, beh, month) %>%
        do(tidy(cor.test(.$value, .$ests, method = "spearman"))))

####################################################
# SECTION 5: NT CHANGE vs. HDRS-17                 #
####################################################
cat("\n===== Section 5: DA x 5-HT change predicting HDRS-17 (Results; Fig. 4a) =====\n")

model_RL_m1 <- lm(HDRS ~ deltaDA_RL * deltaSE_RL, data = data_m1)
model_RL_m3 <- lm(HDRS ~ deltaDA_RL * deltaSE_RL, data = data_m3)
model_RL_m6 <- lm(HDRS ~ deltaDA_RL * deltaSE_RL, data = data_m6)
model_UG_m1 <- lm(HDRS ~ deltaDA_UG * deltaSE_UG, data = data_m1)
model_UG_m3 <- lm(HDRS ~ deltaDA_UG * deltaSE_UG, data = data_m3)
model_UG_m6 <- lm(HDRS ~ deltaDA_UG * deltaSE_UG, data = data_m6)

for (m in c("model_RL_m1", "model_RL_m3", "model_RL_m6", "model_UG_m1", "model_UG_m3", "model_UG_m6")) {
  cat("\n", m, "\n"); print(summary(get(m)))
}

cat("\n===== Section 5: Model comparison, month 6 HDRS-17 (Supplementary Table 21) =====\n")

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

####################################################
# SECTION 6: MULTIPLICITY CORRECTION               #
# (Supplementary Table 27)                         #
####################################################
cat("\n===== Section 6: Multiplicity correction (Supplementary Table 27; Methods) =====\n")

# Family 1: DA increase in RL (test a) and 5-HT increase in UG (test b), from the Section 2 models.
# Family 2: DA x 5-HT interaction predicting month 6 HDRS-17 (test c, model_DAxSE).
test_a_single_task = tidy(rl.lm.DA) %>%
  filter(effect == "fixed", term == "stim1") %>%
  mutate(test = "a_DA_RL_single_task_model")

test_b_single_task = tidy(ug.lm.SE) %>%
  filter(effect == "fixed", term == "stim1") %>%
  mutate(test = "b_5HT_UG_single_task_model")

test_c = tidy(model_DAxSE) %>%
  filter(term == "deltaDA_UG:deltaSE_UG") %>%
  mutate(test = "c_DAx5HT_HDRS_key_interaction")

stopifnot("Family 1 tests not found in the model output" =
            nrow(test_a_single_task) == 1 && nrow(test_b_single_task) == 1 && nrow(test_c) == 1)

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

cat("\nFamily 1 (two tests: DA x RL, 5-HT x UG)\n")
correction_table_family1 = build_correction_table(
  primary_family_table %>% filter(test %in% c("a_DA_RL_single_task_model", "b_5HT_UG_single_task_model"))
)
print(correction_table_family1 %>% as.data.frame())

cat("\nThree-test grouping (Family 1 + Family 2)\n")
correction_table_3test = build_correction_table(primary_family_table)
print(correction_table_3test %>% as.data.frame())

cat("\n===== Section 7: Robustness of Fig. 3f/g correlations (Supplementary Tables 28-29) =====\n")

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

# Supplementary Table 28: bootstrap confidence intervals
print(fig3_boot_ci)

# Sanity check: at least 95% of bootstrap resamples were usable in every panel/month
stopifnot(all(fig3_boot_ci$n_valid >= 0.95 * fig3_boot_ci$B))

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

print(fig3_influence, n = Inf)

# Supplementary Table 29: most influential participant (rank-based Cook's D) per panel and month
fig3_influence_top = fig3_influence %>%
  group_by(panel, month) %>%
  mutate(n = n()) %>%
  filter(cooks_d_rank_rank == 1) %>%
  ungroup() %>%
  transmute(panel, month, n, subject, cooks_d_rank, loo_delta_rho,
            reverses_sign = sign(loo_rho) != sign(rho_obs))
print(fig3_influence_top)

cat("\n===== Section 8: Robustness of the DA x 5-HT -> HDRS-17 interaction (Supplementary Tables 30-31) =====\n")

# --- Supplementary Table 30, Month 6 column: bootstrap ---
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

# --- Supplementary Table 30, Month 6 column: permutation test ---
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

##############################################################
# SUPPLEMENTARY TABLE 31: month 6 DA x 5-HT interaction,     #
# leave-one-out refits and influence diagnostics             #
##############################################################

# Leave-one-out: refit the model with each participant removed in turn (interaction term)
loo_lm = function(data, model_formula) {
  map_dfr(seq_len(nrow(data)), function(i) {
    fit_i = lm(model_formula, data = data[-i, ])
    tidy(fit_i) %>% mutate(idx_removed = data$idx[i])
  })
}

model_DAxSE_loo = loo_lm(data_m6, HDRS ~ deltaDA_UG * deltaSE_UG) %>%
  filter(term == "deltaDA_UG:deltaSE_UG") %>%
  select(idx_removed, estimate, std.error, statistic, p.value) %>%
  arrange(desc(p.value))

# Influence diagnostics on the full-sample model (Cook's D, DFBETAS)
stopifnot(
  "cooks.distance()/dfbetas() row count no longer matches data_m6 -- likely a dropped observation from lm()'s na.action; do not proceed with the bind_cols() below until this is resolved" =
    nrow(data_m6) == length(cooks.distance(model_DAxSE))
)

model_DAxSE_influence = tibble(idx = data_m6$idx, cooks_d = cooks.distance(model_DAxSE)) %>%
  bind_cols(as_tibble(dfbetas(model_DAxSE)) %>% rename_with(~ paste0("dfbetas_", .x))) %>%
  mutate(cooks_d_rank = rank(-cooks_d, ties.method = "min")) %>%
  arrange(cooks_d_rank)

cat("\nSupplementary Table 31: month 6 interaction, leave-one-out refits with influence diagnostics\n")
supp_table_31 = model_DAxSE_loo %>%
  left_join(model_DAxSE_influence %>% select(idx, cooks_d, cooks_d_rank, `dfbetas_deltaDA_UG:deltaSE_UG`),
            by = c("idx_removed" = "idx")) %>%
  transmute(participant_removed = idx_removed,
            loo_estimate = estimate, loo_se = std.error, loo_t = statistic, loo_p = p.value,
            cooks_d, cooks_d_rank, dfbetas_interaction = `dfbetas_deltaDA_UG:deltaSE_UG`)
print(supp_table_31, n = Inf)


cat("\n===== Section 9: Robustness of the Fig. 4b product correlation (Supplementary Tables 32-33) =====\n")

# Fig. 4b plots month 6 HDRS-17 against the product of the UG DA and 5-HT changes
data_m6 = data_m6 %>%
  mutate(product_DAxSE = deltaDA_UG * deltaSE_UG)

##############################################################
# SUPPLEMENTARY TABLE 32: Fig. 4b correlation, full-sample   #
# test, bootstrap and permutation                            #
##############################################################

cat("\nSupplementary Table 32: Fig. 4b correlation\n")
# Full-sample Pearson test (r, t, df, P, analytic 95% CI)
fig4b_cortest = cor.test(data_m6$product_DAxSE, data_m6$HDRS, method = "pearson")
print(fig4b_cortest)

# Percentile bootstrap (2,000 case resamples): 95% CI, median, IQR, % same sign
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

# Permutation test (10,000 shuffles of HDRS-17 across participants)
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

##############################################################
# SUPPLEMENTARY TABLE 33: Fig. 4b correlation, leave-one-    #
# out refits and influence diagnostics                       #
##############################################################

# Leave-one-out: Pearson r recomputed with each participant removed in turn
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

# Influence diagnostics: the correlation fitted as a simple regression (HDRS ~ product), whose slope
# test is identical to the Pearson test; Cook's D and DFBETAS for the product term
fig4b_lm = lm(HDRS ~ product_DAxSE, data = data_m6)

stopifnot(
  "cooks.distance()/dfbetas() row count no longer matches data_m6 -- likely a dropped observation from lm()'s na.action; do not proceed with the bind_cols() below until this is resolved" =
    nrow(data_m6) == length(cooks.distance(fig4b_lm))
)

fig4b_influence = tibble(idx = data_m6$idx, cooks_d = cooks.distance(fig4b_lm)) %>%
  bind_cols(as_tibble(dfbetas(fig4b_lm)) %>% rename_with(~ paste0("dfbetas_", .x))) %>%
  mutate(cooks_d_rank = rank(-cooks_d, ties.method = "min")) %>%
  arrange(cooks_d_rank)

cat("\nSupplementary Table 33: Fig. 4b correlation, leave-one-out refits with influence diagnostics\n")
supp_table_33 = fig4b_loo %>%
  left_join(fig4b_influence %>% select(idx, cooks_d, cooks_d_rank, dfbetas_product_DAxSE),
            by = c("idx_removed" = "idx")) %>%
  transmute(participant_removed = idx_removed, loo_r, loo_p,
            cooks_d, cooks_d_rank, dfbetas_product = dfbetas_product_DAxSE)
print(supp_table_33, n = Inf)
