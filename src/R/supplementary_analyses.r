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
# 2025/09/16    Blair Shevlin                          adapted "statistical_tests.R" to use alternative z-scoring
# 2025/11/14    Blair Shevlin                          revised ANOVA table generation to include all metrics
# 2026/02/23    Blair Shevlin                          exploring behavioral analyses and clinical correlates of NT change

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
cl.df.m6 = cl.df %>% filter(session == "month 6") %>% select(idx,HDRS)

# Load behavioral data
rl.beh = read.csv(file = beh_dir / "rl-data_deid_07-10-25.csv")

ug.beh = read.csv(file = beh_dir / "ug-data_deid_07-10-25.csv")

############################
# Alternative NT Analyses #
###########################
rl.EST.Reward = rl.EST %>%
  filter(event == "Reward") %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))

ug.EST.Offer = ug.EST %>%
  filter(event == "Offer") %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim"))) 

# Regressions

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


# Trial-level analysis

# UG
ug.DA.lm = lmer(   data = ug.EST.pr[ug.EST.pr$nt == "DA",],
                formula = Oz ~ offer_z*stim + trial_z*stim + (1+stim|idx),
                REML=F,
                control = lmerControl(optimizer="nmkbw") )
summary(ug.DA.lm) # trial, trial x stim

ug.SE.lm = lmer(   data = ug.EST.pr[ug.EST.pr$nt == "SE",],
                formula = Oz ~ offer_z*stim + trial_z*stim + (1+stim|idx),
                REML=F,
                control = lmerControl(optimizer="nmkbw") )
summary(ug.SE.lm) # stim

ug.NE.lm = lmer(data = ug.EST.pr[ug.EST.pr$nt == "NE",],
                formula = Oz ~ offer_z*stim + trial_z*stim + (1+stim|idx),
                REML=F,
                control = lmerControl(optimizer="nmkbw") )
summary(ug.NE.lm) # trial, trial x stim


# RL
rl.DA.lm = lmer(data = rl.EST.pr[rl.EST.pr$nt == "DA",],
                formula = Oz ~ rew*stim + trial_z*stim + (1+stim|idx),
                REML=F,
                control = lmerControl(optimizer="nmkbw"))
summary(rl.DA.lm) # Stim, trial, trial x stim

rl.SE.lm = lmer(   data = rl.EST.pr[rl.EST.pr$nt == "SE",],
                formula = Oz ~ rew*stim + trial_z*stim + (1+stim|idx),
                REML=F,
                control = lmerControl(optimizer="nmkbw"))
summary(rl.SE.lm) # Reward, trial, reward x post-stim

rl.NE.lm = lmer(data = rl.EST.pr[rl.EST.pr$nt == "NE",],
                formula = Oz ~ rew*stim + trial_z*stim + (1+stim|idx),
                REML=F,
                control = lmerControl(optimizer="nmkbw"))
summary(rl.NE.lm) # trial, trial x stim


# Now combine for task-comparison

rl.premerge = rl.EST.pr %>%
  select(idx,stim,nt,Oz,trial_z) %>%
  mutate(task = "RL")
ug.premerge = ug.EST.pr %>%
  select(idx,stim,nt,Oz,trial_z) %>%
  mutate(task = "UG")
task.merge = rbind(rl.premerge,ug.premerge)

# Table of Subject-level NT values
task.merge %>% group_by(idx,stim,task, nt) %>%
  summarise(m = mean(Oz), s = sd(Oz)) %>% as.data.frame()

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
  as.data.frame()

# Us
DA.lm = lmer(   data = task.merge[task.merge$nt == "DA",],
                formula = Oz ~ stim * task * trial_z + (1+stim+task|idx),
                control = lmerControl(optimizer="nmkbw"))
summary(DA.lm) # stim x task (and trial)
SE.lm = lmer(   data = task.merge[task.merge$nt == "SE",],
                formula = Oz ~ stim * task * trial_z + (1+stim+task|idx),
                control = lmerControl(optimizer="nmkbw"))
summary(SE.lm) # stim x task (not trial)
NE.lm = lmer(   data = task.merge[task.merge$nt == "NE",],
                formula = Oz ~ stim * task * trial_z + (1+stim+task|idx),
                control = lmerControl(optimizer="nmkbw"))
summary(NE.lm) # stim x task (and trial)


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

anova_summaries %>%
  filter(task == "RL", nt == "DA", term == "stim")
# Significant for all metrics

anova_summaries %>%
  filter(task == "RL", nt == "SE", term == "stim")
# None significant

anova_summaries %>%
  filter(task == "RL", nt == "NE", term == "stim")
# None significant

anova_summaries %>%
  filter(task == "UG", nt == "DA", term == "stim")
# All significant, except Rz, Rz10

anova_summaries %>%
  filter(task == "UG", nt == "SE", term == "stim")
# All significant, except Rz, Rz10

anova_summaries %>%
  filter(task == "UG", nt == "NE", term == "stim")
# None significant


# Look at RL with rew instead of block 

rl.EST.aov.alt = rl.EST.pr %>% group_by(idx,stim,nt,rew_f) %>%
  summarise(oz = mean(Oz), mz = mean(Mz), pz = mean(Pz), rz = mean(Rz), 
            oz10 = mean(Oz10), mz10 = mean(Mz10), pz10 = mean(Pz10), rz10 = mean(Rz10))

rl_anovas_alt <- map2(
  rl_combinations$metric,
  rl_combinations$nt_type,
  ~ run_anova(rl.EST.aov.alt, .x, .y, "~ stim * rew_f")
)

rl_summaries_alt <- pmap_dfr(
  list(
    rl_anovas_alt,
    rep("RL", length(rl_anovas_alt)),
    rep(rl_combinations$nt_type, 1),
    rep(rl_combinations$metric, 1)
  ),
  extract_anova_summary
)

rl_summaries_alt %>%
  filter(task == "RL", nt == "DA", term == "stim")
# Significant for all metrics

rl_summaries_alt %>%
  filter(task == "RL", nt == "SE", term == "stim")
# None significant

rl_summaries_alt %>%
  filter(task == "RL", nt == "NE", term == "stim")
# None significant


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

# Create flextable for publication
create_anova_flextable <- function(formatted_data) {
  
  # Add residual df (should be consistent within task)
  # Assuming df2 = 42 for RL and similar for UG based on your structure
  
  ft <- formatted_data %>%
    flextable() %>%
    theme_vanilla() %>%
    
    # Header formatting
    bold(part = "header") %>%
    bg(bg = "#f0f0f0", part = "header") %>%
    align(align = "center", part = "header") %>%
    
    # Body alignment
    align(j = 1:3, align = "left", part = "body") %>%
    align(j = 4:5, align = "center", part = "body") %>%
    
    # Fonts
    fontsize(size = 10) %>%
    font(fontname = "Arial") %>%
    
    # Column widths
    width(j = 1, width = 1.8) %>%  # Task
    width(j = 2, width = 1.5) %>%  # Neurotransmitter
    width(j = 3, width = 1.3) %>%  # Metric
    width(j = 4, width = 0.8) %>%  # F-statistic
    width(j = 5, width = 1.0) %>%  # p-value
    
    # Padding
    padding(padding.top = 3, padding.bottom = 3) %>%
    valign(valign = "center") %>%
    
    # Merge cells for repeated task/neurotransmitter values
    merge_v(j = c("Task", "Neurotransmitter")) %>%
    
    # Add footnote
    add_footer_lines(
      values = c(
        "Note: Two-way repeated-measures ANOVAs testing stimulation effect (Pre vs. Post).",
        "RL task includes block type (Reward/Mixed/Punishment); UG task includes offer amount bins.",
        "***p < .001, **p < .01, *p < .05.",
        "df1 = 1 for all stimulation effects. Residual df: RL = 42, UG = [check your data]."
      )
    ) %>%
    fontsize(size = 8, part = "footer") %>%
    italic(part = "footer")
  
  return(ft)
}

# Alternative: Wide format table (metrics as columns)
create_wide_format_table <- function(formatted_data) {
  
  wide_data <- formatted_data %>%
    select(Task, Neurotransmitter, Metric, p) %>%
    pivot_wider(
      names_from = Metric,
      values_from = p,
      names_sort = TRUE
    )
  
  ft <- wide_data %>%
    flextable() %>%
    theme_vanilla() %>%
    bold(part = "header") %>%
    bg(bg = "#f0f0f0", part = "header") %>%
    align(align = "center", part = "header") %>%
    align(j = 1:2, align = "left", part = "body") %>%
    align(j = 3:ncol(wide_data), align = "center", part = "body") %>%
    fontsize(size = 9) %>%
    font(fontname = "Times New Roman") %>%
    width(j = 1, width = 1.5) %>%
    width(j = 2, width = 1.3) %>%
    autofit() %>%
    merge_v(j = "Task") %>%
    add_footer_lines(
      values = c(
        "Note: p-values from two-way repeated-measures ANOVAs testing stimulation main effect.",
        "***p < .001, **p < .01, *p < .05"
      )
    ) %>%
    fontsize(size = 8, part = "footer") %>%
    italic(part = "footer")
  
  return(ft)
}

# Main function to generate all tables
generate_anova_summary_tables <- function(anova_summaries, output_dir = "results") {
  
  # Format the data
  formatted_data <- format_anova_results(anova_summaries)
  
  # Create both table versions
  long_table <- create_anova_flextable(formatted_data)
  wide_table <- create_wide_format_table(formatted_data)
  
  # Save long format (recommended for publication)
  long_doc <- read_docx() %>%
    body_add_par("Supplementary Table 13", style = "heading 1") %>%
    body_add_par("Robustness Analysis of Neurotransmitter Quantification Metrics", 
                 style = "heading 2") %>%
    body_add_flextable(long_table)
  
  print(long_doc, target = file.path(output_dir, "Table_13_ANOVA_metrics_long.docx"))
  
  cat("Tables generated successfully!\n")
  cat(sprintf("Long format: %s\n", file.path(output_dir, "Table_13_ANOVA_metrics_long.docx")))
  return(list(
    formatted_data = formatted_data,
    long_table = long_table
  ))
}

anova_res = generate_anova_summary_tables(anova_summaries)


#########################
# Alternative Z-Scoring #
#########################

# Generate summary data - focusing only on reward events (RL) and offer events (UG)
# This time using data that was z-scored using both tasks
rl.EST.Shared= 
  tasks.EST.alt %>% filter(task == "RL") %>%
  filter(event == "Shared")

ug.EST.Shared = 
  tasks.EST.alt %>% filter(task == "UG") %>%
  filter(event == "Shared")

rl.EST.Choice = 
  tasks.EST.alt %>% filter(task == "RL") %>%
  filter(event == "Choice") 

ug.EST.Choice = 
  tasks.EST.alt %>% filter(task == "UG") %>%
  filter(event == "Choice")


# Merge with task
rl.task.OR = 
    rl.beh %>% 
    filter(sess %in% c("Pre-Stim","Post-Stim")) %>%
    mutate(stim = factor(sess, levels = c("Pre-Stim","Post-Stim"))) %>% 
    select(!c(session,sess,score))

ug.task.OR = 
    ug.beh %>% 
    filter(sess %in% c("Pre-Stim","Post-Stim")) %>%
    mutate(stim = factor(sess, levels = c("Pre-Stim","Post-Stim"))) %>% 
    select(!c(session,sess))

rl.EST.Shared = merge(rl.EST.Shared,rl.task.OR, by = c("idx","stim","trial"))
ug.EST.Shared = merge(ug.EST.Shared,ug.task.OR)

rl.EST.Choice = merge(rl.EST.Choice,rl.task.OR, by = c("idx","stim","trial"))
ug.EST.Choice = merge(ug.EST.Choice,ug.task.OR)

# Regression
rl.EST.pr = 
    rl.EST.Shared %>%
    group_by(idx,stim,block,event,nt) %>%
    mutate(prev_rew_raw = ifelse(trial_within_block==1,NA,lag(outcome)),
                        prev_rew = ifelse(is.na(prev_rew_raw),NA,
                                        ifelse(prev_rew_raw == 1,
                                                ifelse(cond %in% c("Mixed", "Positive"),10,0),
                                                ifelse(cond %in% c("Mixed", "Negative"),-10,0)
                                        )
                        )
                ) %>% 
  ungroup() %>%
  mutate(nt = factor(nt, levels = c("SE","DA","NE")),
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
         rew_f = factor(rew, levels=c(0,-10,10)))

ug.EST.pr = 
    ug.EST.Shared %>%
    group_by(idx,stim,event,nt) %>%
      mutate(prev_offer_raw = ifelse(trial==1,NA,lag(offer)),
             prev_offer = ifelse(is.na(prev_offer_raw),NA,
                                 ifelse(prev_offer_raw > offer, "-PE","+PE"))
      ) %>%
    ungroup() %>%
  mutate(nt = factor(nt, levels = c("SE","DA","NE")),
         offer_bin = ifelse(offer <4,"Low",
                            ifelse(offer>6,"High","Middle")),
         offer_bin = factor(offer_bin,
                            levels = c("Low","Middle","High")),
         offer_bin_f = factor(offer_bin, levels = c("Middle","Low","High")),
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
rl.EST.aov = rl.EST.pr %>% 
    group_by(idx,stim,nt,block_type) %>%
  summarise(oz = mean(Mz))
ug.EST.aov = ug.EST.pr %>% group_by(idx,stim,nt,offer_bin_f) %>%
  summarise(oz = mean(Mz))

rl.DA.anova=aov(oz ~ stim * block_type,
                data = rl.EST.aov[rl.EST.aov$nt == "DA",])
rl.SE.anova=aov(oz ~ stim * block_type,
                data = rl.EST.aov[rl.EST.aov$nt == "SE",])

ug.DA.anova=aov(oz ~ stim * offer_bin_f,
                data = ug.EST.aov[ug.EST.aov$nt == "DA",])
ug.SE.anova=aov(oz ~ stim * offer_bin_f,
                data = ug.EST.aov[ug.EST.aov$nt == "SE",])

summary(rl.DA.anova) # Stim effect still significant
summary(rl.SE.anova) # Stim effect still non-significant
summary(ug.DA.anova) # Stim effect NOW non-significant
summary(ug.SE.anova) # Stim effect still significant


# Look at combined task

# Trial-level analysis
DA.lm = lmer(   data = tasks.EST.alt[tasks.EST.alt$event == "Shared" & tasks.EST.alt$nt == "DA",],
                formula = Oz ~ stim * task + (1|idx),
                control = lmerControl(optimizer="bobyqa"))
summary(DA.lm)
SE.lm = lmer(   data = tasks.EST.alt[tasks.EST.alt$event == "Shared" & tasks.EST.alt$nt == "SE",],
                formula = Oz ~ stim * task + (1|idx),
                control = lmerControl(optimizer="bobyqa"))
summary(SE.lm)
NE.lm = lmer(   data = tasks.EST.alt[tasks.EST.alt$event == "Shared" & tasks.EST.alt$nt == "NE",],
                formula = Oz ~ stim * task + (1|idx),
                control = lmerControl(optimizer="bobyqa"))
summary(NE.lm)

# Session/Task-level
tasks.anova = tasks.EST.alt %>% 
    filter(event == "Shared") %>%
    group_by(idx,task, stim, nt) %>%
    summarise(oz = mean(Oz))

DA.anova = aov(oz ~ stim * task,
                data = tasks.anova[tasks.anova$nt == "DA",])
                
SE.anova = aov(oz ~ stim * task,
                data = tasks.anova[tasks.anova$nt == "SE",])

summary(DA.anova)
summary(SE.anova)

TukeyHSD(DA.anova)
TukeyHSD(SE.anova)

################
#   Behavior   #
################

# Trial-level changes
ug.beh.OR = ug.beh %>%
    filter(sess %in% c("Pre-Stim","Post-Stim"), rt < 10) %>%
    mutate( offer_z = scale(offer)[,1],
            logRT = log(rt),
            sess = factor(sess, levels = c("Pre-Stim","Post-Stim")),
            accept = 1* (rej==0))


rl.beh %>% head()

rl.beh.OR = rl.beh %>%
    filter(sess %in% c("Pre-Stim","Post-Stim"), rt < 10) %>%
    mutate( cond = factor(cond, levels = c("Mixed","Negative","Positive")),
            logRT = log(rt),
            sess = factor(sess, levels = c("Pre-Stim","Post-Stim")),
            otpimal = opt,
            rew = ifelse(outcome==1,
                         ifelse(cond=="Negative",0,10),
                         ifelse(cond=="Positive",0,-10)),
            rew_f = factor(rew, levels=c(0,-10,10))
    )

ug.mChoice = glmer( accept ~ offer_z + sess + (1+offer_z | idx),
                    data = ug.beh.OR,
                    family = binomial,
                    control = glmerControl(optimizer='bobyqa'))
ug.mRT = lmer( logRT ~ offer_z + sess + (1+offer_z  | idx),
                    data = ug.beh.OR,
                    control = lmerControl(optimizer='bobyqa'))
ug.mMood = lmer( mood ~ offer_z + sess + (1+offer_z  | idx),
                    data = ug.beh.OR[!is.na(ug.beh.OR$mood), ],
                    control = lmerControl(optimizer='bobyqa'))

rl.mChoice = glmer( otpimal ~ cond + sess + (1+cond  | idx),
                    data = rl.beh.OR,
                    family = binomial,
                    control = glmerControl(optimizer='bobyqa'))
rl.mRT = lmer( logRT ~ cond + sess + (1+cond  | idx),
                    data = rl.beh.OR,
                    control = lmerControl(optimizer='bobyqa'))

summary(ug.mChoice) # no stim effect
summary(ug.mRT) # no stim effect
summary(ug.mMood) # stim effect

summary(rl.mChoice) # no stim effect
summary(rl.mRT) # Stim effects on RT

# Look at specific condition for RL RT effects
rl.mRT.M = lmer( logRT ~ sess + (1 | idx),
                    data = rl.beh.OR[rl.beh.OR$cond == "Mixed",],
                    control = lmerControl(optimizer='bobyqa'))
rl.mRT.N = lmer( logRT ~ sess + (1 | idx),
                    data = rl.beh.OR[rl.beh.OR$cond == "Negative",],
                    control = lmerControl(optimizer='bobyqa'))
rl.mRT.P = lmer( logRT ~ sess + (1 | idx),
                    data = rl.beh.OR[rl.beh.OR$cond == "Positive",],
                    control = lmerControl(optimizer='bobyqa'))

summary(rl.mRT.M) # Faster post-stim
summary(rl.mRT.N) # No stim effect
summary(rl.mRT.P) # Faster post-stim

#######################
# Clinical correlates #
#######################

# Could NT changes be explained by baseline HDRS?

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


rl.reg.DA.HDRS = lmer( Oz ~ stim * HDRS_z + stim * trial_z + (1+stim|idx),
                    data= rl.EST.pr[rl.EST.pr$nt == "DA",],
                    control = lmerControl(optimizer='nmkbw'))
rl.reg.SE.HDRS = lmer( Oz ~ stim * HDRS_z + stim * trial_z + (1+stim|idx),
                    data= rl.EST.pr[rl.EST.pr$nt == "SE",],
                    control = lmerControl(optimizer='nmkbw'))
rl.reg.NE.HDRS = lmer( Oz ~ stim * HDRS_z + stim * trial_z + (1+stim|idx),
                    data= rl.EST.pr[rl.EST.pr$nt == "NE",],
                    control = lmerControl(optimizer='nmkbw'))

ug.reg.DA.HDRS = lmer( Oz ~ stim * HDRS_z + stim * trial_z + (1+stim|idx),
                    data= ug.EST.pr[ug.EST.pr$nt == "DA",],
                    control = lmerControl(optimizer='nmkbw'))
ug.reg.SE.HDRS = lmer( Oz ~ stim * HDRS_z + stim * trial_z + (1+stim|idx),
                    data= ug.EST.pr[ug.EST.pr$nt == "SE",],
                    control = lmerControl(optimizer='nmkbw'))
ug.reg.NE.HDRS = lmer( Oz ~ stim * HDRS_z + stim * trial_z + (1+stim|idx),
                    data= ug.EST.pr[ug.EST.pr$nt == "NE",],
                    control = lmerControl(optimizer='nmkbw'))

summary(rl.reg.DA.HDRS) # Stim only
summary(rl.reg.SE.HDRS) # Stim, HDRS, Stim x HDRS
summary(rl.reg.NE.HDRS) # Stim, HDRS, Stim x HDRS

summary(ug.reg.DA.HDRS) # Stim, Stim x HDRS
summary(ug.reg.SE.HDRS) # Stim, Stim x HDRS
summary(ug.reg.NE.HDRS) # Stim, Stim x HDRS

# Reduce to subject-levels
rl.EST.aov.hdrs = rl.EST.pr %>% group_by(idx,stim,nt,HDRS_z,HDRS) %>%
  summarise(oz = mean(Oz))
ug.EST.aov.hdrs = ug.EST.pr %>% group_by(idx,stim,nt,HDRS_z,HDRS) %>%
  summarise(oz = mean(Oz))

rl.simple.DA.HDRS = lm(oz ~ stim * HDRS_z, data = rl.EST.aov.hdrs[rl.EST.aov.hdrs$nt == "DA",])
ug.simple.DA.HDRS = lm(oz ~ stim * HDRS_z, data = ug.EST.aov.hdrs[ug.EST.aov.hdrs$nt == "DA",])

rl.simple.SE.HDRS = lm(oz ~ stim * HDRS_z, data = rl.EST.aov.hdrs[rl.EST.aov.hdrs$nt == "SE",])
ug.simple.SE.HDRS = lm(oz ~ stim * HDRS_z, data = ug.EST.aov.hdrs[ug.EST.aov.hdrs$nt == "SE",])

rl.simple.NE.HDRS = lm(oz ~ stim * HDRS_z, data = rl.EST.aov.hdrs[rl.EST.aov.hdrs$nt == "NE",])
ug.simple.NE.HDRS = lm(oz ~ stim * HDRS_z, data = ug.EST.aov.hdrs[ug.EST.aov.hdrs$nt == "NE",])

summary(rl.simple.DA.HDRS) # Stim only
summary(ug.simple.DA.HDRS) # Stim only

summary(rl.simple.SE.HDRS) # NS
summary(ug.simple.SE.HDRS) # Stim, HDRS x Stim

summary(rl.simple.NE.HDRS) # NS
summary(ug.simple.NE.HDRS) # HDRS x Stim


# Create individual plots for each neurotransmitter
RL_SE <- rl.EST.aov.hdrs %>%
  filter(nt == "SE") %>%
  ggplot(aes(x = HDRS, y = oz, color = stim)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, size = 2) +
  labs(title = "Reversal learning",
       x = "HDRS",color = element_blank(),
       y = "5-HT [z]") +
  theme_pubr(base_size = 14) +
  scale_color_manual(values = c("Pre-Stim" = "#9ecae1", "Post-Stim" = "#08519c"))  +
  coord_cartesian(ylim = c(-6,6))

RL_DA <- rl.EST.aov.hdrs %>%
  filter(nt == "DA") %>%
  ggplot(aes(x = HDRS, y = oz, color = stim)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, size = 2) +
  labs(#title = "Dopamine (DA)",
       x = "HDRS", color = element_blank(),
       y = "DA [z]") +
  theme_pubr(base_size = 14) +
  scale_color_manual(values = c("Pre-Stim" = "#fcbba1", "Post-Stim" = "#a50f15"))  +
  coord_cartesian(ylim = c(-6,6))

RL_NE <- rl.EST.aov.hdrs %>%
  filter(nt == "NE") %>%
  ggplot(aes(x = HDRS, y = oz, color = stim)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, size = 2) +
  labs(#title = "Norepinephrine (NE)",
       x = "HDRS",color = element_blank(),
       y = "NE [z]") +
  theme_pubr(base_size = 14) +
  scale_color_manual(values = c("Pre-Stim" = "#c7e9c0", "Post-Stim" = "#006d2c"))  +
  coord_cartesian(ylim = c(-6,6))

UG_SE <- ug.EST.aov.hdrs %>%
  filter(nt == "SE") %>%
  ggplot(aes(x = HDRS, y = oz, color = stim)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, size = 2) +
  labs(title = "Ultimatum game",
       x = "HDRS",color = element_blank(),
       y = "5-HT [z]") +
  theme_pubr(base_size = 14) +
  scale_color_manual(values = c("Pre-Stim" = "#9ecae1", "Post-Stim" = "#08519c"))  +
  coord_cartesian(ylim = c(-6,6))

UG_DA <- ug.EST.aov.hdrs %>%
  filter(nt == "DA") %>%
  ggplot(aes(x = HDRS, y = oz, color = stim)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, size = 2) +
  labs(#title = "Ultimatum game",
       x = "HDRS", color = element_blank(),
       y = "DA [z]") +
  theme_pubr(base_size = 14) +
  scale_color_manual(values = c("Pre-Stim" = "#fcbba1", "Post-Stim" = "#a50f15")) +
  coord_cartesian(ylim = c(-6,6))

UG_NE <- ug.EST.aov.hdrs %>%
  filter(nt == "NE") %>%
  ggplot(aes(x = HDRS, y = oz, color = stim)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, size = 2) +
  labs(#title = "Norepinephrine (NE)",
       x = "HDRS",color = element_blank(),
       y = "NE [z]") +
  theme_pubr(base_size = 14) +
  scale_color_manual(values = c("Pre-Stim" = "#c7e9c0", "Post-Stim" = "#006d2c")) +
  coord_cartesian(ylim = c(-6,6))

# Combine plots using patchwork
combined_HDRS_NT <- (RL_DA + RL_SE + RL_NE) / (UG_DA + UG_SE + UG_NE) 

# Display the plot
combined_HDRS_NT

ggsave(res_dir / "Extended-Data_Figure4_9-17-25.png", 
       plot = combined_HDRS_NT,
       device = "png",
       width = 10,          # Width in inches
       height = 6,         # Height in inches  
       dpi = 300)           # Resolution (300 DPI for publication quality)


# Potentially this is driven by HDRS outlier

rl.simple.DA.HDRS.noOutlier = lm(oz ~ stim * HDRS_z, data = rl.EST.aov.hdrs[rl.EST.aov.hdrs$nt == "DA" & rl.EST.aov.hdrs$idx !=6 ,])
ug.simple.DA.HDRS.noOutlier = lm(oz ~ stim * HDRS_z, data = ug.EST.aov.hdrs[ug.EST.aov.hdrs$nt == "DA" & ug.EST.aov.hdrs$idx !=6,])

rl.simple.SE.HDRS.noOutlier = lm(oz ~ stim * HDRS_z, data = rl.EST.aov.hdrs[rl.EST.aov.hdrs$nt == "SE" & rl.EST.aov.hdrs$idx !=6,])
ug.simple.SE.HDRS.noOutlier = lm(oz ~ stim * HDRS_z, data = ug.EST.aov.hdrs[ug.EST.aov.hdrs$nt == "SE" & ug.EST.aov.hdrs$idx !=6,])

rl.simple.NE.HDRS.noOutlier = lm(oz ~ stim * HDRS_z, data = rl.EST.aov.hdrs[rl.EST.aov.hdrs$nt == "NE" & rl.EST.aov.hdrs$idx !=6,])
ug.simple.NE.HDRS.noOutlier = lm(oz ~ stim * HDRS_z, data = ug.EST.aov.hdrs[ug.EST.aov.hdrs$nt == "NE" & ug.EST.aov.hdrs$idx !=6,])

summary(rl.simple.DA.HDRS.noOutlier) # NS
summary(ug.simple.DA.HDRS.noOutlier) # NS

summary(rl.simple.SE.HDRS.noOutlier) # NS
summary(ug.simple.SE.HDRS.noOutlier) # Stim, trending positive HDRS x Stim

summary(rl.simple.NE.HDRS.noOutlier) # NS
summary(ug.simple.NE.HDRS.noOutlier) # NS

# Robust regression suggested by Xiaosi

rl.EST.aov.hdrs = rl.EST.pr %>% group_by(idx,stim,nt,HDRS_z,HDRS) %>%
  summarise(oz = mean(Oz)) %>%
  group_by(idx,nt,HDRS_z,HDRS) %>%
  summarise(oz_diff = oz[stim=="Post-Stim"] - oz[stim=="Pre-Stim"])
ug.EST.aov.hdrs = ug.EST.pr %>% group_by(idx,stim,nt,HDRS_z,HDRS) %>%
  summarise(oz = mean(Oz)) %>%
  group_by(idx,nt,HDRS_z,HDRS) %>%
  summarise(oz_diff = oz[stim=="Post-Stim"] - oz[stim=="Pre-Stim"])

rl.robust.DA.HDRS.noOutlier = lm(oz_diff ~ HDRS_z, data = rl.EST.aov.hdrs[rl.EST.aov.hdrs$nt == "DA",])
ug.robust.DA.HDRS.noOutlier = lm(oz_diff ~ HDRS_z, data = ug.EST.aov.hdrs[ug.EST.aov.hdrs$nt == "DA",])

rl.robust.SE.HDRS.noOutlier = lm(oz_diff ~ HDRS_z, data = rl.EST.aov.hdrs[rl.EST.aov.hdrs$nt == "SE" ,])
ug.robust.SE.HDRS.noOutlier = lm(oz_diff ~ HDRS_z, data = ug.EST.aov.hdrs[ug.EST.aov.hdrs$nt == "SE" ,])

rl.robust.NE.HDRS.noOutlier = lm(oz_diff ~ HDRS_z, data = rl.EST.aov.hdrs[rl.EST.aov.hdrs$nt == "NE" ,])
ug.robust.NE.HDRS.noOutlier = lm(oz_diff ~ HDRS_z, data = ug.EST.aov.hdrs[ug.EST.aov.hdrs$nt == "NE" ,])

summary(rl.robust.DA.HDRS.noOutlier) # NS
summary(ug.robust.DA.HDRS.noOutlier) # NS

summary(rl.robust.SE.HDRS.noOutlier) # NS
summary(ug.robust.SE.HDRS.noOutlier) # NS

summary(rl.robust.NE.HDRS.noOutlier) # NS
summary(ug.robust.NE.HDRS.noOutlier) # NS

##############
# HDRS + Beh #
##############
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

# BASELINE CORRELATIONS

# UG #
######

# Pre-Stim (baseline)
cor.test(ug.cl.baseline.mean$HDRS_z, 
                   ug.cl.baseline.mean$mAccept, 
                   method = "pearson")
# NS
cor.test(ug.cl.baseline.mean$HDRS_z, 
                   ug.cl.baseline.mean$mRT, 
                   pmethod = "pearson")
# NS
cor.test(ug.cl.baseline.mean$HDRS_z, 
                   ug.cl.baseline.mean$mMood, 
                   method = "pearson")
# No difference in baseline mood

# RL #
######
cor.test(rl.cl.baseline.mean$HDRS_z, 
                   rl.cl.baseline.mean$mOpt, 
                   method = "pearson")
# Sig
cor.test(rl.cl.baseline.mean$HDRS_z, 
                   rl.cl.baseline.mean$mRT, 
                   method = "pearson")

# REGRESSION (trial-level)
ug.cl.baseline.pr = ug.cl.baseline %>%
  mutate(accept = 1 * (rej==0), logRT = log(rt))
rl.cl.baseline.pr = rl.cl.baseline %>%
  mutate(logRT = log(rt))

glmer.ug.mChoice = glmer(accept  ~ HDRS_z + offer_z + (1|idx), data = ug.cl.baseline.pr, 
                          family = binomial,
                          control = glmerControl(optimizer='bobyqa') )
lmer.ug.mRT = lmer(logRT ~ HDRS_z + offer_z + (1|idx), data = ug.cl.baseline.pr,
                          control = lmerControl(optimizer='bobyqa'))
lmer.ug.mMood = lmer(mood ~ HDRS_z + offer_z + (1|idx), data = ug.cl.baseline.pr,
                          control = lmerControl(optimizer='bobyqa'))

glmer.rl.mChoice = glmer(opt ~ HDRS_z + cond + (1|idx), data = rl.cl.baseline.pr,
                      family = binomial,
                      control = glmerControl(optimizer='bobyqa') )
lmer.rl.mRT = lmer(logRT ~ HDRS_z + cond + (1|idx), data = rl.cl.baseline.pr,
                          control = lmerControl(optimizer='bobyqa'))

summary(glmer.ug.mChoice) # NS
summary(lmer.ug.mRT) # NS
summary(lmer.ug.mMood) # NS

summary(glmer.rl.mChoice) # Significant
summary(lmer.rl.mRT) # NS                   

# OR #
rl.task.OR = 
    rl.beh %>% 
    filter(sess %in% c("Pre-Stim","Post-Stim"), rt < 10) %>%
    mutate(stim = factor(sess, levels = c("Pre-Stim","Post-Stim"))) %>% 
    select(!c(session,sess,score))

ug.task.OR = 
    ug.beh %>% 
    filter(sess %in% c("Pre-Stim","Post-Stim"), rt < 10) %>%
    mutate(stim = factor(sess, levels = c("Pre-Stim","Post-Stim"))) %>% 
    select(!c(session,sess)) 

ug.cl.OR = merge(ug.task.OR, cl.df.baseline) %>% mutate(offer_z = scale(offer)[,1])
rl.cl.OR = merge(rl.task.OR, cl.df.baseline)

ug.cl.OR.mean = ug.cl.OR %>% 
  group_by(idx,stim, HDRS,HDRS_z) %>%
  summarise(mAccept = mean(rej==0),
            mRT = mean(log(rt)),
            mMood = mean(mood, na.rm = T))

rl.cl.OR.mean = rl.cl.OR %>% 
  group_by(idx,stim, HDRS,HDRS_z) %>%
  summarise(mOpt = mean(opt==1),
            mRT = mean(log(rt)))


# OR CORRELATIONS

# UG #
######

# Pre-Stim (baseline)
cor.test(ug.cl.OR.mean$HDRS_z[ug.cl.OR.mean$stim == "Pre-Stim"], 
                   ug.cl.OR.mean$mAccept[ug.cl.OR.mean$stim == "Pre-Stim"], 
                   method = "spearman")
# NS
cor.test(ug.cl.OR.mean$HDRS_z[ug.cl.OR.mean$stim == "Pre-Stim"], 
                   ug.cl.OR.mean$mRT[ug.cl.OR.mean$stim == "Pre-Stim"], 
                   pmethod = "spearman")
# NS
cor.test(ug.cl.OR.mean$HDRS_z[ug.cl.OR.mean$stim == "Pre-Stim"], 
                   ug.cl.OR.mean$mMood[ug.cl.OR.mean$stim == "Pre-Stim"], 
                   method = "spearman")
# No difference in baseline mood

# RL #
######
cor.test(rl.cl.OR.mean$HDRS_z[rl.cl.OR.mean$stim == "Pre-Stim"], 
                   rl.cl.OR.mean$mOpt[rl.cl.OR.mean$stim == "Pre-Stim"], 
                   method = "spearman")
# NS
cor.test(rl.cl.OR.mean$HDRS_z[rl.cl.OR.mean$stim == "Pre-Stim"], 
                   rl.cl.OR.mean$mRT[rl.cl.OR.mean$stim == "Pre-Stim"], 
                   method = "spearman")
# NS

# REGRESSION (average)
lm.ug.mChoice = lm(mAccept ~ HDRS_z * stim, data = ug.cl.OR.mean)
lm.ug.mRT = lm(mRT ~ HDRS_z * stim, data = ug.cl.OR.mean)
lm.ug.mMood = lm(mMood ~ HDRS_z * stim, data = ug.cl.OR.mean)

lm.rl.mChoice = lm(mOpt ~ HDRS_z * stim, data = rl.cl.OR.mean)
lm.rl.mRT = lm(mRT ~ HDRS_z * stim, data = rl.cl.OR.mean)

summary(lm.ug.mChoice) # NS
summary(lm.ug.mRT) # NS
summary(lm.ug.mMood) # NS
summary(lm.rl.mChoice) # NS
summary(lm.rl.mRT) # NS

# REGRESSION (trial-level)
ug.cl.OR.pr = ug.cl.OR %>%
  mutate(accept = 1 * (rej==0), logRT = log(rt))
rl.cl.OR.pr = rl.cl.OR %>%
  mutate(logRT = log(rt))

glmer.ug.mChoice = glmer(accept  ~ HDRS_z * stim + offer_z + (1|idx), data = ug.cl.OR.pr, 
                          family = binomial,
                          control = glmerControl(optimizer='bobyqa') )
lmer.ug.mRT = lmer(logRT ~ HDRS_z * stim + offer_z + (1|idx), data = ug.cl.OR.pr,
                          control = lmerControl(optimizer='bobyqa'))
lmer.ug.mMood = lmer(mood ~ HDRS_z * stim + offer_z + (1|idx), data = ug.cl.OR.pr,
                          control = lmerControl(optimizer='bobyqa'))

glmer.rl.mChoice = glmer(opt ~ HDRS_z * stim + cond + (1|idx), data = rl.cl.OR.pr,
                      family = binomial,
                      control = glmerControl(optimizer='bobyqa') )
lmer.rl.mRT = lmer(logRT ~ HDRS_z * stim + cond + (1|idx), data = rl.cl.OR.pr,
                          control = lmerControl(optimizer='bobyqa'))

summary(glmer.ug.mChoice) # NS
summary(lmer.ug.mRT) # Interaction between HDRS x Stim (Higher HDRS, fastest increase in RT)
summary(lmer.ug.mMood) # NS

summary(glmer.rl.mChoice) # NS
summary(lmer.rl.mRT) # NS


##################
# Testing Fatigue #
###################

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

# Test for fatigue effects in RT within each task
# If fatigue is present, expect increasing RTs over time

# Reversal Learning RT slopes
rl_rt_model = lmer(log_rt ~ trial_z * stim + (1 + stim|idx),
                   data = rl_rt_analysis,
                   REML = FALSE)

# Ultimatum Game RT slopes  
ug_rt_model = lmer(log_rt ~ trial_z * stim + (1 + stim|idx),
                   data = ug_rt_analysis,  
                   REML = FALSE)

summary(rl_rt_model)
summary(ug_rt_model)

eta_squared(rl_rt_model, partial = TRUE) # small effect
eta_squared(ug_rt_model, partial = TRUE) # small effect

# Choice accuracy analysis (for reversal learning)
rl_choice_analysis = rl.EST.pr %>%
  mutate(trial_z = scale(trial)[,1]) %>%
  group_by(idx, stim, trial, trial_z, rew) %>%
  summarise(optimal_choice = mean(opt, na.rm = TRUE),
            .groups = 'drop')

# Test choice accuracy over time
rl_choice_model = glmer(optimal_choice ~ trial_z * stim + (1 + stim|idx),
                        data = rl_choice_analysis,
                        family = binomial,
                        control = glmerControl(optimizer = "bobyqa"))

summary(rl_choice_model)

# Acceptance rate analysis (for ultimatum game)
ug_choice_analysis = ug.EST.pr %>%
  mutate(trial_z = scale(trial)[,1],
         offer_z = scale(offer)[,1],
         accept = 1*(rej==0)) %>%
  group_by(idx, stim, trial, trial_z, offer_z) %>%
  summarise(accept_rate = mean(accept, na.rm = TRUE),
            .groups = 'drop')

# Test acceptance rate over time
ug_choice_model = glmer(accept_rate ~  trial_z * stim + (1 + stim|idx),
                        data = ug_choice_analysis,
                        family = binomial,
                        control = glmerControl(optimizer = "bobyqa"))

summary(ug_choice_model)

# Visualization of behavioral measures over time

fig.rl.rt = rl_rt_analysis %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         color = ifelse(stim == "Pre-Stim", "#fcbba1", "#a50f15")) %>%
  ggplot(aes(x = trial, y = log_rt, fill = color, group = stim)) +
  theme_pubr(base_size = 14) +
  geom_hline(yintercept = mean(rl_rt_analysis$log_rt, na.rm = TRUE), 
             linetype = "dashed", linewidth = 1, color = "gray") +
  stat_summary(aes(color = color), geom = "line", size = 1, linewidth = 1.5) +
  stat_summary(geom = "ribbon", alpha = 0.35, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_identity() +
  labs(y = "Response time [log sec]", x = "Trial", title = "Reversal earning") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA)) +
  guides(color = "none", fill = "none")

fig.ug.rt = ug_rt_analysis %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         color = ifelse(stim == "Pre-Stim", "#9ecae1", "#08519c")) %>%
  ggplot(aes(x = trial, y = log_rt, fill = color, group = stim)) +
  theme_pubr(base_size = 14) +
  geom_hline(yintercept = mean(ug_rt_analysis$log_rt, na.rm = TRUE), 
             linetype = "dashed", linewidth = 1, color = "gray") +
  stat_summary(aes(color = color), geom = "line", size = 1, linewidth = 1.5) +
  stat_summary(geom = "ribbon", alpha = 0.35, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_identity() +
  labs(y = "Response time [log sec]", x = "Trial", title = "Ultimatum game") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA)) +
  guides(color = "none", fill = "none")

fig.rl.c = rl_choice_analysis %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         color = ifelse(stim == "Pre-Stim", "#fcbba1", "#a50f15")) %>%
  ggplot(aes(x = trial, y = optimal_choice, fill = color, group = stim)) +
  theme_pubr(base_size = 14) +
  geom_hline(yintercept = mean(rl_choice_analysis$optimal_choice, na.rm = TRUE), 
             linetype = "dashed", linewidth = 1, color = "gray") +
  stat_summary(aes(color = color), geom = "line", size = 1, linewidth = 1.5) +
  stat_summary(geom = "ribbon", alpha = 0.35, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_identity() +
  labs(y = "p(Optimal)", x = "Trial", title = "Reversal earning") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA)) +
  guides(color = "none", fill = "none")

fig.ug.c = ug_choice_analysis %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         color = ifelse(stim == "Pre-Stim", "#9ecae1", "#08519c")) %>%
  ggplot(aes(x = trial, y = accept_rate, fill = color, group = stim)) +
  theme_pubr(base_size = 14) +
  geom_hline(yintercept = mean(ug_choice_analysis$accept_rate, na.rm = TRUE), 
             linetype = "dashed", linewidth = 1, color = "gray") +
  stat_summary(aes(color = color), geom = "line", size = 1, linewidth = 1.5) +
  stat_summary(geom = "ribbon", alpha = 0.35, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_identity() +
  labs(y = "p(Accept)", x = "Trial", title = "Ultimatum game") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA)) +
  guides(color = "none", fill = "none")

# Create combined behavioral fatigue figure
behavioral_fatigue_fig = (fig.rl.rt + fig.ug.rt ) / (fig.rl.c + fig.ug.c) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") & 
  theme(legend.position = 'right')

# Extract slope estimates for interpretation
rl_rt_slopes <- emtrends(rl_rt_model, ~ stim, var = "trial_z")
ug_rt_slopes <- emtrends(ug_rt_model, ~ stim, var = "trial_z")
rl_c_slopes <- emtrends(rl_choice_model, ~ stim, var = "trial_z")
ug_c_slopes <- emtrends(ug_choice_model, ~ stim, var = "trial_z")

summary(rl_rt_slopes)
summary(ug_rt_slopes)

# Test if slopes are significantly different from zero (fatigue effect)
test(rl_rt_slopes, null = 0)
test(ug_rt_slopes, null = 0)

test(rl_c_slopes, null = 0)
test(ug_c_slopes, null = 0)


##
# Alternative legend for all NT #

create_prestim_legend <- function() {
  # Dummy data for Pre-Stim legend
  legend_data_pre <- data.frame(
    nt = factor(c("DA", "5-HT", "NE"), levels = c("DA", "5-HT", "NE")),
    color = c("#fcbba1", "#9ecae1", "#c7e9c0"),
    x = 1,
    y = 1
  )
  
  ggplot(legend_data_pre, aes(x = x, y = y, color = nt)) +
    geom_point(size = 0) +  # Invisible points
    geom_line(size = 2) +   # This creates the line legend
    scale_color_manual(
      name = "Pre-Stim",
      values = c("DA" = "#fcbba1", "5-HT" = "#9ecae1", "NE" = "#c7e9c0")
    ) +
    theme_void() +
    theme(
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 10),
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.margin = margin(0, 0, 10, 0)
    ) +
    guides(color = guide_legend(
      override.aes = list(size = 2, alpha = 1),
      title.position = "top",
      title.hjust = 0.5
    ))
}

create_poststim_legend <- function() {
  # Dummy data for Post-Stim legend
  legend_data_post <- data.frame(
    nt = factor(c("DA", "5-HT", "NE"), levels = c("DA", "5-HT", "NE")),
    color = c("#a50f15", "#08519c", "#006d2c"),
    x = 1,
    y = 1
  )
  
  ggplot(legend_data_post, aes(x = x, y = y, color = nt)) +
    geom_point(size = 0) +  # Invisible points
    geom_line(size = 2) +   # This creates the line legend
    scale_color_manual(
      name = "Post-Stim",
      values = c("DA" = "#a50f15", "5-HT" = "#08519c", "NE" = "#006d2c")
    ) +
    theme_void() +
    theme(
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 10),
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.margin = margin(10, 0, 0, 0)
    ) +
    guides(color = guide_legend(
      override.aes = list(size = 2, alpha = 1),
      title.position = "top",
      title.hjust = 0.5
    ))
}

# Extract the legends
prestim_legend <- get_legend(create_prestim_legend())
poststim_legend <- get_legend(create_poststim_legend())

# Stack the legends vertically
stacked_legends <- plot_grid(prestim_legend, poststim_legend, 
                             ncol = 1, align = "v")


ggsave(res_dir / "neurotransmitter_legend.png", 
       plot = stacked_legends,
       device = "png",
       width = 4,          # Width in inches
       height = 4,         # Height in inches  
       dpi = 300)    
############################################################
# Does Baseline HDRS predict context-specific changes in NT?
############################################################
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

task.means = task.merge.baseline %>%
  group_by(idx, task, nt, stim, HDRS_z, HDRS) %>%
  summarize(Oz_mean = mean(Oz, na.rm = TRUE)) %>%
  ungroup()

# Pivot to create difference scores
task.wide = task.means %>%
  pivot_wider(names_from = task, 
              values_from = Oz_mean,
              names_prefix = "Oz_") %>%
  mutate(task_diff = Oz_UG - Oz_RL) 

# Model HDRS as outcome, with task differences as predictors
hdrs.model.DA = lm(HDRS ~ task_diff * stim,
                     data = task.wide[task.wide$nt == "DA",])

hdrs.model.SE = lm(HDRS ~ task_diff * stim,
                     data = task.wide[task.wide$nt == "SE",])

hdrs.model.NE = lm(HDRS ~ task_diff * stim,
                     data = task.wide[task.wide$nt == "NE",])

summary(hdrs.model.DA)
summary(hdrs.model.SE)
summary(hdrs.model.NE)

# What differentiates participants who show the expected pattern vs. those who do not?
task.means = task.merge.baseline %>%
  group_by(idx, task, nt, stim, HDRS_z, HDRS) %>%
  summarize(Oz_mean = mean(Oz, na.rm = TRUE), O_mean = mean(O, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(HDRS_M6 = cl.df.m6$HDRS[match(idx, cl.df.m6$idx)],
         HDRS_M6_Change = HDRS_M6 - HDRS,
         HDRS_M6_z = scale(HDRS_M6)[,1])

# Look at ID and HDRS CHange with OZ change patterns
task.means %>% select(idx, task, nt, stim, Oz_mean, HDRS_M6_Change) %>% distinct() %>%
  filter(nt!="NE") %>%
  group_by(idx, task, nt) %>%
  mutate(Post_minus_Pre = Oz_mean[stim == "Post-Stim"] - Oz_mean[stim == "Pre-Stim"]) %>%
  select(!c(Oz_mean,stim)) %>%
  distinct() %>%
  pivot_wider(names_from = c(task, nt), values_from = Post_minus_Pre, names_prefix = "Oz_")
  


# Pull out the post-stim Oz_mean for each subject
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
  
task.means %>%
  group_by(task,nt,stim,expected_pattern) %>%
  summarise(mean_HDRS = mean(HDRS, na.rm = TRUE),
            mean_Oz = mean(Oz_mean, na.rm = TRUE),
            count = n())

# Visualize HDRS by expected pattern
task.means %>% select(idx, task, nt, HDRS, expected_pattern) %>%
   distinct() %>%
  ggplot(aes(x = expected_pattern, y = HDRS, color = expected_pattern)) +
  geom_point(size = 1.5, position = position_dodge2(width = 0.2)) +
  stat_summary(size = 1, linewidth = 1) +
  facet_grid(task ~ nt) +
  theme_pubr() +
  labs(color = "Pattern") +
  scale_color_manual(values = c("Expected" = "#2ca25f", "Unexpected" = "#de2d26")) +
  theme(legend.position = "bottom") +
  stat_compare_means(method = "t.test", label = "p.signif") +
  labs(y = "Baseline HDRS", x = "Pattern of Change")
# Higher baseline, but Nothing significant!

task.means %>% select(idx, task, nt, HDRS_M6, expected_pattern) %>%
   distinct() %>%
  ggplot(aes(x = expected_pattern, y = HDRS_M6, color = expected_pattern)) +
  geom_point(size = 1.5, position = position_dodge2(width = 0.2)) +
  stat_summary(size = 1, linewidth = 1) +
  facet_grid(task ~ nt) +
  theme_pubr() +
  labs(color = "Pattern") +
  scale_color_manual(values = c("Expected" = "#2ca25f", "Unexpected" = "#de2d26")) +
  theme(legend.position = "bottom") +
  stat_compare_means(method = "t.test", label = "p.signif") +
  labs(y = "Month 6 HDRS", x = "Pattern of Change")

# Visualize O_mean by expected pattern (pre-stim only)
task.means %>% 
  filter(stim == "Pre-Stim") %>%
  select(idx, task, nt, O_mean, expected_pattern) %>%
  distinct() %>%
  ggplot(
  aes(x = expected_pattern, y = O_mean, color = expected_pattern)) +
  stat_summary() +
  facet_grid(task ~ nt) +
  theme_pubr() +
  labs(color = "Pattern") +
  scale_color_manual(values = c("Expected" = "#2ca25f", "Unexpected" = "#de2d26")) +
  theme(legend.position = "bottom") +
  stat_compare_means(method = "t.test", label = "p.signif")

task.means %>% 
  filter(stim == "Post-Stim") %>%
  select(idx, task, nt, O_mean, expected_pattern) %>%
  distinct() %>%
  ggplot(
  aes(x = expected_pattern, y = O_mean, color = expected_pattern)) +
  stat_summary() +
  facet_grid(task ~ nt) +
  theme_pubr() +
  labs(color = "Pattern") +
  scale_color_manual(values = c("Expected" = "#2ca25f", "Unexpected" = "#de2d26")) +
  theme(legend.position = "bottom") +
  stat_compare_means(method = "t.test", label = "p.signif")

task.means %>%
  select(idx, task, expected_pattern, HDRS_M6) %>%
  mutate(expected_pattern = factor(expected_pattern, levels = c("Unexpected", "Expected"))) %>%
  distinct() %>%
  filter(task == "UG") %>%
  glm(HDRS_M6 ~ expected_pattern, data = ., family = gaussian) %>%
  summary()
# No effect - but marginally lower HDRS at 6 months in those with expected pattern (increased DA and SE)

task.means %>%
  select(idx, task, expected_pattern, HDRS) %>%
  mutate(expected_pattern = factor(expected_pattern, levels = c("Unexpected", "Expected"))) %>%
  distinct() %>%
  filter(task == "UG") %>%
  glm(HDRS ~ expected_pattern, data = ., family = gaussian) %>%
  summary()
# No difference in baseline HDRS between those with vs. without the expected pattern

task.means %>%
  select(idx, task, expected_pattern, HDRS, HDRS_M6) %>%
  distinct() %>%
  mutate(HDRS_change = HDRS_M6 - HDRS) %>%
  ggplot(aes(x = expected_pattern, y = HDRS_change, color = expected_pattern)) +
  stat_summary(size = 1, linewidth = 1) +
  geom_point(position = position_dodge2(width = 0.2)) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  facet_wrap(~ task) +
  scale_color_manual(values = c("Expected" = "#2ca25f", "Unexpected" = "#de2d26")) +
  theme_pubr()

task.means %>%
  select(idx, task, expected_pattern, HDRS, HDRS_M6) %>%
  distinct() %>%
  mutate(HDRS_change = HDRS_M6 - HDRS) %>%
  filter(task == "UG") %>%
  glm(HDRS_change ~ expected_pattern, data = ., family = gaussian) %>%
  summary()
# But this is significant! 
# Those with the expected pattern (increased DA and SE) show greater reductions in HDRS over 6 months, compared to those without the expected pattern.


# Look at behavior

ug.EST.pr = merge(ug.EST.pr, post_stim_vals[post_stim_vals$task == "UG",], by = c("idx", "nt"))
rl.EST.pr = merge(rl.EST.pr, post_stim_vals[post_stim_vals$task == "RL",], by = c("idx", "nt"))

ug.beh = ug.EST.pr %>% 
  filter(nt == "DA") %>%
  group_by(idx,expected_pattern,nt,stim) %>%
  summarise(mAccept = mean(rej==0), mRT = mean(log(rt)))

rl.beh = rl.EST.pr %>% 
  filter(nt == "DA") %>%
  group_by(idx,expected_pattern,nt,stim) %>%
  summarise(mOpt = mean(opt==0), mRT = mean(log(rt)))

ug.choice.aov = aov(mAccept ~ expected_pattern * stim, data = ug.beh)
ug.rt.aov = aov(mRT ~ expected_pattern * stim, data = ug.beh)

summary(ug.choice.aov)
summary(ug.rt.aov)

rl.choice.aov = aov(mOpt ~ expected_pattern * stim, data = rl.beh)
rl.rt.aov = aov(mRT ~ expected_pattern * stim, data = rl.beh)

summary(rl.choice.aov)
summary(rl.rt.aov)


ug.EST.pr %>% 
  group_by(expected_pattern,nt,stim) %>%
  summarise(mAccept = mean(rej==0)) %>%
  ggplot(aes(x = expected_pattern, y = mAccept, color = expected_pattern)) +
  stat_summary() +
  facet_grid( ~ nt) +
  theme_pubr() +
  labs(color = "Pattern") +
  scale_color_manual(values = c("Expected" = "#2ca25f", "Unexpected" = "#de2d26")) +
  theme(legend.position = "bottom") +
  stat_compare_means(method = "t.test", label = "p.signif")
  # Those who increased in 5-HT tended to accept more offers

 

ug.EST.pr %>% 
  group_by(expected_pattern,nt,stim) %>%
  summarise(mRT = mean(log(rt))) %>%
  ggplot(aes(x = expected_pattern, y = mRT, color = expected_pattern)) +
  stat_summary() +
  facet_grid( ~ nt) +
  theme_pubr() +
  labs(color = "Pattern") +
  scale_color_manual(values = c("Expected" = "#2ca25f", "Unexpected" = "#de2d26")) +
  theme(legend.position = "bottom") +
  stat_compare_means(method = "t.test", label = "p.signif")
# No significant differences in RT

summary(lmer(log(rt) ~ expected_pattern * stim + offer_z + (1+offer_z|idx), data = ug.EST.pr, family = binomial))  


rl.EST.pr %>% 
  group_by(expected_pattern,nt,stim) %>%
  summarise(mOpt= mean(opt==0)) %>%
  ggplot(aes(x = expected_pattern, y = mOpt, color = expected_pattern)) +
  stat_summary() +
  facet_grid( ~ nt) +
  theme_pubr() +
  labs(color = "Pattern") +
  scale_color_manual(values = c("Expected" = "#2ca25f", "Unexpected" = "#de2d26")) +
  theme(legend.position = "bottom") +
  stat_compare_means(method = "t.test", label = "p.signif")
  # Those who increased in DA tended to choose more optimally

rl.EST.pr %>% 
  group_by(expected_pattern,nt,stim) %>%
  summarise(mRT = mean(log(rt))) %>%
  ggplot(aes(x = expected_pattern, y = mRT, color = expected_pattern)) +
  stat_summary() +
  facet_grid( ~ nt) +
  theme_pubr() +
  labs(color = "Pattern") +
  scale_color_manual(values = c("Expected" = "#2ca25f", "Unexpected" = "#de2d26")) +
  theme(legend.position = "bottom") +
  stat_compare_means(method = "t.test", label = "p.signif")
# No significant differences in RT

## 
## Other clinical varianbles

# junel fe = norethindrone acetate/ethinyl estradiol/ferrous fumarate -> birth control
# Levothyroxine = thyroid hormone replacement
# Lisinopril-hydrochlorothiazide = high blood pressure
# venlafaxine = SNRI
# lisdexamfetamine = stimulant for ADHD
# busprione = for anxiety
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


# Get all unique drugs
all_drugs <- cl.alt %>%
  filter(!is.na(meds)) %>%
  separate_rows(meds, sep = ",") %>%
  mutate(meds = trimws(meds)) %>%
  filter(meds != "") %>%
  pull(meds) %>%
  unique()

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

# Drugs where multiple participants are taking them
multi_drugs = cl.alt.binary %>%
  select(-idx) %>%
  summarise(across(everything(), sum)) %>%
  pivot_longer(everything(), names_to = "drug", values_to = "n_participants") %>%
  filter(n_participants > 1) %>%
  arrange(desc(n_participants))

cl.alt.stim = post_stim_vals %>%
  select(idx,task,expected_pattern) %>%
  left_join(cl.alt.binary, by = "idx") %>%
  mutate(expected_pattern_binary = ifelse(expected_pattern == "Expected", 1, 0))

  

# See if any medications are associated with expected vs. unexpected pattern of change
med_models = list()
for (task in c("UG","RL") ) {
  for (drug in multi_drugs$drug) {
    model <- glm(expected_pattern_binary ~ get(drug), data = cl.alt.stim[cl.alt.stim$task == task,], 
    family = binomial)
    med_models[[paste0(task, "_", drug)]] <- summary(model)
  } 
}

# $UG_levothyroxine, $UG_lamotrigine (almost significant)

# Look at actual values
post_stim_oz <- task.means %>%
  select(idx, task, nt, stim, Oz_mean) %>%
  group_by(idx,task,nt) %>%
  summarise(Oz_change = Oz_mean[stim == "Post-Stim"] - Oz_mean[stim == "Pre-Stim"],
            Oz_pre = Oz_mean[stim == "Pre-Stim"]) %>%
  ungroup()

cl.alt.stim_oz = post_stim_oz %>%
  left_join(cl.alt.binary, by = "idx") %>%
  pivot_wider(names_from = nt, values_from = Oz_change, names_prefix = "Oz_")

cl.alt.stim_oz %>%
 select(idx,task,Oz_DA,Oz_SE,Oz_NE)

med_models_oz = list()
for (task in c("UG","RL") ) {
  for (drug in multi_drugs$drug) {
    model <- glm(get(drug) ~ Oz_SE * Oz_DA, data = cl.alt.stim_oz[cl.alt.stim_oz$task == task,], 
    family = binomial)
    med_models_oz[[paste0(task, "_", drug)]] <- summary(model)
  } 
}
# No pattern of DA/SE increase associated with any medications

# Look at individual NTs
med_models_SE = list()
for (task in c("UG","RL") ) {
  for (drug in multi_drugs$drug) {
    model <- glm(get(drug) ~ Oz_SE, data = cl.alt.stim_oz[cl.alt.stim_oz$task == task,], 
    family = binomial)
    med_models_SE[[paste0(task, "_", drug)]] <- summary(model)
  } 
}
# Nothing significant

med_models_DA = list()
for (task in c("UG","RL") ) {
  for (drug in multi_drugs$drug) {
    model <- glm(get(drug) ~ Oz_DA, data = cl.alt.stim_oz[cl.alt.stim_oz$task == task,], 
    family = binomial)
    med_models_DA[[paste0(task, "_", drug)]] <- summary(model)
  } 
}

med_models_NE = list()
for (task in c("UG","RL") ) {
  for (drug in multi_drugs$drug) {
    model <- glm(get(drug) ~ Oz_NE, data = cl.alt.stim_oz[cl.alt.stim_oz$task == task,], 
    family = binomial)
    med_models_NE[[paste0(task, "_", drug)]] <- summary(model)
  } 
}

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

print(class_counts)


# Define classes with enough participants (n > 1, ideally > 2)
multi_classes <- class_counts %>%
  filter(n_participants > 1) %>%
  pull(class)

# Join with stim data
cl.alt.stim_class <- post_stim_vals %>%
  select(idx, task, expected_pattern) %>%
  left_join(cl.alt.class, by = "idx") %>%
  mutate(expected_pattern_binary = ifelse(expected_pattern == "Expected", 1, 0))

# Join with Oz change data
cl.alt.stim_oz_class <- post_stim_oz %>%
  left_join(cl.alt.class, by = "idx") %>%
  pivot_wider(names_from = nt, values_from = c(Oz_change, Oz_pre), names_prefix = "Oz_")

colnames(cl.alt.stim_oz_class)

# Look at outliers
cl.alt.stim_oz_class %>%
  filter(idx %in% c(4,7)) %>%
  as.data.frame()

#
cl.alt.stim_oz_class %>%
  filter(SNRI == 1, Benzo == 1, Stimulant == 1) %>% as.data.frame()

  cl.alt.stim_oz_class %>%
  filter(Thyroid == 1) %>% as.data.frame()

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

# --- Drug class ~ Oz_SE + Oz_DA (interaction) ---
class_models_oz <- list()
for (task in c("UG", "RL")) {
  for (drug_class in multi_classes) {
    model <- glm(get(drug_class) ~ Oz_change_Oz_SE * Oz_change_Oz_DA,
                 data = cl.alt.stim_oz_class[cl.alt.stim_oz_class$task == task, ],
                 family = binomial)
    class_models_oz[[paste0(task, "_", drug_class)]] <- summary(model)
  }
}

# --- Individual NTs ---
for (nt in c("Oz_change_Oz_DA", "Oz_change_Oz_SE", "Oz_change_Oz_NE","Oz_pre_Oz_SE","Oz_pre_Oz_DA","Oz_pre_Oz_NE")) {
  assign(paste0("class_models_", nt), {
    models <- list()
    for (task in c("UG", "RL")) {
      for (drug_class in multi_classes) {
        model <- glm(get(drug_class) ~ get(nt),
                     data = cl.alt.stim_oz_class[cl.alt.stim_oz_class$task == task, ],
                     family = binomial)
        models[[paste0(task, "_", drug_class)]] <- summary(model)
      }
    }
    models
  })
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

# Apply to each model list
results_pattern <- extract_glm_results(class_models_pattern)
results_oz      <- extract_glm_results(class_models_oz)
results_SE      <- extract_glm_results(class_models_Oz_SE)
results_DA      <- extract_glm_results(class_models_Oz_DA)
results_NE      <- extract_glm_results(class_models_Oz_NE)

# View top hits
results_pattern %>% filter(p_value < 0.1)
results_oz %>% filter(p_value < 0.1) # NS

# Collect all results into one tidy table
nt_vars <- c("Oz_change_Oz_DA", "Oz_change_Oz_SE", "Oz_change_Oz_NE",
             "Oz_pre_Oz_SE",    "Oz_pre_Oz_DA",    "Oz_pre_Oz_NE")

all_class_results <- map_dfr(nt_vars, function(nt) {
  model_list <- get(paste0("class_models_", nt))
  extract_glm_results(model_list, nt_var = nt)
}) %>%
  # Parse nt_var into type (change vs pre) and neurotransmitter
  mutate(
    stim_period = case_when(
      str_detect(nt, "Oz_change") ~ "change",
      str_detect(nt, "Oz_pre")   ~ "pre",
      TRUE ~ NA_character_
    ),
    neurotransmitter = str_extract(nt, "(DA|SE|NE)$")
  ) %>%
  select(model, task, predictor, stim_period, neurotransmitter, estimate, se, z, p_value)

# Quick views
all_class_results %>% filter(p_value < 0.1)

# Split by stim period
all_class_results %>%
  select(!model) %>%
  filter(stim_period == "pre") %>%
  arrange(p_value)

all_class_results %>%
  select(!model) %>%
  filter(stim_period == "change") %>%
  arrange(p_value)

# Explore what's driving the change in Expected/Unexpected pattern - is it the absolute change in DA/SE, the pre-stim levels, or some combination?
cl.alt.stim_oz_class <- cl.alt.stim_oz_class %>%
  mutate(
    DA_SE_product = Oz_change_Oz_DA * Oz_change_Oz_SE,
    DA_SE_min     = pmin(Oz_change_Oz_DA, Oz_change_Oz_SE),  # the "weakest link" - how much did BOTH increase
    DA_SE_diff    = Oz_change_Oz_DA - Oz_change_Oz_SE        # relative dominance
  )

new_vars <- c("DA_SE_product", "DA_SE_min", "DA_SE_diff")

new_models <- map_dfr(new_vars, function(var) {
  models <- list()
  for (task in c("UG", "RL")) {
    for (drug_class in multi_classes) {
      model <- glm(get(drug_class) ~ get(var),
                   data = cl.alt.stim_oz_class[cl.alt.stim_oz_class$task == task, ],
                   family = binomial)
      models[[paste0(task, "_", drug_class)]] <- summary(model)
    }
  }
  extract_glm_results(models, nt_var = var)
}) %>%
  arrange(p_value)

new_models %>%
  arrange(p_value)

post_stim_vals %>%
  filter(task == "UG", nt == "DA") %>%
  dplyr::select(idx, task, expected_pattern) %>%
  left_join(cl.alt.class %>% dplyr::select(idx, Stimulant, SNRI, Thyroid), 
            by = "idx") %>% 
  dplyr::summarise(
    expected             = sum(expected_pattern == "Expected"),
    expected_SNRI        = sum(expected_pattern == "Expected" & SNRI == 1),
    expected_Stimulant   = sum(expected_pattern == "Expected" & Stimulant == 1),
    expected_thyroid     = sum(expected_pattern == "Expected" & Thyroid == 1),
    unexpected_SNRI      = sum(expected_pattern == "Unexpected" & SNRI == 1),
    unexpected_Stimulant = sum(expected_pattern == "Unexpected" & Stimulant == 1),
    unexpected_thyroid   = sum(expected_pattern == "Unexpected" & Thyroid == 1),
    unexpected_any = sum(expected_pattern == "Unexpected" & (SNRI == 1 | Stimulant == 1 | Thyroid == 1))
  ) %>% as.data.frame()


post_stim_vals %>%
  filter(task == "UG", nt == "DA") %>%
  dplyr::select(idx, task, expected_pattern,DA_change,SE_change) %>%
  left_join(cl.alt.class %>% dplyr::select(idx, Stimulant, SNRI, Thyroid), 
            by = "idx") %>% 
  dplyr::summarise(
    DA_Decr_SNRI        = sum(DA_change == "Decreased" & SNRI == 1),
    SE_Decr_SNRI       = sum(SE_change == "Decreased" & SNRI == 1),
    DA_Decr_Stimulant   = sum(DA_change == "Decreased" & Stimulant == 1),
    SE_Decr_Stimulant   = sum(SE_change == "Decreased" & Stimulant == 1),
    DA_Decr_thyroid     = sum(DA_change == "Decreased" & Thyroid == 1),
    SE_Decr_thyroid     = sum(SE_change == "Decreased" & Thyroid == 1),
    DA_Decr_any = sum(DA_change == "Decreased" & (SNRI == 1 | Stimulant == 1 | Thyroid == 1)),
    SE_Decr_any = sum(SE_change == "Decreased" & (SNRI == 1 | Stimulant == 1 | Thyroid == 1))
  ) %>% as.data.frame()
# Not seemingly particular to one NT


# Summary: chronic monoamine-elevating medications may disrupt the coordinated DA+SE response 
# to subthreshold stimulation without reducing the overall magnitude of NT change

# Assess Medication class × expected pattern × outcome

# Requires one row per participant per task
med_df <- task.means %>%
  filter(task == "UG") %>%  
  mutate(expected_pattern = factor(expected_pattern, levels = c("Unexpected", "Expected"))) %>%
  dplyr::select(idx, expected_pattern, HDRS, HDRS_M6) %>%
  distinct() %>%
  mutate(HDRS_change = HDRS_M6 - HDRS) %>%
  left_join(cl.alt.class %>% dplyr::select(idx, Stimulant, SNRI, Thyroid), by = "idx")

lm(HDRS_change ~ expected_pattern * Stimulant, data = med_df) %>% summary() # NS
lm(HDRS_change ~ expected_pattern * SNRI, data = med_df) %>% summary() #NS 
lm(HDRS_change ~ expected_pattern * Thyroid, data = med_df) %>% summary() 
# Significant interaction - those on thyroid meds show less improvement if they have the expected pattern (increased DA and SE) 
#- but small n (n=3 on thyroid meds)

lm(HDRS_M6 ~ expected_pattern * Stimulant, data = med_df) %>% summary() # NS
lm(HDRS_M6 ~ expected_pattern * SNRI, data = med_df) %>% summary() #NS 
lm(HDRS_M6 ~ expected_pattern * Thyroid, data = med_df) %>% summary() # Very significant

model <- lm(HDRS_M6 ~ expected_pattern * Thyroid, data = med_df)
ggpredict(model, terms = c("expected_pattern", "Thyroid [0,1]")) %>%
  plot() +
  scale_color_manual(
    values = c("0" = "#2ca25f", "1" = "#de2d26"),
    labels = c("0" = "No Thyroid Medication", "1" = "Thyroid Medication")
  ) +
  labs(
    x = "NT Response Pattern",
    y = "Month 6 HDRS",
    color = "Thyroid Medication",
    title = "Thyroid Medication × NT Pattern Interaction on M6 Depression Severity"
  ) +
  theme_pubr()

# Combined medication burden score
med_df <- med_df %>%
  mutate(med_burden = Stimulant + SNRI + Thyroid)  # 0, 1, 2, or 3

# Interaction models
lm(HDRS_change ~ expected_pattern * med_burden, data = med_df) %>% summary()
# NS

# Look at other therepeutic experiences
cl.alt.therapies = data.frame(
  idx = c(1:10),
  ECT = c(1, 1, 1, 1, 1, 1, 0, 1, 1, 1),
  therapy = c(0,1,1,1,1,1,1,1,1,1),
  TMS = c(1,0,1,0,1,1,1,1,0,1),
  ketamine = c(1,1,1,1,1,1,1,0,1,1)
)

cl.alt.stim_therapies <- post_stim_vals %>%
  dplyr::select(idx, task, expected_pattern) %>% distinct() %>%
  left_join(cl.alt.therapies, by = "idx") %>%
  mutate(expected_pattern_binary = ifelse(expected_pattern == "Expected", 1, 0))

cl.alt.oz_therapies <- post_stim_oz %>%
  left_join(cl.alt.therapies, by = "idx") %>%
  pivot_wider(names_from = nt, values_from = c(Oz_change, Oz_pre), names_prefix = "Oz_")

glm(ECT ~ expected_pattern_binary, data = cl.alt.stim_therapies[cl.alt.stim_therapies$task == "UG",], family = binomial) %>% summary()
glm(therapy ~ expected_pattern_binary, data = cl.alt.stim_therapies[cl.alt.stim_therapies$task == "UG",], family = binomial) %>% summary()
glm(TMS ~ expected_pattern_binary, data = cl.alt.stim_therapies[cl.alt.stim_therapies$task == "UG",], family = binomial) %>% summary()
glm(ketamine ~ expected_pattern_binary, data = cl.alt.stim_therapies[cl.alt.stim_therapies$task == "UG",], family = binomial) %>% summary()
# Nothing related to therapy

lm(ECT ~ Oz_change_Oz_SE * Oz_change_Oz_DA, data = cl.alt.oz_therapies[cl.alt.oz_therapies$task == "UG",]) %>% summary()
lm(therapy ~ Oz_change_Oz_SE * Oz_change_Oz_DA, data = cl.alt.oz_therapies[cl.alt.oz_therapies$task == "UG",]) %>% summary()
lm(TMS ~ Oz_change_Oz_SE * Oz_change_Oz_DA, data = cl.alt.oz_therapies[cl.alt.oz_therapies$task == "UG",]) %>% summary()

# signicant effect of ketamine
lm(ketamine ~ Oz_change_Oz_SE * Oz_change_Oz_DA, data = cl.alt.oz_therapies[cl.alt.oz_therapies$task == "UG",]) %>% summary()

lm(ECT ~ Oz_pre_Oz_SE * Oz_pre_Oz_DA, data = cl.alt.oz_therapies[cl.alt.oz_therapies$task == "UG",]) %>% summary()
lm(therapy ~ Oz_pre_Oz_SE * Oz_pre_Oz_DA, data = cl.alt.oz_therapies[cl.alt.oz_therapies$task == "UG",]) %>% summary()
lm(TMS ~ Oz_pre_Oz_SE * Oz_pre_Oz_DA, data = cl.alt.oz_therapies[cl.alt.oz_therapies$task == "UG",]) %>% summary()
lm(ketamine ~ Oz_pre_Oz_SE * Oz_pre_Oz_DA, data = cl.alt.oz_therapies[cl.alt.oz_therapies$task == "UG",]) %>% summary()
# Slightly related to baseline

# BUT CANOT INTERPRET BECAUSE ONLY ONE PARTICIPANT WITHOUT PRIOR KETAMINE

# Check overlap between ketamine and the three drug classes
cl.alt.therapies %>%
  left_join(cl.alt.class %>% dplyr::select(idx, Stimulant, SNRI, Thyroid), by = "idx") %>%
  dplyr::select(idx, ketamine, Stimulant, SNRI, Thyroid)

# How many patients are on ketamine AND each drug class?
cl.alt.therapies %>%
  left_join(cl.alt.class %>% dplyr::select(idx, Stimulant, SNRI, Thyroid), by = "idx") %>%
  dplyr::summarise(
    n_ketamine         = sum(ketamine),
    ket_and_stimulant  = sum(ketamine == 1 & Stimulant == 1),
    ket_and_SNRI       = sum(ketamine == 1 & SNRI == 1),
    ket_and_thyroid    = sum(ketamine == 1 & Thyroid == 1),
    no_ket_stimulant   = sum(ketamine == 0 & Stimulant == 1),
    no_ket_SNRI        = sum(ketamine == 0 & SNRI == 1),
    no_ket_thyroid     = sum(ketamine == 0 & Thyroid == 1),
    thyroid = sum(Thyroid == 1),
  )
