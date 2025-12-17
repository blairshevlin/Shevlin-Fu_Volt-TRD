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
# 2025/11/4    Blair Shevlin                          adapted code from supplementary_analyses.r to generate figure supplement 2

rm(list = ls())

library(tidyverse)
library(fs)
library(here)
library(patchwork)
library(ggpubr)
library(RColorBrewer)
library(lmerTest)
library(emmeans)
library(cowplot)

# Paths
dir = path(here())
nt_dir = dir / "data" / "nt" / "processed"
beh_dir = dir / "data" / "behavior"
res_dir = dir / "results" # Updated

# Subject-level NT Data
load(nt_dir / "UG_RL_NT-Continuous_9-23-25.RData")

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
  select(idx,stim,nt,Oz,trial_z) %>%
  mutate(task = "RL")
ug.premerge = ug.EST.pr %>%
  select(idx,stim,nt,Oz,trial_z) %>%
  mutate(task = "UG")

task.merge = rbind(rl.premerge,ug.premerge) %>% 
            mutate(task = factor(task))

sess.lm.DA = lmer(Oz ~ stim * task * trial_z +
                    (1+stim+task|idx),
               data = task.merge[task.merge$nt == "DA",],
               REML=F,
               control = lmerControl(optimizer = "nmkbw"))
sess.lm.SE = lmer(Oz ~ stim * task * trial_z +
                    (1+stim+task|idx),
                  data = task.merge[task.merge$nt == "SE",],
                  REML=F,
                  control = lmerControl(optimizer = "nmkbw"))
sess.lm.NE = lmer(Oz ~ stim * task * trial_z +
                    (1+stim+task|idx),
                  data = task.merge[task.merge$nt == "NE",],
                  REML=F,
                  control = lmerControl(optimizer = "nmkbw"))

summary(sess.lm.DA)
summary(sess.lm.SE)
summary(sess.lm.NE)

# Simple slopes analysis for stim effects
da_stim <- emmeans(sess.lm.DA, ~ stim * task)
se_stim <- emmeans(sess.lm.SE, ~ stim * task)
ne_stim <- emmeans(sess.lm.NE, ~ stim * task)

pairs(da_stim, by = "task")
pairs(se_stim, by = "task")
pairs(ne_stim, by = "task")


# Simple slopes analysis for the three-way interactions

# 1. Simple slopes for DA model
# Get slopes for each stim*task combination
da_slopes <- emtrends(sess.lm.DA, ~ stim * task, var = "trial_z")
summary(da_slopes)

# Test differences between slopes
pairs(da_slopes, by = "stim")  # Compare tasks within each stim condition
pairs(da_slopes, by = "task")  # Compare stim within each task

# 2. Simple slopes for SE model  
se_slopes <- emtrends(sess.lm.SE, ~ stim * task, var = "trial_z")
summary(se_slopes)
pairs(se_slopes, by = "stim") # 
pairs(se_slopes, by = "task") # NS

# 3. Simple slopes for NE model
ne_slopes <- emtrends(sess.lm.NE, ~ stim * task, var = "trial_z")
summary(ne_slopes)
pairs(ne_slopes, by = "stim") #
pairs(ne_slopes, by = "task")

# 4. Effect sizes at early vs late trials
trial_q = as.numeric(quantile(task.merge$trial_z, probs = c(.25, .75)))

# Early trials
early_effects_DA <- emmeans(sess.lm.DA, ~ stim * task, at = list(`trial_z` = min(trial_q)))
early_effects_SE <- emmeans(sess.lm.SE, ~ stim * task, at = list(`trial_z` = min(trial_q)))
early_effects_NE <- emmeans(sess.lm.NE, ~ stim * task, at = list(`trial_z` = min(trial_q)))

pairs(early_effects_DA, by = "stim")
pairs(early_effects_SE, by = "stim")
pairs(early_effects_NE, by = "stim")

pairs(early_effects_DA, by = "task")
pairs(early_effects_SE, by = "task")
pairs(early_effects_NE, by = "task")

# Late trials
late_effects_DA <- emmeans(sess.lm.DA, ~ stim * task, at = list(`trial_z` = max(trial_q)))
late_effects_SE <- emmeans(sess.lm.SE, ~ stim * task, at = list(`trial_z` = max(trial_q)))
late_effects_NE <- emmeans(sess.lm.NE, ~ stim * task, at = list(`trial_z` = max(trial_q)))

pairs(late_effects_DA, by = "stim")
pairs(late_effects_SE, by = "stim")
pairs(late_effects_NE, by = "stim")

# Compare within task periods directly
task.merge = task.merge %>%
  group_by(idx, nt, task, stim) %>%
  mutate(period = ifelse(trial_z <= min(trial_q), "Early", 
                         ifelse(trial_z >= max(trial_q), "Late", "Middle")))


RL.lm.DA = lmer(Oz ~ stim * period +
                    (1+stim|idx),
               data = task.merge[task.merge$task == "RL" &  task.merge$nt == "DA" & task.merge$period != "Middle",],
               REML=F,
               control = lmerControl(optimizer = "bobyqa"))
RL.lm.SE = lmer(Oz ~ stim * period +
                    (1+stim|idx),
               data = task.merge[task.merge$task == "RL" &  task.merge$nt == "SE" & task.merge$period != "Middle",],
                  REML=F,
                  control = lmerControl(optimizer = "bobyqa"))
RL.lm.NE = lmer(Oz ~ stim * period +
                    (1+stim|idx),
               data = task.merge[task.merge$task == "RL" &  task.merge$nt == "NE" & task.merge$period != "Middle",],
                  REML=F,
                  control = lmerControl(optimizer = "bobyqa"))

UG.lm.DA = lmer(Oz ~ stim * period +
                    (1+stim|idx),
               data = task.merge[task.merge$task == "UG" &  task.merge$nt == "DA" & task.merge$period != "Middle",],
               REML=F,
               control = lmerControl(optimizer = "bobyqa"))
UG.lm.SE = lmer(Oz ~ stim * period +
                    (1+stim|idx),
               data = task.merge[task.merge$task == "UG" &  task.merge$nt == "SE" & task.merge$period != "Middle",],
                  REML=F,
                  control = lmerControl(optimizer = "bobyqa"))
UG.lm.NE = lmer(Oz ~ stim * period +
                    (1+stim|idx),
               data = task.merge[task.merge$task == "UG" &  task.merge$nt == "NE" & task.merge$period != "Middle",],
                  REML=F,
                  control = lmerControl(optimizer = "bobyqa"))

summary(RL.lm.DA)
summary(RL.lm.SE)
summary(RL.lm.NE)

summary(UG.lm.DA)
summary(UG.lm.SE)
summary(UG.lm.NE)

period_RL_effects_DA <- emmeans(RL.lm.DA, ~ period * stim)
period_RL_effects_SE <- emmeans(RL.lm.SE, ~ period * stim)
period_RL_effects_NE <- emmeans(RL.lm.NE, ~ period * stim)
period_UG_effects_DA <- emmeans(UG.lm.DA, ~ period * stim)
period_UG_effects_SE <- emmeans(UG.lm.SE, ~ period * stim)
period_UG_effects_NE <- emmeans(UG.lm.NE, ~ period * stim)

pairs(period_RL_effects_DA, by = "stim")
pairs(period_RL_effects_SE, by = "stim")
pairs(period_RL_effects_NE, by = "stim")
pairs(period_UG_effects_DA, by = "stim")
pairs(period_UG_effects_SE, by = "stim")
pairs(period_UG_effects_NE, by = "stim")


# Supplementary figure
# Extract slope estimates from  emtrends results
da_slopes_df = da_slopes %>% as.data.frame()
se_slopes_df = se_slopes %>% as.data.frame()
ne_slopes_df = ne_slopes %>% as.data.frame()

slope_data <- data.frame(
  nt = rep(c("DA", "5-HT", "NE"), each = 4),
  task = rep(c("RL", "UG", "RL", "UG"), 3),
  stim = rep(c("Pre-Stim", "Pre-Stim", "Post-Stim", "Post-Stim"), 3),
  slope = c(
    # DA slopes from  emtrends output
    da_slopes_df$trial_z.trend[da_slopes_df$stim == "Pre-Stim" & da_slopes_df$task == "RL"],
    da_slopes_df$trial_z.trend[da_slopes_df$stim == "Pre-Stim" & da_slopes_df$task == "UG"],
    da_slopes_df$trial_z.trend[da_slopes_df$stim == "Post-Stim" & da_slopes_df$task == "RL"],
    da_slopes_df$trial_z.trend[da_slopes_df$stim == "Post-Stim" & da_slopes_df$task == "UG"],
    se_slopes_df$trial_z.trend[se_slopes_df$stim == "Pre-Stim" & se_slopes_df$task == "RL"],
    se_slopes_df$trial_z.trend[se_slopes_df$stim == "Pre-Stim" & se_slopes_df$task == "UG"],
    se_slopes_df$trial_z.trend[se_slopes_df$stim == "Post-Stim" & se_slopes_df$task == "RL"],
    se_slopes_df$trial_z.trend[se_slopes_df$stim == "Post-Stim" & se_slopes_df$task == "UG"],
    ne_slopes_df$trial_z.trend[ne_slopes_df$stim == "Pre-Stim" & ne_slopes_df$task == "RL"],
    ne_slopes_df$trial_z.trend[ne_slopes_df$stim == "Pre-Stim" & ne_slopes_df$task == "UG"],
    ne_slopes_df$trial_z.trend[ne_slopes_df$stim == "Post-Stim" & ne_slopes_df$task == "RL"],
    ne_slopes_df$trial_z.trend[ne_slopes_df$stim == "Post-Stim" & ne_slopes_df$task == "UG"]
  ),
  lower = c(
    da_slopes_df$lower.CL[da_slopes_df$stim == "Pre-Stim" & da_slopes_df$task == "RL"],
    da_slopes_df$lower.CL[da_slopes_df$stim == "Pre-Stim" & da_slopes_df$task == "UG"],
    da_slopes_df$lower.CL[da_slopes_df$stim == "Post-Stim" & da_slopes_df$task == "RL"],
    da_slopes_df$lower.CL[da_slopes_df$stim == "Post-Stim" & da_slopes_df$task == "UG"],
    se_slopes_df$lower.CL[se_slopes_df$stim == "Pre-Stim" & se_slopes_df$task == "RL"],
    se_slopes_df$lower.CL[se_slopes_df$stim == "Pre-Stim" & se_slopes_df$task == "UG"],
    se_slopes_df$lower.CL[se_slopes_df$stim == "Post-Stim" & se_slopes_df$task == "RL"],
    se_slopes_df$lower.CL[se_slopes_df$stim == "Post-Stim" & se_slopes_df$task == "UG"],
    ne_slopes_df$lower.CL[ne_slopes_df$stim == "Pre-Stim" & ne_slopes_df$task == "RL"],
    ne_slopes_df$lower.CL[ne_slopes_df$stim == "Pre-Stim" & ne_slopes_df$task == "UG"],
    ne_slopes_df$lower.CL[ne_slopes_df$stim == "Post-Stim" & ne_slopes_df$task == "RL"],
    ne_slopes_df$lower.CL[ne_slopes_df$stim == "Post-Stim" & ne_slopes_df$task == "UG"]
  ),
  upper = c(
    da_slopes_df$upper.CL[da_slopes_df$stim == "Pre-Stim" & da_slopes_df$task == "RL"],
    da_slopes_df$upper.CL[da_slopes_df$stim == "Pre-Stim" & da_slopes_df$task == "UG"],
    da_slopes_df$upper.CL[da_slopes_df$stim == "Post-Stim" & da_slopes_df$task == "RL"],
    da_slopes_df$upper.CL[da_slopes_df$stim == "Post-Stim" & da_slopes_df$task == "UG"],
    se_slopes_df$upper.CL[se_slopes_df$stim == "Pre-Stim" & se_slopes_df$task == "RL"],
    se_slopes_df$upper.CL[se_slopes_df$stim == "Pre-Stim" & se_slopes_df$task == "UG"],
    se_slopes_df$upper.CL[se_slopes_df$stim == "Post-Stim" & se_slopes_df$task == "RL"],
    se_slopes_df$upper.CL[se_slopes_df$stim == "Post-Stim" & se_slopes_df$task == "UG"],
    ne_slopes_df$upper.CL[ne_slopes_df$stim == "Pre-Stim" & ne_slopes_df$task == "RL"],
    ne_slopes_df$upper.CL[ne_slopes_df$stim == "Pre-Stim" & ne_slopes_df$task == "UG"],
    ne_slopes_df$upper.CL[ne_slopes_df$stim == "Post-Stim" & ne_slopes_df$task == "RL"],
    ne_slopes_df$upper.CL[ne_slopes_df$stim == "Post-Stim" & ne_slopes_df$task == "UG"]
  )
  
) %>%
  mutate(nt = factor(nt, levels = c("DA","5-HT","NE")),
         stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))


# Create early vs late comparison plot
early_late_data <- task.merge %>%
  mutate(period = ifelse(trial_z <= min(trial_q), "Early", 
                         ifelse(trial_z >= max(trial_q), "Late", "Middle"))) %>%
  filter(period != "Middle") %>%
  group_by(idx, nt, task, stim, period) %>%
  summarise(m = mean(Oz))

# Panel A: Conceptual diagram showing expected vs observed patterns
conceptual_data <- data.frame(
  scenario = rep(c("If Order Effect", "If Task-Specific"), each = 4),
  task = rep(c("RL", "UG"), 2),
  time = rep(c("Early", "Late"), each = 2),
  value = c(
    # Order effect scenario: both tasks decline similarly (parallel slopes)
    1, 1,      # RL early, UG early (start same)
    0.5, 0.5,  # RL late, UG late (end same - parallel decline)
    # Task-specific scenario: divergent patterns (opposite slopes)
    0.5, 0.8,  # RL early, UG early (start different)
    0.2, 1.5   # RL late, UG late (diverge further - opposite slopes)
  )
) %>%
  expand_grid(nt = "DA") %>%
  mutate(
    time_num = ifelse(time == "Early", 1, 2),
    scenario = factor(scenario, levels = c("If Order Effect", "If Task-Specific"))
  )

# Panel A: Slope Comparison Plot
panel_a <- 
  slope_data %>%
  mutate(
    color = ifelse(stim == "Pre-Stim",
                   case_when(
                     nt == "5-HT" ~ "#9ecae1", 
                     nt == "DA" ~ "#fcbba1", 
                     nt == "NE" ~ "#c7e9c0"
                   ),
                   case_when(
                     nt == "5-HT" ~ "#08519c",
                     nt == "DA" ~ "#a50f15",
                     nt == "NE" ~ "#006d2c"
                   )),
    stim = factor(stim, levels = c("Pre-Stim","Post-Stim"))
  ) %>%
  ggplot(aes(x = task, y = slope)) +
  geom_col(position = position_dodge2(width = .5), aes(group = stim,fill = color)) +
  geom_errorbar(position = position_dodge2(width = .5), 
                linewidth = 0.75,
                #width = 0.7,
                aes(ymin = lower, ymax = upper))+
  facet_wrap(~nt, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Temporal slope \n(ΔEstimate/trial)", x = "Task", 
       fill = element_blank()) +
  theme_pubr(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color = "black",size = 15),
    panel.spacing = unit(0.2, "cm")
  ) +
  scale_fill_identity()

# Panel B & C: Early vs Late with NT-specific colors
panel_b <- 
  early_late_data %>% 
  filter(task == "UG") %>%
  mutate(    nt = recode(nt,
                         "SE" = "5-HT"),
             nt = factor(nt, levels = c("DA","5-HT","NE")),
    color = ifelse(stim == "Pre-Stim",
                   case_when(
                     nt == "5-HT" ~ "#9ecae1", 
                     nt == "DA" ~ "#fcbba1", 
                     nt == "NE" ~ "#c7e9c0"
                   ),
                   case_when(
                     nt == "5-HT" ~ "#08519c",
                     nt == "DA" ~ "#a50f15",
                     nt == "NE" ~ "#006d2c"
                   ))
  ) %>%
  ggplot(aes(x = period, y = m, color = color)) +
  geom_point(aes(group = interaction(idx, stim)), alpha = 0.3) +
  geom_line(aes(group = interaction(idx, stim)), alpha = 0.3) +
  stat_summary(aes(group = stim), size = 1 ,linewidth = 1) +
  stat_summary(geom = "line", aes(group = stim), size = 2) +
  facet_wrap(~ nt) +
  scale_color_identity() +
  labs(y = "Estimate [z]", x = "Within task period", 
       title = "Ultimatum game") +
  theme_pubr(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 14),
    panel.spacing = unit(0.2, "cm"),
    axis.text.x = element_text(angle = 15, vjust = 1, hjust=0.5)
  )
panel_c <- early_late_data %>% 
  filter(task == "RL") %>%
  mutate(
    nt = recode(nt,
                "SE" = "5-HT"),
    nt = factor(nt, levels = c("DA","5-HT","NE")),
    color = ifelse(stim == "Pre-Stim",
                   case_when(
                     nt == "5-HT" ~ "#9ecae1", 
                     nt == "DA" ~ "#fcbba1", 
                     nt == "NE" ~ "#c7e9c0"
                   ),
                   case_when(
                     nt == "5-HT" ~ "#08519c",
                     nt == "DA" ~ "#a50f15",
                     nt == "NE" ~ "#006d2c"
                   ))
  ) %>%
  ggplot(aes(x = period, y = m, color = color)) +
  geom_point(aes(group = interaction(idx, stim)), alpha = 0.3) +
  geom_line(aes(group = interaction(idx, stim)), alpha = 0.3) +
  stat_summary(aes(group = stim), size = 1,linewidth = 1) +
  stat_summary(geom = "line", aes(group = stim), size = 2) +
  facet_wrap(~ nt) +
  scale_color_identity() +
  labs(y = "Estimate [z]", x = "Within task period",title="Reversal learning") +
  theme_pubr(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 14),
    panel.spacing = unit(0.2, "cm"),
    axis.text.x = element_text(angle = 15, vjust = 1, hjust=0.5)

  )

 create_prestim_legend <- function() {
  # Dummy data for Pre-Stim legend
  legend_data_pre <- data.frame(
    nt = factor(c("DA", "5-HT", "NE"), levels = c("DA", "5-HT", "NE")),
    color = c("#fcbba1", "#9ecae1","#c7e9c0"),
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
    nt = factor(c("DA", "5-HT", "NE"), levels = c("DA", "5-HT","NE")),
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

# Combine panels with your existing legend system
supp_fig <- plot_grid(
  panel_a,
  plot_grid(panel_b, panel_c, ncol = 1),
  stacked_legends, 
  ncol = 3,
  rel_widths = c(1.2, 1.2, 0.3),
  align = "h"
)

full_supp_temporal_fig = (panel_c + panel_b) + #panel_a / (panel_c + panel_b) +
  plot_annotation(tag_levels = "a")

final_plot <- plot_grid(full_supp_temporal_fig, stacked_legends, 
                        ncol = 2, rel_widths = c(4, 1),
                        align = "h", axis = "tb")


ggsave(res_dir / "figureS2.png", 
       plot = final_plot,
       device = "png",
       width = 13,         # Width in inches
       height = 6,         # Height in inches  
       dpi = 300)          # Resolution (300 DPI for publication quality)
