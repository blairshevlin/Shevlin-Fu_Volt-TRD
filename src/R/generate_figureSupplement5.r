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
# 2025/11/05    Blair Shevlin                          finalized code for rebuttal

rm(list = ls())

library(tidyverse)
library(fs)
library(here)
library(patchwork)
library(ggpubr)
library(RColorBrewer)
library(lmerTest)
library(cowplot)

# Paths
dir = path(here())
nt_dir = dir / "data" / "nt" / "processed"
beh_dir = dir / "data" / "behavior"
clin_dir = dir / "data" / "clinical"
res_dir = dir / "results" # Updated

# Load behavioral data
rl.beh = read.csv(file = beh_dir / "rl-data_deid_07-10-25.csv")

ug.beh = read.csv(file = beh_dir / "ug-data_deid_07-10-25.csv")

# Specify colors
NT_colors = data.frame(id = c("DA","5-HT"),
                       color = c("#cb181d","#2171b5"))

# Subject-level NT Data
load(nt_dir / "UG_RL_NT-Continuous_9-23-25.RData")

# Generate summary data - focusing only on reward events (RL) and offer events (UG)
# This time using data that was z-scored using both tasks

tasks.EST.alt = tasks.EST.alt %>% 
  group_by(task,nt,event,idx) %>%
  mutate(trial_z = scale(trial)[,1]) %>%
  ungroup()


rl.EST.Shared= 
  tasks.EST.alt %>% filter(task == "RL") %>%
  filter(event == "Shared") %>%
  mutate(trial_z = scale(trial)[,1])

ug.EST.Shared = 
  tasks.EST.alt %>% filter(task == "UG") %>%
  filter(event == "Shared") %>%
  mutate(trial_z = scale(trial)[,1])

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


# Re-parameterize
rl.EST.pr = 
    rl.EST.Shared %>%
    group_by(idx,stim,block,event,nt) %>%
    mutate(prev_rew_raw = ifelse(trial_within_block==1,NA,lag(outcome)),
                        prev_rew = ifelse(is.na(prev_rew_raw),NA,
                                        ifelse(prev_rew_raw == 1,
                                                ifelse(cond %in% c("Mixed", "Positive"),10,0),
                                                ifelse(cond %in% c("Mixed", "Negative"),-10,0)
                                        )
                        ),
                 stim = factor(stim, levels = c("Pre-Stim","Post-Stim"))                   
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
                                 ifelse(prev_offer_raw > offer, "-PE","+PE")),
            stim = factor(stim, levels = c("Pre-Stim","Post-Stim"))
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
rl.DA.lm = lmer(   data = tasks.EST.alt[tasks.EST.alt$task == "RL" & tasks.EST.alt$event == "Shared" & tasks.EST.alt$nt == "DA",],
                formula = Oz ~ stim * trial_z + (1+stim|idx),
                control = lmerControl(optimizer="bobyqa"))
summary(rl.DA.lm)

rl.SE.lm = lmer(   data = tasks.EST.alt[tasks.EST.alt$task == "RL" & tasks.EST.alt$event == "Shared" & tasks.EST.alt$nt == "SE",],
                formula = Oz ~ stim * trial_z + (1+stim|idx),
                control = lmerControl(optimizer="bobyqa"))
summary(rl.SE.lm)

rl.NE.lm = lmer(   data = tasks.EST.alt[tasks.EST.alt$task == "RL" & tasks.EST.alt$event == "Shared" & tasks.EST.alt$nt == "NE",],
                formula = Oz ~ stim * trial_z  + (1+stim|idx),
                control = lmerControl(optimizer="bobyqa"))
summary(rl.NE.lm)

ug.DA.lm = lmer(   data = tasks.EST.alt[tasks.EST.alt$task == "UG" &tasks.EST.alt$event == "Shared" & tasks.EST.alt$nt == "DA",],
                formula = Oz ~ stim * trial_z + (1+stim|idx),
                control = lmerControl(optimizer="bobyqa"))
summary(ug.DA.lm)
ug.SE.lm = lmer(   data = tasks.EST.alt[tasks.EST.alt$task == "UG" &tasks.EST.alt$event == "Shared" & tasks.EST.alt$nt == "SE",],
                formula = Oz ~ stim * trial_z + (1+stim|idx),
                control = lmerControl(optimizer="bobyqa"))
summary(ug.SE.lm)
ug.NE.lm = lmer(   data = tasks.EST.alt[tasks.EST.alt$task == "UG" & tasks.EST.alt$event == "Shared" & tasks.EST.alt$nt == "NE",],
                formula = Oz ~ stim * trial_z + (1+stim|idx),
                control = lmerControl(optimizer="bobyqa"))
summary(ug.NE.lm)


DA.lm = lmer(   data = tasks.EST.alt[tasks.EST.alt$event == "Shared" & tasks.EST.alt$nt == "DA",],
                formula = Oz ~ stim * task + (1+stim|idx),
                control = lmerControl(optimizer="bobyqa"))
summary(DA.lm)
SE.lm = lmer(   data = tasks.EST.alt[tasks.EST.alt$event == "Shared" & tasks.EST.alt$nt == "SE",],
                formula = Oz ~ stim * task + (1+stim|idx),
                control = lmerControl(optimizer="bobyqa"))
summary(SE.lm)
NE.lm = lmer(   data = tasks.EST.alt[tasks.EST.alt$event == "Shared" & tasks.EST.alt$nt == "NE",],
                formula = Oz ~ stim * task + (1+stim|idx),
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

# Recreate Figure 2 using new approach

ug.trial = ug.EST.Shared %>%
  filter(nt %in% c("DA","SE")) %>%
  pivot_longer(cols = c("Oz","Rz","Pz","Mz"),
               names_to = "nt_metric", values_to = "nt_val") %>%
  mutate(offer_bin = ifelse(offer < 4, "$1-3",
                            ifelse(offer < 7, "$4-6","$7-9")),
         offer_bin = factor(offer_bin,
                            levels = c("$1-3","$4-6","$7-9"))) %>%
  group_by(idx,stim,nt,nt_metric,trial,offer,offer_bin) %>%
  summarise(mTrial = mean(nt_val)) %>% ungroup() %>%
  mutate(offerz = scale(offer)[,1])

rl.trial = rl.EST.Shared %>%
  filter(nt %in% c("DA","SE")) %>%
  pivot_longer(cols = c("Oz","Rz","Pz","Mz"),
               names_to = "nt_metric", values_to = "nt_val") %>%
  mutate(block_type = factor(cond, levels = c("Mixed","Negative","Positive"))
  ) %>%
  group_by(idx,stim,nt,nt_metric,trial,block_type,trial_within_block) %>%
  summarise(mTrial = mean(nt_val)) %>% ungroup()

fig.ug.trial = 
  ug.trial %>%
  filter(nt_metric == "Oz") %>%
  group_by(idx,stim,nt,trial) %>%
  summarise(m = mean(mTrial)) %>%
  mutate(nt = factor(nt, levels = c("DA","SE"), labels = c("DA","5-HT")),
         stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         color = ifelse(stim == "Pre-Stim",
                        ifelse(nt == "5-HT","#9ecae1", "#fcbba1"),
                        ifelse(nt == "5-HT","#08519c", "#a50f15"))
  ) %>%
  ggplot(   aes(x = trial, y = m, fill = color, group = interaction(nt,stim))) +
  theme_pubr(base_size = 14) +
  facet_wrap(~nt) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  stat_summary( aes(color = color),geom = "line", size = 1, linewidth = 1.5) +
  stat_summary(
    geom = "ribbon", alpha = 0.35,show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_identity() +
  # Common styling
  labs(y = "Estimate [z]", x = "Trial", title = "Ultimatum game") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        strip.text = element_text(color = "black"),
        panel.spacing = unit(0.2, "cm")) +
  coord_cartesian(ylim=c(-10,10))+
  guides(color = "none",
         fill = "none",
  ) +
  theme(legend.key = element_rect(color = NA)) 

fig.rl.trial = 
rl.trial %>%
  filter(nt_metric == "Oz") %>%
  group_by(idx,stim,nt,trial) %>%
  summarise(m = mean(mTrial)) %>%
  mutate(nt = factor(nt, levels = c("DA","SE","NE"), labels = c("DA","5-HT","NE")),
         stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         color = ifelse(stim == "Pre-Stim",
                        ifelse(nt == "5-HT","#9ecae1", "#fcbba1"),
                        ifelse(nt == "5-HT","#08519c","#a50f15"))
  ) %>%
  ggplot(   aes(x = trial, y = m, fill = color, group = interaction(nt,stim))) +
  theme_pubr(base_size = 14) +
  facet_wrap(~nt) +
  scale_color_identity() +
  scale_fill_identity() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  stat_summary( aes(color = color),geom = "line", size = 1, linewidth = 1.5) +
  stat_summary(
    geom = "ribbon", alpha = 0.35,
    show.legend = FALSE
  ) +
  # Common styling
  labs(y = "Estimate [z]", x = "Trial", title = "Reversal learning") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        strip.text = element_text(color = "black"),
        panel.spacing = unit(0.2, "cm")) +
  coord_cartesian(ylim=c(-10,10))+
  guides(color = "none",
         fill = "none",
  ) +
  theme(legend.key = element_rect(color = NA)) 

# Create participant avg plots
rl.avg =   rl.trial %>%
  filter(nt_metric == "Oz") %>%
  group_by(idx,stim,nt) %>%
  summarise(m = mean(mTrial)) %>%
  mutate(nt = factor(nt, levels = c("DA","SE"), labels = c("DA","5-HT")),
         stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         color = ifelse(stim == "Pre-Stim",
                        ifelse(nt == "5-HT","#9ecae1","#fcbba1"),
                        ifelse(nt == "5-HT","#08519c", "#a50f15"))
  ) %>% as.data.frame()

ug.avg =   ug.trial %>%
  filter(nt_metric == "Oz") %>%
  group_by(idx,stim,nt) %>%
  summarise(m = mean(mTrial)) %>%
  mutate(nt = factor(nt, levels = c("DA","SE"), labels = c("DA","5-HT")),
         stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         color = ifelse(stim == "Pre-Stim",
                        ifelse(nt == "5-HT","#9ecae1","#fcbba1"),
                        ifelse(nt == "5-HT","#08519c", "#a50f15"))
  ) %>% as.data.frame()

offset_data_RL <- rl.avg %>%
  mutate(
    # For points: Baseline to right, Month 6 to left
    point_x_pos = ifelse(stim == "Pre-Stim", 1 + 0.1, 2 - 0.1),
    # For boxplots: Baseline to left, Month 6 to right (opposite of points)
    box_x_pos = ifelse(stim == "Pre-Stim", 1 - 0.1, 2 + 0.1),
    # Convert sess_fig to numeric for proper plotting
    sess_num = ifelse(stim == "Pre-Stim", 1, 2)
  )
offset_data_UG <- ug.avg %>%
  mutate(
    # For points: Baseline to right, Month 6 to left
    point_x_pos = ifelse(stim == "Pre-Stim", 1 + 0.1, 2 - 0.1),
    # For boxplots: Baseline to left, Month 6 to right (opposite of points)
    box_x_pos = ifelse(stim == "Pre-Stim", 1 - 0.1, 2 + 0.1),
    # Convert sess_fig to numeric for proper plotting
    sess_num = ifelse(stim == "Pre-Stim", 1, 2)
  )

fig.rl.avg = 
  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  geom_line(data = offset_data_RL,
    aes(x = point_x_pos,y = m, group = idx),
    position = position_dodge2(width = .2),
    size = 1.5, alpha = .15
  ) +
  geom_point(data = offset_data_RL, 
    aes(x = point_x_pos, y = m, group = idx, color = color),
    position = position_dodge2(width = .2), 
    shape = 1,
    stroke = 2,
    alpha = 0.5,
    size = 3
  ) +
  stat_summary(geom = "line", data = offset_data_RL,
    aes(x = box_x_pos, y = m, group = nt),
    linewidth = 1.1)+
  stat_summary(data = offset_data_RL,
    aes(x = box_x_pos,y = m, group = stim, color = color),
    size = 1.1,
    #shape = 15,
    linewidth = 1.1) +
  facet_wrap(~ nt) +
  scale_color_identity() +
  labs(y = "Estimate [z]", x = element_blank(), 
       title = "Reversal learning") +
  theme_pubr(base_size = 14) +
    scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Pre-Stim", "Post-Stim"),
    limits = c(0.5, 2.5)
  ) +
  coord_cartesian(ylim=c(-10,10))+
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 14),
    panel.spacing = unit(0.2, "cm"),
    axis.text.x = element_text(angle = 15, vjust = 0.95, hjust=0.5)
  )

fig.ug.avg = 
  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  geom_line(data = offset_data_UG,
    aes(x = point_x_pos,y = m, group = idx),
    position = position_dodge2(width = .2), 
    size = 1.5, 
    alpha = .15) +
  geom_point(data = offset_data_UG, 
    aes(x = point_x_pos, y = m, group = idx, color = color),
    position = position_dodge2(width = .2), 
    shape = 1,
    stroke = 2,
    alpha = 0.5,
    size = 3) +
  stat_summary(geom = "line", data = offset_data_UG,
    aes(x = box_x_pos, y = m, group = nt),
    linewidth = 1.1
  )+
  stat_summary(data = offset_data_UG,aes(x = box_x_pos,y = m, group = stim, color = color), 
  size = 1.1 ,linewidth = 1.1, 
  #shape = 15
  ) +
  facet_wrap(~ nt) +
  scale_color_identity() +
   coord_cartesian(ylim=c(-10,10))+
  labs(y = "Estimate [z]", x = element_blank(), 
       title = "Ultimatum game") +
  theme_pubr(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 14),
    panel.spacing = unit(0.2, "cm"),
    axis.text.x = element_text(angle = 15, vjust = 0.95, hjust=0.5)
  )+
    scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Pre-Stim", "Post-Stim"),
    limits = c(0.5, 2.5)
  )

# Task feature data

offset_data_UG_task <- ug.trial %>%
  filter(nt_metric == "Oz") %>%
  group_by(idx,stim,nt,offer_bin) %>%
  summarise(m = mean(mTrial)) %>%
  mutate(nt = factor(nt, levels = c("DA","SE"), labels = c("DA","5-HT")),
         stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         color = ifelse(stim == "Pre-Stim",
                        ifelse(nt == "5-HT","#9ecae1", "#fcbba1"),
                        ifelse(nt == "5-HT","#08519c", "#a50f15")),
         # Create numeric positions for offsetting
         offer_num = as.numeric(offer_bin),
         # Points offset inward, group means offset outward
        # point_x_pos = ifelse(stim == "Pre-Stim", offer_num - 0.1, offer_num + 0.1),
         #group_x_pos = ifelse(stim == "Pre-Stim", offer_num - 0.3, offer_num + 0.3)
         point_x_pos = ifelse(stim == "Pre-Stim", offer_num - 0.15, offer_num + 0.15),
         group_x_pos = ifelse(stim == "Pre-Stim", offer_num - 0.35, offer_num + 0.35)
  )

  offset_data_RL_task <- rl.trial %>%
  filter(nt_metric == "Oz") %>%
  group_by(idx,stim,nt,block_type) %>%
  summarise(m = mean(mTrial)) %>%
  mutate(nt = factor(nt, levels = c("DA","SE"), labels = c("DA","5-HT")),
         stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         block_type = factor(block_type,
                             levels = c("Positive","Mixed","Negative"),
                             labels = c("Reward","Mixed","Punishment")),
         color = ifelse(stim == "Pre-Stim",
                        ifelse(nt == "5-HT","#9ecae1","#fcbba1"),
                        ifelse(nt == "5-HT","#08519c", "#a50f15")),
         # Create numeric positions for offsetting
         block_num = as.numeric(block_type),
         # Points offset inward, group means offset outward
         point_x_pos = ifelse(stim == "Pre-Stim", block_num - 0.15, block_num + 0.15),
         group_x_pos = ifelse(stim == "Pre-Stim", block_num - 0.35, block_num + 0.35) #was 0.3
  )


fig.rl.task = 
  ggplot() +
  theme_pubr(base_size = 14) +
  facet_wrap(~nt) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  # Individual points with inward offset
  geom_point(data = offset_data_RL_task,
             aes(x = point_x_pos, y = m, color = color),
             alpha = .75, size =  1.5, stroke = 1.5, shape = 1) + #size was 2
  # Group means (diamonds) with outward offset
  stat_summary(data = offset_data_RL_task,
               aes(x = group_x_pos, y = m, color = color, group = stim),
               size = 1.25, linewidth = 1.25,
               position = position_identity()) +
  scale_color_identity() +
  scale_fill_identity() +
  # Use original block_type labels for x-axis
  scale_x_continuous(
    breaks = c(1, 2, 3),
    labels = c("Reward", "Mixed", "Punishment"),
    limits = c(0.5, 3.5)
  ) +
  labs(y = "Estimate [z]", x = "Block type", title = "Reversal learning") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        strip.text = element_text(color = "black"),
        panel.spacing = unit(0.2, "cm"),
        axis.text.x = element_text(angle = 15, vjust = 0.9, hjust=0.6)) +
   coord_cartesian(ylim=c(-10,10))+
  guides(color = "none", fill = "none") +
  theme(legend.key = element_rect(color = NA))

fig.ug.task = 
  ggplot() +
  theme_pubr(base_size = 14) +
  facet_wrap(~nt) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  geom_point(data = offset_data_UG_task, 
    aes(x = point_x_pos, y = m, group = idx, color = color),
    shape = 1,
    stroke = 1.75,
    alpha = 0.75,
    size = 1.5) +
  stat_summary(data = offset_data_UG_task,
  aes(x = group_x_pos,y = m, group = stim, color = color), 
  size = 1.25 ,linewidth = 1.25, 
  #shape = 15
  ) +
  scale_color_identity() +
  scale_fill_identity() +
  # Use original offer_bin labels for x-axis
  scale_x_continuous(
    breaks = c(1, 2, 3),
    labels = c("$1-3", "$4-6", "$7-9"),
    limits = c(0.5, 3.5)
  ) +
  labs(y = "Estimate [z]", x = "Offer", title = "Ultimatum game") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        strip.text = element_text(color = "black"),
        panel.spacing = unit(0.2, "cm")) +
   coord_cartesian(ylim=c(-10,10))+
  guides(color = "none", fill = "none") +
  theme(legend.key = element_rect(color = NA))

create_prestim_legend <- function() {
  # Dummy data for Pre-Stim legend
  legend_data_pre <- data.frame(
    nt = factor(c("DA", "5-HT"), levels = c("DA", "5-HT")),
    color = c("#fcbba1", "#9ecae1"),
    x = 1,
    y = 1
  )
  
  ggplot(legend_data_pre, aes(x = x, y = y, color = nt)) +
    geom_point(size = 0) +  # Invisible points
    geom_line(size = 2) +   # This creates the line legend
    scale_color_manual(
      name = "Pre-Stim",
      values = c("DA" = "#fcbba1", "5-HT" = "#9ecae1")
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
    nt = factor(c("DA", "5-HT"), levels = c("DA", "5-HT")),
    color = c("#a50f15", "#08519c"),
    x = 1,
    y = 1
  )
  
  ggplot(legend_data_post, aes(x = x, y = y, color = nt)) +
    geom_point(size = 0) +  # Invisible points
    geom_line(size = 2) +   # This creates the line legend
    scale_color_manual(
      name = "Post-Stim",
      values = c("DA" = "#a50f15", "5-HT" = "#08519c")
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

 main_plot <-(fig.rl.trial | fig.ug.trial) / 
             (fig.rl.avg | fig.ug.avg) / 
             (fig.rl.task | fig.ug.task) + 
  plot_annotation(tag_levels = "a")

final_plot <- plot_grid(main_plot, stacked_legends, 
                        ncol = 2, rel_widths = c(4, 1),
                        align = "h", axis = "tb")

ggsave(res_dir / "figureS5.png", 
       plot = final_plot,
       device = "png",
       width = 12,          # Width in inches
       height = 10,         # Height in inches  
       dpi = 300)           # Resolution (300 DPI for publication quality)


