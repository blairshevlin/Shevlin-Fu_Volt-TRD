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
# 2025/10/08    Blair Shevlin                           include updated NT data, OR session
# 2025/11/04    Blair Shevlin                           minor edits to finalize
# 2025/11/14    Blair Shevlin                           split into seperate figures

rm(list = ls())

library(tidyverse)
library(fs)
library(here)
library(patchwork)
library(RColorBrewer)
library(GGally)
library(purrr)
library(broom)
library(caTools)
library(rstatix)
library(scales)
library(ggrepel)
library(cowplot)
library(ggpubr)

# Paths
dir = path(here())
nt_dir = dir / "data" / "nt" / "processed"
beh_dir = dir / "data" / "behavior"
res_dir = dir / "results" # Updated
clin_dir = dir / "data" / "clinical"

# Load clinical data
cl.df = read.csv(clin_dir / "clinical-data_deid_07-10-25.csv" ) 

# Subject-level NT Data
load(nt_dir / "UG_RL_NT-Continuous_9-23-25.RData")

# Extract relevant dataframes
rl.EST.Reward = rl.EST %>%
  filter(event == "Reward") %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))

ug.EST.Offer = ug.EST %>%
  filter(event == "Offer") %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim"))) 

# IDs
ids_final = unique(ug.EST.Offer$idx)

# Specify colors
NT_colors = data.frame(id = c("DA","5-HT"),
                       color = c("#cb181d","#2171b5"))

# Load behavioral data
rl.beh = read.csv(file = beh_dir / "rl-data_deid_07-10-25.csv")

ug.beh = read.csv(file = beh_dir / "ug-data_deid_07-10-25.csv")

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

  
# HDRS ~ Beh

cl.change = 
  cl.df %>%
  filter(!session %in% c("fmri","post stim")) %>%
  mutate(sess = recode(session, 
        "pre stim" = "Baseline",
        "week 1" = "Week 1", 
        "month 1" = "Month 1",
        "month 2" = "Month 2",
        "month 3" = "Month 3",
        "month 4" = "Month 4",
        "month 5" = "Month 5",
        "month 6" = "Month 6")) %>%
  select(!c(MADRS,PVSS,SHAPS,session)) %>%
  group_by(idx) %>%
  mutate(HDRS_d = HDRS - HDRS[sess == "Baseline"])

rl.beh.change = rl.beh %>%
  filter(rt < 10) %>%
  group_by(idx,sess) %>%
  summarise(mChoice = mean(opt),
            mRT = mean(rt),
            mLogRT = mean(log(rt))) %>% 
  mutate(idx = factor(idx),
         sess = recode(sess,"fMRI" = "Baseline")) 

ug.beh.change = 
  ug.beh.means = ug.beh %>% 
  group_by(idx,sess) %>%
  filter(rt < 10) %>%
  summarise(mChoice = mean(rej==0),
            mRT = mean(rt),
            mLogRT = mean(log(rt)),
            mMood = mean(mood, na.rm = T)) %>%   
  mutate(idx = factor(idx),
         sess = recode(sess,"fMRI" = "Baseline")) %>%
  as.data.frame()


rl.beh.cl = merge(rl.beh.change, cl.change) 
ug.beh.cl = merge(ug.beh.change, cl.change) 

rl.beh.cl.l = 
  rl.beh.cl %>% 
  pivot_longer(cols =c("mChoice","mRT","mLogRT"), names_to = "beh")

ug.beh.cl.l = 
  ug.beh.cl %>% 
  pivot_longer(cols =c("mChoice","mRT","mLogRT", "mMood"), names_to = "beh")

# Run correlations and tidy the results

# Change scores
rl.beh.cl.l %>%
  filter(!sess == "Baseline") %>%
  group_by(beh,sess) %>%
  do(tidy(cor.test(.$value, .$HDRS_d, method = "spearman"))) %>%
  # Filter for significant results (p < 0.05)
  filter(p.value < .05)

  rl.beh.cl.l %>%
  filter(!sess == "Baseline") %>%
  group_by(beh,sess) %>%
  do(tidy(cor.test(.$value, .$HDRS_d, method = "spearman"))) 

ug.beh.cl.l %>%
  filter(!sess == "Baseline") %>%
  group_by(beh,sess) %>%
  do(tidy(cor.test(.$value, .$HDRS_d, method = "spearman"))) %>% as.data.frame() %>%
  # Filter for significant results (p < 0.05)
  filter(p.value < .05)

# Raw scores
rl.beh.cl.l %>%
  filter(!sess == "Baseline") %>%
  group_by(beh,sess) %>%
  do(tidy(cor.test(.$value, .$HDRS, method = "spearman"))) %>%
  # Filter for significant results (p < 0.05)
  filter(p.value < .05)

ug.beh.cl.l %>%
  filter(!sess == "Baseline") %>%
  group_by(beh,sess) %>%
  do(tidy(cor.test(.$value, .$HDRS, method = "spearman"))) %>%
  # Filter for significant results (p < 0.05)
  filter(p.value < .05)


# Plot change HDRS with beh

HDRS.rl.choice = 
  rl.beh.cl.l %>%
  filter(sess %in% c("Month 1", "Month 3", "Month 6"), beh == "mChoice") %>%
  ggplot(aes(x = HDRS_d, y = value)) +
  theme_pubr(base_size = 16) +
  facet_wrap(~sess, nrow = 1) +
  stat_smooth(method = "lm", 
              alpha = .25,
              linetype = "solid",
              linewidth = 2.1,
              color = "purple", 
              fill = "purple" ) +
  geom_point(size = 3,stroke = 1.5) +
  labs(title = "P(optimal) [RL]",
       x = "HDRS Change from Baseline",
       y = "Change from Baseline",
       shape = element_blank(),
       color = element_blank()) +
  scale_x_continuous(breaks=c(-20,-10,0)) +
  theme(
    panel.spacing = unit(1, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = NA, fill = NA)) +
  stat_cor(method="spearman",aes(label =paste0(cut(..p.., 
                                                   breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                   labels = c("'****'", "'***'", "'**'", "'*'", "'ns'"))
  )),label.x = -15, size = c(7,7,7), label.y = 1)  

HDRS.rl.rt = 
  rl.beh.cl.l %>%
  filter(sess %in% c("Month 1", "Month 3", "Month 6"), beh == "mLogRT") %>%
  ggplot(aes(x = HDRS_d, y = value)) +
  theme_pubr(base_size = 16) +
  facet_wrap(~sess, nrow = 1) +
  stat_smooth(method = "lm", 
              alpha = .25,
              linetype = "solid",
              linewidth = 2.1,
              color = "purple", 
              fill = "purple" ) +
  geom_point(size = 3,stroke = 1.5) +
  labs(title = "Response time [RL]",
       x = "HDRS Change from Baseline",
       y = "Change from Baseline",
       shape = element_blank(),
       color = element_blank()) +
  scale_x_continuous(breaks=c(-20,-10,0)) +
  theme(
    panel.spacing = unit(1, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = NA, fill = NA)) +
  stat_cor(method="spearman",aes(label =paste0(cut(..p.., 
                                                   breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                   labels = c("'****'", "'***'", "'**'", "'*'", "'ns'"))
  )),label.x = -15, size = c(7,7,7), label.y = 1)  

HDRS.ug.mood = ug.beh.cl.l %>%
  filter(sess %in% c("Month 1", "Month 3", "Month 6"), beh == "mMood") %>%
  ggplot(aes(x = HDRS_d, y = value)) +
  theme_pubr(base_size = 16) +
  facet_wrap(~sess, nrow = 1) +
  stat_smooth(method = "lm", 
              alpha = .25,
              linetype = "solid",
              linewidth = 2.1,
              color = "purple", 
              fill = "purple" ) +
  geom_point(size = 3,stroke = 1.5) +
  labs(title = "Mood [UG]",
       x = "HDRS Change from Baseline",
       y = "Change from Baseline",
       shape = element_blank(),
       color = element_blank()) +
  scale_x_continuous(breaks=c(-20,-10,0)) +
  theme(
    panel.spacing = unit(1, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = NA, fill = NA)) +
  stat_cor(method="spearman",aes(label =paste0(cut(..p.., 
                                                   breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                   labels = c("'****'", "'***'", "'**'", "'*'", "'ns'"))
  )),label.x = -15, size = c(7,12,7), label.y = 85)

HDRS.ug.choice = 
  ug.beh.cl.l %>%
  filter(sess %in% c("Month 1", "Month 3", "Month 6"), beh == "mChoice") %>%
  ggplot(aes(x = HDRS_d, y = value)) +
  theme_pubr(base_size = 16) +
  facet_wrap(~sess, nrow = 1) +
  stat_smooth(method = "lm", 
              alpha = .25,
              linetype = "solid",
              linewidth = 2.1,
              color = "purple", 
              fill = "purple" ) +
  geom_point(size = 3,stroke = 1.5) +
  labs(title = "P(accept) [UG]",
       x = "HDRS Change from Baseline",
       y = "Change from Baseline",
       shape = element_blank(),
       color = element_blank()) +
  scale_x_continuous(breaks=c(-20,-10,0)) +
  theme(
    panel.spacing = unit(1, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = NA, fill = NA)) +
  stat_cor(method="spearman",aes(label =paste0(cut(..p.., 
                                                   breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                   labels = c("'****'", "'***'", "'**'", "'*'", "'ns'"))
  )),label.x = -15, size = c(12,12,7), label.y = 1.25)  

HDRS.ug.rt = 
  ug.beh.cl.l %>%
  filter(sess %in% c("Month 1", "Month 3", "Month 6"), beh == "mLogRT") %>%
  ggplot(aes(x = HDRS_d, y = value)) +
  theme_pubr(base_size = 16) +
  facet_wrap(~sess, nrow = 1) +
  stat_smooth(method = "lm", 
              alpha = .25,
              linetype = "solid",
              linewidth = 2.1,
              color = "purple", 
              fill = "purple" ) +
  geom_point(size = 3,stroke = 1.5) +
  labs(title = "Response time [UG]",
       x = "HDRS Change from Baseline",
       y = "Change from Baseline",
       shape = element_blank(),
       color = element_blank()) +
  scale_x_continuous(breaks=c(-20,-10,0)) +
  theme(
    panel.spacing = unit(1, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = NA, fill = NA)) +
  stat_cor(method="spearman",aes(label =paste0(cut(..p.., 
                                                   breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                   labels = c("'****'", "'***'", "'**'", "'*'", "'ns'"))
  )),label.x = -15, size = c(7,7,7), label.y = 1)    


figure4 = ( HDRS.rl.rt + HDRS.rl.choice)/(HDRS.ug.rt + HDRS.ug.choice + HDRS.ug.mood) +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = "collect",nrow = 2) &
  theme(legend.position = 'right') 


ggsave(res_dir / "figureS4_raw.png", 
       plot = figure4,
       device = "png",
       width = 15,          # Width in inches
       height = 8,         # Height in inches  
       dpi = 300)           # Resolution (300 DPI for publication quality)

