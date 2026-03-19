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
# 2026/03/19    Blair Shevlin                           refactored to read from source data CSV

# Reads source data from data/figures/figureSupplement3_source_data.csv

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

source_data <- read.csv(here::here("data/figures/figureSupplement3_source_data.csv"), stringsAsFactors = FALSE)

# Re-apply factor orderings
sessions = c("Baseline", "Pre-Stim", "Post-Stim", "DBS","Month 1","Month 3","Month 6")
source_data$sess <- factor(source_data$sess, levels = sessions)
source_data$idx  <- factor(source_data$idx)

# Split panels
rl.beh.means.plot <- source_data %>% filter(panel == "RL_choice")
rl.beh.logrt.plot <- source_data %>% filter(panel == "RL_RT")
ug.beh.means.plot <- source_data %>% filter(panel == "UG_choice")
ug.beh.logrt.plot <- source_data %>% filter(panel == "UG_RT")
ug.mood.means.plot <- source_data %>% filter(panel == "UG_mood")

rl.c.sess =
  rl.beh.means.plot %>%
  mutate(sess= factor(sess, levels = sessions)) %>%
  ggplot(aes(x = sess, y = mChoice))+
  theme_pubr(base_size = 14) +
  geom_vline(xintercept = "DBS", linewidth = 2) +
    geom_boxplot(data = rl.beh.means.plot,
               linewidth = 1.1,outlier.alpha = 0,
               show.legend = F) +
  geom_point(data = rl.beh.means.plot,
             size=2,color = "purple",alpha=.5,
             position = position_dodge2(width = .3),
             stroke = 1.75,
             aes( group = idx
             ),shape = 18
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x =  element_blank(), y = "p(optimal)",
       title = "Choice [RL]",
       color = element_blank(),
       shape = element_blank()) +
  scale_shape_manual(values = c(16,4)) +
  stat_compare_means(
    label="p.signif",
    comparisons = list(
      c("Baseline", "Post-Stim"),
      c("Pre-Stim", "Post-Stim") ),
    method = "wilcox",
    paired = TRUE,
    bracket.size = 0.5,
    tip.length = 0.02,
    step.increase = 0.08,label.x = 2.5
  )

ug.c.sess =
  ug.beh.means.plot %>%
  ggplot(aes(x = sess, y = mChoice))+
  theme_pubr(base_size = 14) +
  geom_vline(xintercept = "DBS", linewidth = 2) +
  geom_boxplot(data = ug.beh.means.plot,
               linewidth = 1.1,outlier.alpha = 0,
               show.legend = F) +
  geom_point(data = ug.beh.means.plot,
             size=2,color = "purple",alpha=.5,
             position = position_dodge2(width = .3),
             stroke = 1.75,aes( group = idx
             ),shape = 18
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x =  element_blank(), y = "p(accept)",
       title = "Choice [UG]",
       color = element_blank(),
       shape = element_blank()) +
  scale_shape_manual(values = c(16,4)) +
  stat_compare_means(
    label="p.signif",
    comparisons = list(
      c("Baseline", "Post-Stim"),
      c("Pre-Stim", "Post-Stim") ),
    method = "wilcox",    paired = TRUE,
    bracket.size = 0.5,
    tip.length = 0.02,
    step.increase = 0.08,label.x = 2.5
  )

ug.logrt.sess =
  ug.beh.logrt.plot %>%
  ggplot(aes(x = sess, y = mLogRT))+
  theme_pubr(base_size = 14) +
  geom_vline(xintercept = "DBS", linewidth = 2) +
  geom_boxplot(data = ug.beh.logrt.plot,
               linewidth = 1.1,outlier.alpha = 0,
               show.legend = F) +
  geom_point(data = ug.beh.logrt.plot,
             size=2,color = "purple",alpha=.5,
             position = position_dodge2(width = .3),
             stroke = 1.75,aes( group = idx
             ),shape = 18
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x =  element_blank(), y = "log sec",
       color = element_blank(),title = "Response time [UG]",
       shape = element_blank()) +
  scale_shape_manual(values = c(16,4)) +
  stat_compare_means(
    label="p.signif",
    comparisons = list(
      c("Baseline", "Post-Stim"),
      c("Pre-Stim", "Post-Stim") ),
    method = "wilcox",    paired = TRUE,

    bracket.size = 0.5,
    tip.length = 0.02,
    step.increase = 0.08,label.x =2.5
  )

rl.logrt.sess =
  rl.beh.logrt.plot %>%
  ggplot(aes(x = sess, y = mLogRT))+
  theme_pubr(base_size = 14) +
  geom_vline(xintercept = "DBS", linewidth = 2) +
  geom_boxplot(data = rl.beh.logrt.plot,
               linewidth = 1.1,outlier.alpha = 0,
               show.legend = F) +
  geom_point(data = rl.beh.logrt.plot,
             size=2,color = "purple",alpha=.5,
             position = position_dodge2(width = .3),
             stroke = 1.75,aes( group = idx
             ),shape = 18
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x =  element_blank(), y = "log sec",
       title = "Response time [RL]",
       color = element_blank(),
       shape = element_blank()) +
  scale_shape_manual(values = c(16,4)) +
  stat_compare_means(
    label="p.signif",
    comparisons = list(
      c("Baseline", "Post-Stim"),
      c("Pre-Stim", "Post-Stim") ),
    method = "wilcox",    paired = TRUE,

    bracket.size = 0.5,
    tip.length = 0.02,
    step.increase = 0.08,label.x = 2.5
  )

ug.mood.sess =
  ug.mood.means.plot %>%
  ggplot(aes(x = sess, y = mMood))+
  theme_pubr(base_size = 14) +
  geom_vline(xintercept = "DBS", linewidth = 2) +
  geom_boxplot(data = ug.mood.means.plot,
               linewidth = 1.1,outlier.alpha = 0,
               show.legend = F) +
  geom_point(data = ug.mood.means.plot,
             size=2,color = "purple",alpha=.5,
             position = position_dodge2(width = .3),
             stroke = 1.75,aes( group = idx
             ),shape = 18
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x =  element_blank(), y = "rating",
       title = "Mood [UG]",
       color = element_blank(),
       shape = element_blank()) +
  scale_shape_manual(values = c(16,4)) +
  stat_compare_means(
    label= "p.signif",
    comparisons = list(
      c("Baseline", "Post-Stim"),
      c("Pre-Stim", "Post-Stim") ),
      method = "wilcox",    paired = TRUE,

      bracket.size = 0.5,
      tip.length = 0.02,
      step.increase = 0.08
    )


figure3 = rl.logrt.sess+rl.c.sess+ ug.logrt.sess+ug.c.sess + ug.mood.sess  +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect",nrow = 1) & theme(legend.position = "none")


ggsave(here::here("results/from_source_data/figureS3_raw.png"),
       plot = figure3,
       device = "png",
       width = 15,          # Width in inches
       height = 4,         # Height in inches
       dpi = 300)           # Resolution (300 DPI for publication quality)
