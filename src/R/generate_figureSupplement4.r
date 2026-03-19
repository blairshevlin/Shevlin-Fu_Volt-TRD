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

# Reads source data from data/figures/figureSupplement4_source_data.csv

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

source_data <- read.csv(here::here("data/figures/figureSupplement4_source_data.csv"), stringsAsFactors = FALSE)

# Re-apply factor orderings
source_data$sess <- factor(source_data$sess, levels = c("Month 1", "Month 3", "Month 6"))
source_data$idx  <- factor(source_data$idx)

# Split panels
panel_A <- source_data %>% filter(panel == "A")  # RL choice vs HDRS change
panel_B <- source_data %>% filter(panel == "B")  # RL RT vs HDRS change
panel_C <- source_data %>% filter(panel == "C")  # UG mood vs HDRS change
panel_D <- source_data %>% filter(panel == "D")  # UG choice vs HDRS change
panel_E <- source_data %>% filter(panel == "E")  # UG RT vs HDRS change

HDRS.rl.choice =
  panel_A %>%
  ggplot(aes(x = HDRS_change, y = beh_change)) +
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
  panel_B %>%
  ggplot(aes(x = HDRS_change, y = beh_change)) +
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

HDRS.ug.mood = panel_C %>%
  ggplot(aes(x = HDRS_change, y = beh_change)) +
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
  panel_D %>%
  ggplot(aes(x = HDRS_change, y = beh_change)) +
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
  panel_E %>%
  ggplot(aes(x = HDRS_change, y = beh_change)) +
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


ggsave(here::here("results/from_source_data/figureS4_raw.png"),
       plot = figure4,
       device = "png",
       width = 15,          # Width in inches
       height = 8,         # Height in inches
       dpi = 300)           # Resolution (300 DPI for publication quality)
