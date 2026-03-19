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
# 2026/03/19   Blair Shevlin                          refactored to read from source data CSV

# Reads source data from data/figures/figureSupplement2_source_data.csv

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

# Read source data
source_data <- read.csv(here::here("data/figures/figureSupplement2_source_data.csv"), stringsAsFactors = FALSE)

# Re-apply factor orderings
source_data$stim <- factor(source_data$stim, levels = c("Pre-Stim", "Post-Stim"))
source_data$nt   <- factor(source_data$nt,   levels = c("DA", "5-HT", "NE"))

# Panel A: temporal slope estimates (group-level, model-derived)
slope_data <- source_data %>%
  filter(panel == "A") %>%
  mutate(
    color = ifelse(stim == "Pre-Stim",
                   case_when(
                     nt == "5-HT" ~ "#9ecae1",
                     nt == "DA"   ~ "#fcbba1",
                     nt == "NE"   ~ "#c7e9c0"
                   ),
                   case_when(
                     nt == "5-HT" ~ "#08519c",
                     nt == "DA"   ~ "#a50f15",
                     nt == "NE"   ~ "#006d2c"
                   )),
    stim = factor(stim, levels = c("Pre-Stim", "Post-Stim"))
  )

# Panel B: UG early vs late (subject-level)
early_late_UG <- source_data %>%
  filter(panel == "B") %>%
  mutate(
    nt = factor(nt, levels = c("DA", "5-HT", "NE")),
    color = ifelse(stim == "Pre-Stim",
                   case_when(
                     nt == "5-HT" ~ "#9ecae1",
                     nt == "DA"   ~ "#fcbba1",
                     nt == "NE"   ~ "#c7e9c0"
                   ),
                   case_when(
                     nt == "5-HT" ~ "#08519c",
                     nt == "DA"   ~ "#a50f15",
                     nt == "NE"   ~ "#006d2c"
                   ))
  )

# Panel C: RL early vs late (subject-level)
early_late_RL <- source_data %>%
  filter(panel == "C") %>%
  mutate(
    nt = factor(nt, levels = c("DA", "5-HT", "NE")),
    color = ifelse(stim == "Pre-Stim",
                   case_when(
                     nt == "5-HT" ~ "#9ecae1",
                     nt == "DA"   ~ "#fcbba1",
                     nt == "NE"   ~ "#c7e9c0"
                   ),
                   case_when(
                     nt == "5-HT" ~ "#08519c",
                     nt == "DA"   ~ "#a50f15",
                     nt == "NE"   ~ "#006d2c"
                   ))
  )

# Panel A: Slope Comparison Plot
panel_a <-
  slope_data %>%
  ggplot(aes(x = task, y = slope)) +
  geom_col(position = position_dodge2(width = .5), aes(group = stim, fill = color)) +
  geom_errorbar(position = position_dodge2(width = .5),
                linewidth = 0.75,
                aes(ymin = slope_lower, ymax = slope_upper)) +
  facet_wrap(~nt, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Temporal slope \n(\u0394Estimate/trial)", x = "Task",
       fill = element_blank()) +
  theme_pubr(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(color = "black", size = 15),
    panel.spacing    = unit(0.2, "cm")
  ) +
  scale_fill_identity()

# Panel B: UG early vs late
panel_b <-
  early_late_UG %>%
  ggplot(aes(x = trial_period, y = Oz, color = color)) +
  geom_point(aes(group = interaction(idx, stim)), alpha = 0.3) +
  geom_line(aes(group = interaction(idx, stim)), alpha = 0.3) +
  stat_summary(aes(group = stim), size = 1, linewidth = 1) +
  stat_summary(geom = "line", aes(group = stim), size = 2) +
  facet_wrap(~nt) +
  scale_color_identity() +
  labs(y = "Estimate [z]", x = "Within task period",
       title = "Ultimatum game") +
  theme_pubr(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(color = "black", size = 14),
    panel.spacing    = unit(0.2, "cm"),
    axis.text.x      = element_text(angle = 15, vjust = 1, hjust = 0.5)
  )

# Panel C: RL early vs late
panel_c <-
  early_late_RL %>%
  ggplot(aes(x = trial_period, y = Oz, color = color)) +
  geom_point(aes(group = interaction(idx, stim)), alpha = 0.3) +
  geom_line(aes(group = interaction(idx, stim)), alpha = 0.3) +
  stat_summary(aes(group = stim), size = 1, linewidth = 1) +
  stat_summary(geom = "line", aes(group = stim), size = 2) +
  facet_wrap(~nt) +
  scale_color_identity() +
  labs(y = "Estimate [z]", x = "Within task period", title = "Reversal learning") +
  theme_pubr(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(color = "black", size = 14),
    panel.spacing    = unit(0.2, "cm"),
    axis.text.x      = element_text(angle = 15, vjust = 1, hjust = 0.5)
  )

create_prestim_legend <- function() {
  legend_data_pre <- data.frame(
    nt    = factor(c("DA", "5-HT", "NE"), levels = c("DA", "5-HT", "NE")),
    color = c("#fcbba1", "#9ecae1", "#c7e9c0"),
    x = 1,
    y = 1
  )

  ggplot(legend_data_pre, aes(x = x, y = y, color = nt)) +
    geom_point(size = 0) +
    geom_line(size = 2) +
    scale_color_manual(
      name   = "Pre-Stim",
      values = c("DA" = "#fcbba1", "5-HT" = "#9ecae1", "NE" = "#c7e9c0")
    ) +
    theme_void() +
    theme(
      legend.title      = element_text(face = "bold", size = 12),
      legend.text       = element_text(size = 10),
      legend.key.width  = unit(1, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.margin     = margin(0, 0, 10, 0)
    ) +
    guides(color = guide_legend(
      override.aes   = list(size = 2, alpha = 1),
      title.position = "top",
      title.hjust    = 0.5
    ))
}

create_poststim_legend <- function() {
  legend_data_post <- data.frame(
    nt    = factor(c("DA", "5-HT", "NE"), levels = c("DA", "5-HT", "NE")),
    color = c("#a50f15", "#08519c", "#006d2c"),
    x = 1,
    y = 1
  )

  ggplot(legend_data_post, aes(x = x, y = y, color = nt)) +
    geom_point(size = 0) +
    geom_line(size = 2) +
    scale_color_manual(
      name   = "Post-Stim",
      values = c("DA" = "#a50f15", "5-HT" = "#08519c", "NE" = "#006d2c")
    ) +
    theme_void() +
    theme(
      legend.title      = element_text(face = "bold", size = 12),
      legend.text       = element_text(size = 10),
      legend.key.width  = unit(1, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.margin     = margin(10, 0, 0, 0)
    ) +
    guides(color = guide_legend(
      override.aes   = list(size = 2, alpha = 1),
      title.position = "top",
      title.hjust    = 0.5
    ))
}

prestim_legend  <- get_legend(create_prestim_legend())
poststim_legend <- get_legend(create_poststim_legend())

stacked_legends <- plot_grid(prestim_legend, poststim_legend,
                             ncol = 1, align = "v")

full_supp_temporal_fig <- (panel_c + panel_b) +
  plot_annotation(tag_levels = "a")

final_plot <- plot_grid(full_supp_temporal_fig, stacked_legends,
                        ncol = 2, rel_widths = c(4, 1),
                        align = "h", axis = "tb")

ggsave(here::here("results/from_source_data/figureS2.png"),
       plot   = final_plot,
       device = "png",
       width  = 13,
       height = 6,
       dpi    = 300)
