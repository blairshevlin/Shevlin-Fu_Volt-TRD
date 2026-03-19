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
# 2025/06/24    Blair Shevlin                           wrote original code
# 2025/07/24    Blair Shevlin                           updated to use new NT data
# 2025/11/04    Blair Shevlin                           updated to use revised NT data
# 2026/03/19    Blair Shevlin                           refactored to read from source data CSV

# Reads source data from data/figures/figureExtended1_source_data.csv

rm(list = ls())

library(tidyverse)
library(fs)
library(patchwork)
library(ggpubr)
library(here)
library(RColorBrewer)
library(lmerTest)
library(caTools)
library(cowplot)
library(rstatix)
library(scales)
library(ggrepel)
library(effectsize)
library(emmeans)

# Read source data
source_data <- read.csv(here::here("data/figures/figureExtended1_source_data.csv"), stringsAsFactors = FALSE)

# Re-apply factor orderings
source_data$stim <- factor(source_data$stim, levels = c("Pre-Stim", "Post-Stim"))

# IDs
ids_final <- unique(source_data$idx)

# Specify colors
id_colors <- data.frame(id = ids_final,
                        color = brewer.pal(n = 10, name = "Paired"),
                        shape = c(0:9))

NT_colors <- data.frame(id = c("DA", "5-HT", "NE"),
                        color = c("#cb181d", "#2171b5", "darkgreen"))

# Panel A: RL subject-level averages
rl.avg <- source_data %>%
  filter(panel == "A") %>%
  mutate(color = ifelse(stim == "Pre-Stim", "#c7e9c0", "#006d2c"))

# Panel B: UG subject-level averages
ug.avg <- source_data %>%
  filter(panel == "B") %>%
  mutate(color = ifelse(stim == "Pre-Stim", "#c7e9c0", "#006d2c"))

# Calculate t-statistics
rl.NE.t <- t.test(rl.avg$Oz[rl.avg$stim == "Pre-Stim"],
                  rl.avg$Oz[rl.avg$stim == "Post-Stim"],
                  paired = TRUE)
ug.NE.t <- t.test(ug.avg$Oz[ug.avg$stim == "Pre-Stim"],
                  ug.avg$Oz[ug.avg$stim == "Post-Stim"],
                  paired = TRUE)
print(rl.NE.t)
print(ug.NE.t)

offset_data_RL <- rl.avg %>%
  mutate(
    point_x_pos = ifelse(stim == "Pre-Stim", 1 + 0.1, 2 - 0.1),
    box_x_pos   = ifelse(stim == "Pre-Stim", 1 - 0.1, 2 + 0.1),
    sess_num    = ifelse(stim == "Pre-Stim", 1, 2)
  )

offset_data_UG <- ug.avg %>%
  mutate(
    point_x_pos = ifelse(stim == "Pre-Stim", 1 + 0.1, 2 - 0.1),
    box_x_pos   = ifelse(stim == "Pre-Stim", 1 - 0.1, 2 + 0.1),
    sess_num    = ifelse(stim == "Pre-Stim", 1, 2)
  )

fig.rl.avg <-
  ggplot() +
  theme_pubr(base_size = 14) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "gray",
             linewidth = 1) +
  stat_summary(
    data = offset_data_RL,
    aes(x = box_x_pos, y = Oz, group = stim, color = color),
    linewidth = 1.5, size = 2, shape = 18, stroke = 2,
  ) +
  geom_line(
    data = offset_data_RL,
    aes(x = point_x_pos, y = Oz, group = idx),
    linewidth = 1.5,
    alpha = .25) +
  geom_point(
    data = offset_data_RL,
    aes(x = point_x_pos, y = Oz, group = idx, color = color),
    size = 4,
    stroke = 2,
    alpha = .5,
    show.legend = F, shape = 1
  ) +
  coord_cartesian(ylim = c(-7, 7)) +
  scale_color_identity() +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Pre-Stim", "Post-Stim"),
    limits = c(0.5, 2.5)
  ) +
  labs(x = element_blank(), y = "Estimate [z]", title = "Reversal learning") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        strip.text = element_text(color = "black"),
        panel.spacing = unit(0.2, "cm")
  )

fig.ug.avg <-
  ggplot() +
  theme_pubr(base_size = 14) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "gray",
             linewidth = 1) +
  stat_summary(
    data = offset_data_UG,
    aes(x = box_x_pos, y = Oz, group = stim, color = color),
    linewidth = 1.5, size = 2, shape = 18, stroke = 2,
  ) +
  geom_line(
    data = offset_data_UG,
    aes(x = point_x_pos, y = Oz, group = idx),
    linewidth = 1.5,
    alpha = .25) +
  geom_point(
    data = offset_data_UG,
    aes(x = point_x_pos, y = Oz, group = idx, color = color),
    size = 4,
    stroke = 2,
    alpha = .5, shape = 1,
    show.legend = F
  ) +
  coord_cartesian(ylim = c(-7, 7)) +
  scale_color_identity() +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Pre-Stim", "Post-Stim"),
    limits = c(0.5, 2.5)
  ) +
  labs(x = element_blank(), y = "Estimate [z]", title = "Ultimatum game") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        strip.text = element_text(color = "black"),
        panel.spacing = unit(0.2, "cm")
  )

# Panel C: RL trial-level trajectories
rl.trial.plot <- source_data %>%
  filter(panel == "C") %>%
  mutate(stim  = factor(stim, levels = c("Pre-Stim", "Post-Stim")),
         color = ifelse(stim == "Pre-Stim", "#c7e9c0", "#006d2c"))

fig.rl.trial <-
  rl.trial.plot %>%
  ggplot(aes(x = trial, y = Oz, fill = color, group = stim)) +
  theme_pubr(base_size = 14) +
  scale_color_identity() +
  scale_fill_identity() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  stat_summary(aes(color = color), geom = "line", size = 1, linewidth = 1.5) +
  stat_summary(
    geom = "ribbon", alpha = 0.35,
    show.legend = FALSE
  ) +
  labs(y = "NE Estimate [z]", x = "Trial", title = "Reversal learning") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        strip.text = element_text(color = "black"),
        panel.spacing = unit(0.2, "cm")) +
  coord_cartesian(ylim = c(-7, 7)) +
  theme(legend.key = element_rect(color = NA))

# Panel D: UG trial-level trajectories
ug.trial.plot <- source_data %>%
  filter(panel == "D") %>%
  mutate(stim  = factor(stim, levels = c("Pre-Stim", "Post-Stim")),
         color = ifelse(stim == "Pre-Stim", "#c7e9c0", "#006d2c"))

fig.ug.trial <-
  ug.trial.plot %>%
  ggplot(aes(x = trial, y = Oz, fill = color, group = stim)) +
  theme_pubr(base_size = 14) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  stat_summary(aes(color = color), geom = "line", size = 1, linewidth = 1.5) +
  stat_summary(
    geom = "ribbon", alpha = 0.35, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_identity() +
  labs(y = "NE Estimate [z]", x = "Trial", title = "Ultimatum game") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        strip.text = element_text(color = "black"),
        panel.spacing = unit(0.2, "cm")) +
  coord_cartesian(ylim = c(-7, 7)) +
  theme(legend.key = element_rect(color = NA))

main_plot <- ((fig.rl.avg / fig.ug.avg) | (fig.rl.trial / fig.ug.trial)) +
  plot_annotation(tag_levels = 'a')

create_legend <- function() {
  legend_data <- data.frame(
    stim  = factor(c("Pre-Stim", "Post-Stim"), levels = c("Pre-Stim", "Post-Stim")),
    color = c("#c7e9c0", "#006d2c"),
    x = 1,
    y = 1
  )

  ggplot(legend_data, aes(x = x, y = y, color = stim)) +
    geom_point(size = 0) +
    geom_line(size = 2) +
    scale_color_manual(
      name = element_blank(),
      values = c("Pre-Stim" = "#c7e9c0", "Post-Stim" = "#006d2c")
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

ne_legend <- get_legend(create_legend())

final_plot <- plot_grid(main_plot, ne_legend,
                        ncol = 2, rel_widths = c(4, 0.5),
                        align = "h", axis = "tb")

ggsave(here::here("results/from_source_data/Extended-Data_Figure1.png"),
       plot   = final_plot,
       device = "tiff",
       width  = 12,
       height = 8,
       dpi    = 300)
