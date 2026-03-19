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
# 2025/09/16    Blair Shevlin                           updated to use revised NT data
# 2025/11/10    Blair Shevlin                           updated figure aesthetics
# 2026/03/19    Blair Shevlin                           refactored to read from source data CSV

# Reads source data from data/figures/figure2_source_data.csv

rm(list = ls())

library(tidyverse)
library(fs)
library(here)
library(patchwork)
library(ggpubr)
library(RColorBrewer)
library(GGally)
library(purrr)
library(broom)
library(caTools)
library(rstatix)
library(scales)
library(ggrepel)
library(cowplot)

# Paths
res_dir = here::here("results/from_source_data")

# Specify colors
NT_colors = data.frame(id = c("DA","5-HT"),
                       color = c("#cb181d","#2171b5"))

# Load source data
source_data <- read.csv(here::here("data/figures/figure2_source_data.csv"), stringsAsFactors = FALSE)

# Re-apply factor orderings
source_data$stim <- factor(source_data$stim, levels = c("Pre-Stim","Post-Stim"))
source_data$nt   <- factor(source_data$nt,   levels = c("DA","5-HT"))

# ---- Panel A: trial-level NT traces ----------------------------------------
panel_A <- source_data %>% filter(panel == "A")

panel_A <- panel_A %>%
  mutate(color = ifelse(stim == "Pre-Stim",
                        ifelse(nt == "5-HT","#9ecae1","#fcbba1"),
                        ifelse(nt == "5-HT","#08519c","#a50f15")))

fig.rl.trial <-
  panel_A %>%
  filter(task == "RL") %>%
  ggplot(aes(x = trial, y = Oz, fill = color, group = interaction(nt, stim))) +
  theme_pubr(base_size = 14) +
  facet_wrap(~nt) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  stat_summary(aes(color = color), geom = "line", size = 1, linewidth = 1.5) +
  stat_summary(geom = "ribbon", alpha = 0.35, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_identity() +
  labs(y = "Estimate [z]", x = "Trial", title = "Reversal learning") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        strip.text = element_text(color = "black"),
        panel.spacing = unit(0.2, "cm")) +
  coord_cartesian(ylim = c(-7, 7)) +
  guides(color = "none", fill = "none") +
  theme(legend.key = element_rect(color = NA))

fig.ug.trial <-
  panel_A %>%
  filter(task == "UG") %>%
  ggplot(aes(x = trial, y = Oz, fill = color, group = interaction(nt, stim))) +
  theme_pubr(base_size = 14) +
  facet_wrap(~nt) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  stat_summary(aes(color = color), geom = "line", size = 1, linewidth = 1.5) +
  stat_summary(geom = "ribbon", alpha = 0.35, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_identity() +
  labs(y = "Estimate [z]", x = "Trial", title = "Ultimatum game") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        strip.text = element_text(color = "black"),
        panel.spacing = unit(0.2, "cm")) +
  coord_cartesian(ylim = c(-7, 7)) +
  guides(color = "none", fill = "none") +
  theme(legend.key = element_rect(color = NA))

# ---- Panel B: subject-level averages by stim --------------------------------
panel_B <- source_data %>% filter(panel == "B")

panel_B <- panel_B %>%
  mutate(color = ifelse(stim == "Pre-Stim",
                        ifelse(nt == "5-HT","#9ecae1","#fcbba1"),
                        ifelse(nt == "5-HT","#08519c","#a50f15")),
         point_x_pos = ifelse(stim == "Pre-Stim", 1 + 0.1, 2 - 0.1),
         box_x_pos   = ifelse(stim == "Pre-Stim", 1 - 0.1, 2 + 0.1),
         sess_num    = ifelse(stim == "Pre-Stim", 1, 2))

offset_data_RL <- panel_B %>% filter(task == "RL")
offset_data_UG <- panel_B %>% filter(task == "UG")

create_prestim_legend <- function() {
  legend_data_pre <- data.frame(
    nt    = factor(c("DA","5-HT"), levels = c("DA","5-HT")),
    color = c("#fcbba1","#9ecae1"),
    x = 1, y = 1
  )
  ggplot(legend_data_pre, aes(x = x, y = y, color = nt)) +
    geom_point(size = 0) +
    geom_line(size = 2) +
    scale_color_manual(name = "Pre-Stim", values = c("DA" = "#fcbba1","5-HT" = "#9ecae1")) +
    theme_void() +
    theme(legend.title = element_text(face = "bold", size = 12),
          legend.text  = element_text(size = 10),
          legend.key.width  = unit(1, "cm"),
          legend.key.height = unit(0.5, "cm"),
          legend.margin = margin(0, 0, 10, 0)) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1),
                                title.position = "top", title.hjust = 0.5))
}

create_poststim_legend <- function() {
  legend_data_post <- data.frame(
    nt    = factor(c("DA","5-HT"), levels = c("DA","5-HT")),
    color = c("#a50f15","#08519c"),
    x = 1, y = 1
  )
  ggplot(legend_data_post, aes(x = x, y = y, color = nt)) +
    geom_point(size = 0) +
    geom_line(size = 2) +
    scale_color_manual(name = "Post-Stim", values = c("DA" = "#a50f15","5-HT" = "#08519c")) +
    theme_void() +
    theme(legend.title = element_text(face = "bold", size = 12),
          legend.text  = element_text(size = 10),
          legend.key.width  = unit(1, "cm"),
          legend.key.height = unit(0.5, "cm"),
          legend.margin = margin(10, 0, 0, 0)) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1),
                                title.position = "top", title.hjust = 0.5))
}

prestim_legend  <- get_legend(create_prestim_legend())
poststim_legend <- get_legend(create_poststim_legend())
stacked_legends <- plot_grid(prestim_legend, poststim_legend, ncol = 1, align = "v")

fig.rl.avg <-
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  geom_line(data = offset_data_RL,
            aes(x = point_x_pos, y = Oz, group = idx),
            position = position_dodge2(width = .2), size = 1.5, alpha = .15) +
  geom_point(data = offset_data_RL,
             aes(x = point_x_pos, y = Oz, group = idx, color = color),
             position = position_dodge2(width = .2), shape = 1, stroke = 2, alpha = 0.5, size = 3) +
  stat_summary(geom = "line", data = offset_data_RL,
               aes(x = box_x_pos, y = Oz, group = nt), linewidth = 1.1) +
  stat_summary(data = offset_data_RL,
               aes(x = box_x_pos, y = Oz, group = stim, color = color),
               size = 1.1, linewidth = 1.1) +
  facet_wrap(~ nt) +
  scale_color_identity() +
  labs(y = "Estimate [z]", x = element_blank(), title = "Reversal learning") +
  theme_pubr(base_size = 14) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Pre-Stim","Post-Stim"), limits = c(0.5, 2.5)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black", size = 14),
        panel.spacing = unit(0.2, "cm"),
        axis.text.x = element_text(angle = 15, vjust = 0.95, hjust = 0.5))

fig.ug.avg <-
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  geom_line(data = offset_data_UG,
            aes(x = point_x_pos, y = Oz, group = idx),
            position = position_dodge2(width = .2), size = 1.5, alpha = .15) +
  geom_point(data = offset_data_UG,
             aes(x = point_x_pos, y = Oz, group = idx, color = color),
             position = position_dodge2(width = .2), shape = 1, stroke = 2, alpha = 0.5, size = 3) +
  stat_summary(geom = "line", data = offset_data_UG,
               aes(x = box_x_pos, y = Oz, group = nt), linewidth = 1.1) +
  stat_summary(data = offset_data_UG,
               aes(x = box_x_pos, y = Oz, group = stim, color = color),
               size = 1.1, linewidth = 1.1) +
  facet_wrap(~ nt) +
  scale_color_identity() +
  labs(y = "Estimate [z]", x = element_blank(), title = "Ultimatum game") +
  theme_pubr(base_size = 14) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Pre-Stim","Post-Stim"), limits = c(0.5, 2.5)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black", size = 14),
        panel.spacing = unit(0.2, "cm"),
        axis.text.x = element_text(angle = 15, vjust = 0.95, hjust = 0.5))

# ---- Panel C: subject-level by condition/offer bin --------------------------
panel_C <- source_data %>% filter(panel == "C")

panel_C <- panel_C %>%
  mutate(color = ifelse(stim == "Pre-Stim",
                        ifelse(nt == "5-HT","#9ecae1","#fcbba1"),
                        ifelse(nt == "5-HT","#08519c","#a50f15")))

offset_data_RL_task <- panel_C %>%
  filter(task == "RL") %>%
  mutate(condition = factor(condition,
                            levels = c("Positive","Mixed","Negative"),
                            labels = c("Reward","Mixed","Punishment")),
         block_num     = as.numeric(condition),
         point_x_pos   = ifelse(stim == "Pre-Stim", block_num - 0.15, block_num + 0.15),
         group_x_pos   = ifelse(stim == "Pre-Stim", block_num - 0.35, block_num + 0.35))

offset_data_UG_task <- panel_C %>%
  filter(task == "UG") %>%
  mutate(condition = factor(condition, levels = c("$1-3","$4-6","$7-9")),
         offer_num     = as.numeric(condition),
         point_x_pos   = ifelse(stim == "Pre-Stim", offer_num - 0.15, offer_num + 0.15),
         group_x_pos   = ifelse(stim == "Pre-Stim", offer_num - 0.35, offer_num + 0.35))

fig.rl.task <-
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  theme_pubr(base_size = 14) +
  facet_wrap(~nt) +
  geom_point(data = offset_data_RL_task,
             aes(x = point_x_pos, y = Oz, color = color),
             alpha = .75, size = 1.5, stroke = 1.5, shape = 1) +
  stat_summary(data = offset_data_RL_task,
               aes(x = group_x_pos, y = Oz, color = color, group = stim),
               size = 1.25, linewidth = 1.25, position = position_identity()) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_x_continuous(breaks = c(1, 2, 3),
                     labels = c("Reward","Mixed","Punishment"),
                     limits = c(0.5, 3.5)) +
  labs(y = "Estimate [z]", x = "Block type", title = "Reversal learning") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        strip.text = element_text(color = "black"),
        panel.spacing = unit(0.2, "cm"),
        axis.text.x = element_text(angle = 15, vjust = 0.9, hjust = 0.6)) +
  coord_cartesian(ylim = c(-7.5, 7.5)) +
  guides(color = "none", fill = "none") +
  theme(legend.key = element_rect(color = NA))

fig.ug.task <-
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  theme_pubr(base_size = 14) +
  facet_wrap(~nt) +
  geom_point(data = offset_data_UG_task,
             aes(x = point_x_pos, y = Oz, group = idx, color = color),
             shape = 1, stroke = 1.75, alpha = 0.75, size = 1.5) +
  stat_summary(data = offset_data_UG_task,
               aes(x = group_x_pos, y = Oz, group = stim, color = color),
               size = 1.25, linewidth = 1.25) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_x_continuous(breaks = c(1, 2, 3),
                     labels = c("$1-3","$4-6","$7-9"),
                     limits = c(0.5, 3.5)) +
  labs(y = "Estimate [z]", x = "Offer", title = "Ultimatum game") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        strip.text = element_text(color = "black"),
        panel.spacing = unit(0.2, "cm")) +
  coord_cartesian(ylim = c(-7.5, 7.5)) +
  guides(color = "none", fill = "none") +
  theme(legend.key = element_rect(color = NA))

# ---- Assemble ---------------------------------------------------------------
main_plot <-
  (fig.rl.trial | fig.ug.trial) /
  (fig.rl.avg   | fig.ug.avg) /
  (fig.rl.task  | fig.ug.task) +
  plot_annotation(tag_levels = "a")

final_plot <- plot_grid(main_plot, stacked_legends,
                        ncol = 2, rel_widths = c(4, 1),
                        align = "h", axis = "tb")

ggsave(file.path(res_dir, "figure2.png"),
       plot = final_plot, device = "png",
       width = 12, height = 10, dpi = 300)

panel_a <- plot_grid((fig.rl.trial | fig.ug.trial) +
                       plot_annotation(tag_levels = "a"), stacked_legends,
                     ncol = 2, rel_widths = c(4, 1),
                     align = "h", axis = "tb")

ggsave(file.path(res_dir, "figure2_panelA.png"),
       plot = panel_a, device = "png",
       width = 12, height = 4, dpi = 300)
