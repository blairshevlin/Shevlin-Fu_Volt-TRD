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
# 2025/10/08    Blair Shevlin                           wrote original code
# 2026/03/19    Blair Shevlin                           refactored to read from source data CSV

rm(list = ls())

library(tidyverse)
library(fs)
library(patchwork)
library(ggpubr)
library(here)
library(RColorBrewer)
library(caTools)
library(cowplot)
library(rstatix)
library(scales)
library(ggrepel)
library(grid)
library(gridExtra)

# Paths
res_dir = here::here("results/from_source_data")

# Color scheme
NT_colors = data.frame(id = c("DA","SE"),
                       color = c("#cb181d","#2171b5"))

# Load source data
source_data <- read.csv(here::here("data/figures/figureExtended2_source_data.csv"), stringsAsFactors = FALSE)

fig.data <- source_data %>% filter(panel == "C")

window_sum <- fig.data %>%
  filter(in_auc_window == TRUE) %>%
  pull(z_score) %>%
  sum()

nM_range <- range(fig.data$nM_smoothed)
z_range  <- range(fig.data$z_score)
scale_factor      <- diff(nM_range) / diff(z_range)
intercept_adjust  <- nM_range[1] - z_range[1] * scale_factor

# Panel 1: Raw nM estimates
p1 <- ggplot(fig.data, aes(x = time_sample, y = nM_smoothed)) +
  theme_pubr(base_size = 14) +
  geom_vline(xintercept = 0, linewidth = 1.1, linetype = "dashed", color = "black") +
  geom_line(linewidth = 1.5, color = NT_colors$color[NT_colors$id == "SE"]) +
  labs(x = "time relative to outcome reveal [ms]",
       y = "5-HT [nM]",
       title = "raw model estimate") +
  coord_cartesian(xlim = c(-250, 750)) +
  annotate("text", x = -40, y = max(fig.data$nM_smoothed), label = "Outcome\nReveal",
           size = 3.5, hjust = 1, vjust = 1.2)

p2 <- ggplot(fig.data, aes(x = time_sample)) +
  theme_pubr(base_size = 14) +
  geom_hline(yintercept = 0 * scale_factor + intercept_adjust, linewidth = 1.1, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linewidth = 1.1, linetype = "dashed", color = "black") +
  geom_line(aes(y = z_score * scale_factor + intercept_adjust), linewidth = 1.5, color = NT_colors$color[NT_colors$id == "SE"]) +
  labs(x = "time relative to outcome reveal [ms]",
       y = "5-HT [nM]",
       title = "z-scored estimate") +
  coord_cartesian(xlim = c(-250, 750)) +
  scale_y_continuous(
    name = "5-HT [nM]",
    sec.axis = sec_axis(~ (. - intercept_adjust) / scale_factor, name = "5-HT [z]")
  ) +
  theme(axis.title.y.right = element_text(color = NT_colors$color[NT_colors$id == "SE"], angle = 90),
        axis.text.y.right  = element_text(color = NT_colors$color[NT_colors$id == "SE"]))

p3 <- ggplot(fig.data, aes(x = time_sample)) +
  theme_pubr(base_size = 14) +
  geom_hline(yintercept = 0 * scale_factor + intercept_adjust, linewidth = 1.1, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linewidth = 1.1, linetype = "dashed", color = "black") +
  geom_line(aes(y = z_score * scale_factor + intercept_adjust), linewidth = 1.5, color = NT_colors$color[NT_colors$id == "SE"]) +
  geom_ribbon(data = subset(fig.data, in_auc_window == TRUE),
              aes(ymax = z_score * scale_factor + intercept_adjust,
                  ymin = 0 * scale_factor + intercept_adjust),
              fill = NT_colors$color[NT_colors$id == "SE"],
              alpha = 0.25) +
  labs(x = "time relative to outcome reveal [ms]",
       y = "5-HT [nM]",
       title = "sum over 500 ms window") +
  coord_cartesian(xlim = c(-250, 750)) +
  scale_y_continuous(
    name = "5-HT [nM]",
    sec.axis = sec_axis(~ (. - intercept_adjust) / scale_factor, name = "5-HT [z]")
  ) +
  theme(axis.title.y.right = element_text(color = NT_colors$color[NT_colors$id == "SE"], angle = 90),
        axis.text.y.right  = element_text(color = NT_colors$color[NT_colors$id == "SE"]))

# Create simple arrow plots
arrow1 <- ggplot() +
  theme_void() +
  annotate("segment", x = 0.2, xend = 0.8, y = 0.5, yend = 0.5,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           size = 1.5) +
  xlim(0, 1) + ylim(0, 1)

arrow2 <- ggplot() +
  theme_void() +
  annotate("segment", x = 0.2, xend = 0.8, y = 0.5, yend = 0.5,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           size = 1.5) +
  xlim(0, 1) + ylim(0, 1)

# Combine with patchwork
final_figure <- p1 + arrow1 + p2 + arrow2 + p3 +
  plot_layout(widths = c(4, 0.5, 4, 0.5, 4), axis_titles = "collect")

# Save the figure
ggsave(file.path(res_dir, "Extended-Data_Figure2.png"),
       plot = final_figure,
       device = "png",
       width = 12,
       height = 4,
       dpi = 300)

ggsave(file.path(res_dir, "Extended-Data_Figure2.pdf"),
       plot = final_figure,
       device = "pdf",
       width = 12,
       height = 4,
       dpi = 300)
