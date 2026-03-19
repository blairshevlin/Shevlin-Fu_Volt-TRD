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
# 2025/07/10    Blair Shevlin                           wrote original code
# 2025/07/14    Blair Shevlin                           updated to use new NT data
# 2025/09/24    Blair Shevlin                           updated to generate subject-level avgs by task
# 2026/03/19    Blair Shevlin                           refactored to read from source data CSV

# Reads source data from data/figures/figure4_source_data.csv

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

# Paths
res_dir = here::here("results/from_source_data")

# Load source data
source_data <- read.csv(here::here("data/figures/figure4_source_data.csv"), stringsAsFactors = FALSE)

# Re-apply factor orderings
source_data$outcome_category <- factor(source_data$outcome_category,
                                       levels = c("Non-Responder","Responder","Remission"))

# ---- Panel A: HDRS trajectories --------------------------------------------
panel_A <- source_data %>%
  filter(panel == "A") %>%
  mutate(idx      = factor(idx),
         sess_fig = factor(session,
                           levels = c("Baseline","DBS","Week 1","Month 1","Month 2",
                                      "Month 3","Month 4","Month 5","Month 6")))

offset_data <- panel_A %>%
  filter(sess_fig %in% c("Baseline","Month 6")) %>%
  mutate(point_x_pos = ifelse(sess_fig == "Baseline", 1 + 0.1, 2 - 0.1),
         box_x_pos   = ifelse(sess_fig == "Baseline", 1 - 0.1, 2 + 0.1),
         sess_num    = ifelse(sess_fig == "Baseline", 1, 2))

fig.hdrs.change <-
  ggplot() +
  theme_pubr(base_size = 14) +
  geom_hline(yintercept = 8, linetype = "dashed", linewidth = 0.75) +
  annotate("text", x = 0.5, y = 9.5, label = "Clinical Remission", size = 5, hjust = 0) +
  geom_boxplot(data = offset_data,
               aes(x = box_x_pos, y = HDRS, group = sess_fig),
               linewidth = 1.1, outlier.alpha = 0, width = 0.15, show.legend = FALSE) +
  geom_line(data = offset_data,
            aes(x = point_x_pos, y = HDRS, group = idx),
            linewidth = 1.5, alpha = .25, color = "black") +
  geom_point(data = offset_data,
             aes(x = point_x_pos, y = HDRS, group = idx),
             size = 3, stroke = 1.5, alpha = .35, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, 32)) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Baseline","Month 6"), limits = c(0.5, 2.5)) +
  labs(x = element_blank(), y = "HDRS-17", title = "Month 6 Symptom Change")

# ---- Panel B: DA vs 5-HT scatter -------------------------------------------
panel_B <- source_data %>%
  filter(panel == "B") %>%
  mutate(nt_pattern = factor(nt_pattern,
                             levels = c("Both Increase","Both Decrease","DA\u2191/5-HT\u2193","DA\u2193/5-HT\u2191")))

# Recompute change_magnitude for size aesthetic (HDRS_change stored in source data)
panel_B <- panel_B %>%
  mutate(change_magnitude = HDRS_change)

fig.profiles <-
  ggplot(panel_B, aes(x = deltaDA, y = deltaSE)) +
  annotate("rect", xmin = 0,    xmax = Inf,  ymin = 0,    ymax = Inf,  alpha = 0.5, fill = "#2166ac") +
  annotate("rect", xmin = -Inf, xmax = 0,    ymin = -Inf, ymax = 0,    alpha = 0.5, fill = "#762a83") +
  annotate("rect", xmin = 0,    xmax = Inf,  ymin = -Inf, ymax = 0,    alpha = 0.5, fill = "#d73027") +
  annotate("rect", xmin = -Inf, xmax = 0,    ymin = 0,    ymax = Inf,  alpha = 0.5, fill = "#f4a582") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_point(aes(size = change_magnitude, fill = nt_pattern),
             shape = 21, color = "black", stroke = 1.5, alpha = 0.8) +
  scale_fill_manual(values = c("Both Increase" = "#2166ac", "Both Decrease" = "#762a83",
                               "DA\u2191/5-HT\u2193" = "#d73027", "DA\u2193/5-HT\u2191" = "#f4a582"),
                    name = "Change Pattern") +
  scale_size_continuous(range = c(12, 3), name = "\u0394 HDRS-17",
                        breaks = c(-10, -15, -20),
                        guide = guide_legend(title.position = "top", title.hjust = 0.5,
                                             override.aes = list(fill = "black", alpha = 0.5))) +
  labs(x = "Change in DA (\u0394DA)", y = "Change in 5-HT (\u03945-HT)",
       title = "Individual Patient Neurochemical Profiles") +
  annotate("text", x =  6, y =  5, label = "Both \u2191",       size = 4, fontface = "bold") +
  annotate("text", x = -6, y = -5, label = "Both \u2193",       size = 4, fontface = "bold") +
  annotate("text", x =  6, y = -5, label = "DA\u2191/5-HT\u2193", size = 4, fontface = "bold") +
  annotate("text", x = -6, y =  5, label = "DA\u2193/5-HT\u2191", size = 4, fontface = "bold") +
  theme_pubr(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, color = "grey40"),
        legend.position = "right")

# ---- Panel C: synergy vs HDRS ----------------------------------------------
panel_C <- source_data %>%
  filter(panel == "C") %>%
  mutate(change_magnitude = HDRS_change_m6)

fig.synergy <-
  ggplot(panel_C, aes(x = synergy_score, y = HDRS_m6)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "black", alpha = 0.2) +
  geom_point(aes(fill = nt_pattern, size = change_magnitude),
             shape = 21, color = "black", stroke = 1.5, alpha = 0.8) +
  geom_text(x = max(panel_C$synergy_score), y = max(panel_C$HDRS_m6),
            label = paste0("Pearson's r = ", round(cor(panel_C$synergy_score, panel_C$HDRS_m6), 3),
                           "\np = ", round(cor.test(panel_C$synergy_score, panel_C$HDRS_m6)$p.value, 3)),
            hjust = 1.1, vjust = 1.1, size = 5, color = "grey40") +
  geom_hline(yintercept = 8, linetype = "dashed", linewidth = 0.75) +
  geom_text(x = min(panel_C$synergy_score), y = 9, label = "Clinical Remission", size = 5, hjust = 0) +
  scale_fill_manual(values = c("Both Increase" = "#2166ac", "Both Decrease" = "#762a83",
                               "DA\u2191/5-HT\u2193" = "#d73027", "DA\u2193/5-HT\u2191" = "#f4a582"),
                    name = "Change Pattern") +
  scale_size_continuous(range = c(12, 3), name = "\u0394 HDRS-17",
                        breaks = c(-10, -15, -20),
                        guide = guide_legend(title.position = "top", title.hjust = 0.5,
                                             override.aes = list(fill = "black", alpha = 0.5))) +
  labs(x = "\u0394DA \u00d7 \u0394SE",
       y = "Month 6 HDRS-17",
       title = "Neurochemical Profiles Predict Treatment Response") +
  theme_pubr(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, color = "grey40"),
        legend.position = "right")

# ---- Assemble ---------------------------------------------------------------
fig.4 <- (fig.hdrs.change + fig.profiles + fig.synergy) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.title    = element_text(size = 16, face = "bold"),
                                plot.subtitle = element_text(size = 12, color = "grey40"))) +
  plot_layout(guides = "collect")

ggsave(file.path(res_dir, "figure4.png"),
       plot = fig.4, device = "png",
       width = 16, height = 5, dpi = 300)
