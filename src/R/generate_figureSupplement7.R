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
# 2026/05/06    Blair Shevlin                           Remove old panel A

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
source_data <- read.csv(here::here("data/figures/figureSupplement7_source_data.csv"), stringsAsFactors = FALSE)

# Re-apply factor orderings
source_data$outcome_category <- factor(source_data$outcome_category,
                                       levels = c("Non-Responder","Responder","Remission"))

# Custom rounding function
midround <- function(x, base = 1) {
  return(base * round(x / base))
}

# ---- Panel A: RL DA vs 5-HT scatter -------------------------------------------
panel_A <- source_data %>%
  filter(panel == "A") %>%
  mutate(nt_pattern = factor(nt_pattern,
                             levels = c("Both Increase","Both Decrease","DA\u2191/5-HT\u2193","DA\u2193/5-HT\u2191")))

# Recompute change_magnitude for size aesthetic (HDRS_change stored in source data)
panel_A <- panel_A %>%
  mutate(change_magnitude = midround(HDRS_change, base = 5))

fig.profiles.rl <-
  ggplot(panel_A, aes(x = deltaDA, y = deltaSE)) +
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
                        breaks = c(-5, -10, -15),
                        guide = guide_legend(title.position = "top", title.hjust = 0.5,
                                             override.aes = list(fill = "black", alpha = 0.5))) +
  labs(x = "Change in DA (\u0394DA)", y = "Change in 5-HT (\u03945-HT)") +
  annotate("text", x =  6, y =  5, label = "Both \u2191",       size = 4, fontface = "bold") +
  annotate("text", x = -6, y = -5, label = "Both \u2193",       size = 4, fontface = "bold") +
  annotate("text", x =  6, y = -5, label = "DA\u2191/5-HT\u2193", size = 4, fontface = "bold") +
  annotate("text", x = -6, y =  5, label = "DA\u2193/5-HT\u2191", size = 4, fontface = "bold") +
  theme_pubr(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, color = "grey40"),
        legend.position = "right")

# ---- Panel B: RL synergy vs HDRS ----------------------------------------------
panel_B <- source_data %>%
  filter(panel == "B", !is.na(synergy_score)) %>%
  mutate(change_magnitude = midround(HDRS_change_m6, base = 5))

fig.synergy.rl <-
  ggplot(panel_B, aes(x = synergy_score, y = HDRS_m6)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "black", alpha = 0.2) +
  geom_point(aes(fill = nt_pattern, size = change_magnitude),
             shape = 21, color = "black", stroke = 1.5, alpha = 0.8) +
  geom_text(x = max(panel_B$synergy_score), y = max(panel_B$HDRS_m6),
            label = paste0("Pearson's r = ", round(cor(panel_B$synergy_score, panel_B$HDRS_m6), 3),
                           "\np = ", round(cor.test(panel_B$synergy_score, panel_B$HDRS_m6)$p.value, 3)),
            hjust = 1.1, vjust = 1.1, size = 5, color = "grey40") +
  geom_hline(yintercept = 7, linetype = "dashed", linewidth = 0.75) +
  geom_text(x = min(panel_B$synergy_score), y = 6, label = "Clinical Remission", size = 5, hjust = 0) +
  scale_fill_manual(values = c("Both Increase" = "#2166ac", "Both Decrease" = "#762a83",
                               "DA\u2191/5-HT\u2193" = "#d73027", "DA\u2193/5-HT\u2191" = "#f4a582"),
                    name = "Change Pattern") +
  scale_size_continuous(range = c(12, 3), name = "\u0394 HDRS-17",
                        breaks = c(-5, -10, -15),
                        guide = guide_legend(title.position = "top", title.hjust = 0.5,
                                             override.aes = list(fill = "black", alpha = 0.5))) +
  labs(x = "\u0394DA \u00d7 \u0394SE",
       y = "Month 6 HDRS-17") +
  theme_pubr(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, color = "grey40"),
        legend.position = "right")

# ---- Assemble ---------------------------------------------------------------
fig.S7 <- (fig.profiles.rl + fig.synergy.rl) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(plot.tag = element_text(size = 22, face = "bold"),
        plot.title    = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 12, color = "grey40")) 

ggsave(file.path(res_dir, "figureS7.png"),
       plot = fig.S7, device = "png",
       width = 16, height = 8, dpi = 300)
