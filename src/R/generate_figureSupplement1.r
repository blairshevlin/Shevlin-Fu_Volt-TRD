
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
# 2025/09/23    Blair Shevlin                           wrote original code
# 2025/11/25    Blair Shevlin                           adapted for raw model outputs
# 2026/03/19    Blair Shevlin                           refactored to read from source data CSV

# Reads source data from data/figures/figureSupplement1_source_data.csv

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

# Function for calculating standard error
se <- function(x) {
  sd(x) / sqrt(length(x))
}

# Read source data
source_data <- read.csv(here::here("data/figures/figureSupplement1_source_data.csv"), stringsAsFactors = FALSE)

# Re-apply factor orderings
source_data$stim <- factor(source_data$stim, levels = c("Pre-Stim", "Post-Stim"))
source_data$nt   <- factor(source_data$nt,   levels = c("DA", "5-HT"))

# Panels A & B: subject-level means (RL and UG)
offset_data_RL_raw <- source_data %>%
  filter(panel == "A") %>%
  mutate(
    point_x_pos = ifelse(stim == "Pre-Stim", 1 + 0.1, 2 - 0.1),
    box_x_pos   = ifelse(stim == "Pre-Stim", 1 - 0.1, 2 + 0.1)
  )

offset_data_UG_raw <- source_data %>%
  filter(panel == "B") %>%
  mutate(
    point_x_pos = ifelse(stim == "Pre-Stim", 1 + 0.1, 2 - 0.1),
    box_x_pos   = ifelse(stim == "Pre-Stim", 1 - 0.1, 2 + 0.1)
  )

# Panels C & D: correlation data (RL and UG, Pre-Stim only)
rl.cl.df <- source_data %>% filter(panel == "C")
ug.cl.df <- source_data %>% filter(panel == "D")

# Print summary statistics for the manuscript
cat("Summary Statistics for Raw Data:\n")
cat("Reversal Learning:\n")
print(offset_data_RL_raw %>%
        group_by(nt, stim) %>%
        summarise(
          mean_estimate  = mean(Oz, na.rm = TRUE),
          sd_estimate    = sd(Oz, na.rm = TRUE),
          .groups = 'drop'
        ))

cat("\nUltimatum Game:\n")
print(offset_data_UG_raw %>%
        group_by(nt, stim) %>%
        summarise(
          mean_estimate  = mean(Oz, na.rm = TRUE),
          sd_estimate    = sd(Oz, na.rm = TRUE),
          .groups = 'drop'
        ))

# Change statistics
rl.wilcox <- offset_data_RL_raw %>%
  group_by(nt) %>%
  summarise(
    wilcox_test = list(wilcox.test(
      Oz[stim == "Pre-Stim"],
      Oz[stim == "Post-Stim"],
      paired = TRUE
    )),
    .groups = 'drop'
  ) %>%
  mutate(tidy_results = map(wilcox_test, tidy)) %>%
  unnest(tidy_results) %>%
  select(-wilcox_test)

ug.wilcox <- offset_data_UG_raw %>%
  group_by(nt) %>%
  summarise(
    wilcox_test = list(wilcox.test(
      Oz[stim == "Pre-Stim"],
      Oz[stim == "Post-Stim"],
      paired = TRUE
    )),
    .groups = 'drop'
  ) %>%
  mutate(tidy_results = map(wilcox_test, tidy)) %>%
  unnest(tidy_results) %>%
  select(-wilcox_test)

print(rl.wilcox)
print(ug.wilcox)

# Run correlations and print results
cat("\nCorrelation Results:\n")
for (nt_type in c("DA", "5-HT")) {
  for (task_label in c("Reversal Learning", "Ultimatum Game")) {
    if (task_label == "Reversal Learning") {
      data_subset <- rl.cl.df %>% filter(nt == nt_type)
    } else {
      data_subset <- ug.cl.df %>% filter(nt == nt_type)
    }

    cor_test <- cor.test(data_subset$Oz, data_subset$HDRS, method = "spearman")

    cat(sprintf("%s - %s:\n", task_label, nt_type))
    cat(sprintf("  Pearson's r: %.3f\n", cor_test$estimate))
    cat(sprintf("  p-value: %.4f\n", cor_test$p.value))
    cat(sprintf("  n: %d\n\n", nrow(data_subset)))
  }
}

# Panel A: RL subject means plot
fig.rl.raw <-
  ggplot() +
  facet_wrap(~nt) +
  theme_pubr(base_size = 14) +
  geom_errorbar(
    data = offset_data_RL_raw,
    aes(x = point_x_pos, y = Oz,
        ymin = Oz - se,  ymax = Oz + se,
        color = color),
    width = 0.05, alpha = 0.6, linewidth = 0.8
  ) +
  geom_line(
    data = offset_data_RL_raw,
    aes(x = point_x_pos, y = Oz, group = idx),
    linewidth = 1, alpha = 0.4, color = "gray50"
  ) +
  geom_point(
    data = offset_data_RL_raw,
    aes(x = point_x_pos, y = Oz, color = color),
    size = 2.5, alpha = 0.7
  ) +
  stat_summary(
    data = offset_data_RL_raw,
    aes(x = box_x_pos, y = Oz, color = color),
    fun = mean,
    fun.min = function(x) mean(x) - se(x),
    fun.max = function(x) mean(x) + se(x),
    geom = "errorbar",
    linewidth = 2, width = 0.1, alpha = 1
  ) +
  stat_summary(
    data = offset_data_RL_raw,
    aes(x = box_x_pos, y = Oz, color = color),
    fun = mean, geom = "point",
    size = 6, shape = 18, alpha = 1
  ) +
  scale_color_identity() +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Pre-Stim", "Post-Stim"),
    limits = c(0.5, 2.5)
  ) +
  labs(x = element_blank(),
       y = "Estimate [Z]",
       title = "Reversal learning") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border     = element_rect(colour = NA, fill = NA),
    strip.text       = element_text(color = "black"),
    panel.spacing    = unit(0.2, "cm"),
    axis.text.x      = element_text(angle = 15, vjust = 1, hjust = 0.5)
  )

# Panel B: UG subject means plot
fig.ug.raw <-
  ggplot() +
  facet_wrap(~nt) +
  theme_pubr(base_size = 14) +
  geom_errorbar(
    data = offset_data_UG_raw,
    aes(x = point_x_pos, y = Oz,
        ymin = Oz - se, ymax = Oz + se,
        color = color),
    width = 0.05, alpha = 0.6, linewidth = 0.8
  ) +
  geom_line(
    data = offset_data_UG_raw,
    aes(x = point_x_pos, y = Oz, group = idx),
    linewidth = 1, alpha = 0.4, color = "gray50"
  ) +
  geom_point(
    data = offset_data_UG_raw,
    aes(x = point_x_pos, y = Oz, color = color),
    size = 2.5, alpha = 0.7
  ) +
  stat_summary(
    data = offset_data_UG_raw,
    aes(x = box_x_pos, y = Oz, color = color),
    fun = mean,
    fun.min = function(x) mean(x) - se(x),
    fun.max = function(x) mean(x) + se(x),
    geom = "errorbar",
    linewidth = 2, width = 0.1, alpha = 1
  ) +
  stat_summary(
    data = offset_data_UG_raw,
    aes(x = box_x_pos, y = Oz, color = color),
    fun = mean, geom = "point",
    size = 6, shape = 18, alpha = 1
  ) +
  scale_color_identity() +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Pre-Stim", "Post-Stim"),
    limits = c(0.5, 2.5)
  ) +
  labs(x = element_blank(),
       y = "Estimate [Z]",
       title = "Ultimatum game") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border     = element_rect(colour = NA, fill = NA),
    strip.text       = element_text(color = "black"),
    panel.spacing    = unit(0.2, "cm"),
    axis.text.x      = element_text(angle = 15, vjust = 1, hjust = 0.5)
  )

# Panel C: RL correlation with HDRS
fig.rl.cor <-
  ggplot(rl.cl.df,
         aes(x = HDRS, y = Oz, color = nt)) +
  geom_point(size = 6) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  facet_wrap(~nt) +
  theme_pubr(base_size = 14) +
  scale_color_manual(values = c("DA" = "#fcbba1", "5-HT" = "#9ecae1")) +
  labs(x = "HDRS", y = "Pre-Stim\nEstimate [Z]", title = "Reversal learning") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border     = element_rect(colour = NA, fill = NA),
    strip.text       = element_text(color = "black"),
    panel.spacing    = unit(0.2, "cm"),
    legend.position  = "none"
  ) +
  stat_cor(color = "black", method = "spearman", label.x.npc = "left", label.y = 200, size = 4)

# Panel D: UG correlation with HDRS
fig.ug.cor <-
  ggplot(ug.cl.df,
         aes(x = HDRS, y = Oz, color = nt)) +
  geom_point(size = 6) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  facet_wrap(~nt) +
  theme_pubr(base_size = 14) +
  scale_color_manual(values = c("DA" = "#fcbba1", "5-HT" = "#9ecae1")) +
  labs(x = "HDRS", y = "Pre-Stim\nEstimate [Z]", title = "Ultimatum game") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border     = element_rect(colour = NA, fill = NA),
    strip.text       = element_text(color = "black"),
    panel.spacing    = unit(0.2, "cm"),
    legend.position  = "none"
  ) +
  stat_cor(color = "black", method = "spearman", label.x.npc = "left", label.y = 70, size = 4)

create_prestim_legend <- function() {
  legend_data_pre <- data.frame(
    nt    = factor(c("DA", "5-HT"), levels = c("DA", "5-HT")),
    color = c("#fcbba1", "#9ecae1"),
    x = 1,
    y = 1
  )

  ggplot(legend_data_pre, aes(x = x, y = y, color = nt)) +
    geom_point(size = 0) +
    geom_line(size = 2) +
    scale_color_manual(
      name   = "Pre-Stim",
      values = c("DA" = "#fcbba1", "5-HT" = "#9ecae1")
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
    nt    = factor(c("DA", "5-HT"), levels = c("DA", "5-HT")),
    color = c("#a50f15", "#08519c"),
    x = 1,
    y = 1
  )

  ggplot(legend_data_post, aes(x = x, y = y, color = nt)) +
    geom_point(size = 0) +
    geom_line(size = 2) +
    scale_color_manual(
      name   = "Post-Stim",
      values = c("DA" = "#a50f15", "5-HT" = "#08519c")
    ) +
    theme_void() +
    theme(
      legend.title      = element_text(face = "bold", size = 12),
      legend.text       = element_text(size = 10),
      legend.key.width  = unit(1, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.margin     = margin(10, 0, 10, 0)
    ) +
    guides(color = guide_legend(
      override.aes   = list(size = 2, alpha = 1),
      title.position = "top",
      title.hjust    = 0.5
    ))
}

prestim_legend  <- get_legend(create_prestim_legend())
poststim_legend <- get_legend(create_poststim_legend())

stacked_legends <-
  plot_grid(prestim_legend, poststim_legend,
            ncol = 1, align = "v")

main_plot_raw <- (fig.rl.raw | fig.ug.raw) / (fig.rl.cor | fig.ug.cor) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.title    = element_text(size = 16, face = "bold"),
                                plot.subtitle = element_text(size = 12)))

final_plot_raw <- plot_grid(main_plot_raw, stacked_legends,
                            ncol = 2, rel_widths = c(4, 1),
                            align = "h", axis = "tb")

ggsave(here::here("results/from_source_data/figureS1.png"),
       plot   = final_plot_raw,
       device = "png",
       width  = 14,
       height = 10,
       dpi    = 300)
