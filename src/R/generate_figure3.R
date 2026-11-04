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
# 2025/10/08    Blair Shevlin                           include updated NT data
# 2026/03/19    Blair Shevlin                           refactored to read from source data CSV

# Reads source data from data/figures/figure3_source_data.csv

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
res_dir = here::here("results/from_source_data")

# Specify colors
NT_colors = data.frame(id = c("DA","5-HT"),
                       color = c("#cb181d","#2171b5"))

# Load source data
source_data <- read.csv(here::here("data/figures/figure3_source_data.csv"), stringsAsFactors = FALSE)

# Re-apply factor orderings for timeline panels
sess_levels <- c("Baseline","Pre-Stim","DBS","Post-Stim","Week 1",
                 "Month 1","Month 2","Month 3","Month 4","Month 5","Month 6")

# ---- Top panels: behavioral timelines (A-E) ---------------------------------

# RL choice (panel = "RL_choice")
rl.beh.means <- source_data %>%
  filter(panel == "RL_choice") %>%
  mutate(idx  = factor(idx),
         sess = factor(sess, levels = sess_levels))

rl.c.sess <-
  rl.beh.means %>%
  filter(sess %in% c("Baseline","DBS","Month 1","Month 3","Month 6")) %>%
  ggplot(aes(x = sess, y = value)) +
  theme_pubr(base_size = 7) +
  geom_vline(xintercept = "DBS", linewidth = 0.5, color = "grey70") +
  geom_boxplot(data = rl.beh.means[rl.beh.means$sess %in% c("Baseline","DBS","Month 1","Month 3","Month 6"),],
               linewidth = 0.4, outlier.alpha = 0, show.legend = FALSE) +
  geom_point(data = rl.beh.means[rl.beh.means$sess %in% c("Baseline","Month 1","Month 3","Month 6"),],
             size = 0.8, color = "purple", alpha = .5,
             position = position_dodge2(width = .3), stroke = 0.6, aes(group = idx)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = element_blank(), y = "p(optimal)", title = "Choice [RL]",
       color = element_blank(), shape = element_blank()) +
  scale_shape_manual(values = c(16, 4)) +
  stat_compare_means(label = "p.signif",
                     comparisons = list(c("Baseline","Month 1"), c("Baseline","Month 3"), c("Baseline","Month 6")),
                     method = "t.test", paired = TRUE,  # rows are sorted by idx in every session
                     bracket.size = 0.3, tip.length = 0.02, size = 2.5,
                     step.increase = 0.12)

# RL RT (panel = "RL_RT")
rl.rt.means <- source_data %>%
  filter(panel == "RL_RT") %>%
  mutate(idx  = factor(idx),
         sess = factor(sess, levels = sess_levels))

rl.logrt.sess <-
  rl.rt.means %>%
  filter(sess %in% c("Baseline","DBS","Month 1","Month 3","Month 6")) %>%
  ggplot(aes(x = sess, y = value)) +
  theme_pubr(base_size = 7) +
  geom_vline(xintercept = "DBS", linewidth = 0.5, color = "grey70") +
  geom_boxplot(data = rl.rt.means[rl.rt.means$sess %in% c("Baseline","DBS","Month 1","Month 3","Month 6"),],
               linewidth = 0.4, outlier.alpha = 0, show.legend = FALSE) +
  geom_point(data = rl.rt.means[rl.rt.means$sess %in% c("Baseline","Month 1","Month 3","Month 6"),],
             size = 0.8, color = "purple", alpha = .5,
             position = position_dodge2(width = .3), stroke = 0.6, aes(group = idx)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = element_blank(), y = "log sec", title = "Response time [RL]",
       color = element_blank(), shape = element_blank()) +
  scale_shape_manual(values = c(16, 4)) +
  stat_compare_means(label = "p.signif",
                     comparisons = list(c("Baseline","Month 1"), c("Baseline","Month 3"), c("Baseline","Month 6")),
                     method = "t.test", paired = TRUE,  # rows are sorted by idx in every session
                     bracket.size = 0.3, tip.length = 0.02, size = 2.5,
                     step.increase = 0.12)

# UG choice (panel = "UG_choice")
ug.beh.means <- source_data %>%
  filter(panel == "UG_choice") %>%
  mutate(idx  = factor(idx),
         sess = factor(sess, levels = sess_levels))

ug.c.sess <-
  ug.beh.means %>%
  filter(sess %in% c("Baseline","DBS","Month 1","Month 3","Month 6")) %>%
  ggplot(aes(x = sess, y = value)) +
  theme_pubr(base_size = 7) +
  geom_vline(xintercept = "DBS", linewidth = 0.5, color = "grey70") +
  geom_boxplot(data = ug.beh.means[ug.beh.means$sess %in% c("Baseline","DBS","Month 1","Month 3","Month 6"),],
               linewidth = 0.4, outlier.alpha = 0, show.legend = FALSE) +
  geom_point(data = ug.beh.means[ug.beh.means$sess %in% c("Baseline","Month 1","Month 3","Month 6"),],
             size = 0.8, color = "purple", alpha = .5,
             position = position_dodge2(width = .3), stroke = 0.6, aes(group = idx)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = element_blank(), y = "p(accept)", title = "Choice [UG]",
       color = element_blank(), shape = element_blank()) +
  scale_shape_manual(values = c(16, 4)) +
  stat_compare_means(label = "p.signif",
                     comparisons = list(c("Baseline","Month 1"), c("Baseline","Month 3"), c("Baseline","Month 6")),
                     method = "t.test", paired = TRUE,  # rows are sorted by idx in every session
                     bracket.size = 0.3, tip.length = 0.02, size = 2.5,
                     step.increase = 0.12)

# UG RT (panel = "UG_RT")
ug.rt.means <- source_data %>%
  filter(panel == "UG_RT") %>%
  mutate(idx  = factor(idx),
         sess = factor(sess, levels = sess_levels))

ug.logrt.sess <-
  ug.rt.means %>%
  filter(sess %in% c("Baseline","DBS","Month 1","Month 3","Month 6")) %>%
  ggplot(aes(x = sess, y = value)) +
  theme_pubr(base_size = 7) +
  geom_vline(xintercept = "DBS", linewidth = 0.5, color = "grey70") +
  geom_boxplot(data = ug.rt.means[ug.rt.means$sess %in% c("Baseline","DBS","Month 1","Month 3","Month 6"),],
               linewidth = 0.4, outlier.alpha = 0, show.legend = FALSE) +
  geom_point(data = ug.rt.means[ug.rt.means$sess %in% c("Baseline","Month 1","Month 3","Month 6"),],
             size = 0.8, color = "purple", alpha = .5,
             position = position_dodge2(width = .3), stroke = 0.6, aes(group = idx)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = element_blank(), y = "log sec", title = "Response time [UG]",
       color = element_blank(), shape = element_blank()) +
  scale_shape_manual(values = c(16, 4)) +
  stat_compare_means(label = "p.signif",
                     comparisons = list(c("Baseline","Month 1"), c("Baseline","Month 3"), c("Baseline","Month 6")),
                     method = "t.test", paired = TRUE,  # rows are sorted by idx in every session
                     bracket.size = 0.3, tip.length = 0.02, size = 2.5,
                     step.increase = 0.12)

# UG mood (panel = "UG_mood")
ug.mood.means <- source_data %>%
  filter(panel == "UG_mood") %>%
  mutate(idx  = factor(idx),
         sess = factor(sess, levels = sess_levels))

ug.mood.sess <-
  ug.mood.means %>%
  filter(sess %in% c("Baseline","DBS","Month 1","Month 3","Month 6")) %>%
  ggplot(aes(x = sess, y = value)) +
  theme_pubr(base_size = 7) +
  geom_vline(xintercept = "DBS", linewidth = 0.5, color = "grey70") +
  geom_boxplot(data = ug.mood.means[ug.mood.means$sess %in% c("Baseline","DBS","Month 1","Month 3","Month 6"),],
               linewidth = 0.4, outlier.alpha = 0, show.legend = FALSE) +
  geom_point(data = ug.mood.means[ug.mood.means$sess %in% c("Baseline","Month 1","Month 3","Month 6"),],
             size = 0.8, color = "purple", alpha = .5,
             position = position_dodge2(width = .3), stroke = 0.6, aes(group = idx)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = element_blank(), y = "rating", title = "Mood [UG]",
       color = element_blank(), shape = element_blank()) +
  scale_shape_manual(values = c(16, 4)) +
  stat_compare_means(label = "p.signif",
                     comparisons = list(c("Baseline","Month 1"), c("Baseline","Month 3"), c("Baseline","Month 6")),
                     method = "t.test", paired = TRUE,  # rows are sorted by idx in every session
                     bracket.size = 0.3, tip.length = 0.02, size = 2.5,
                     step.increase = 0.12)

figure3_top <- rl.logrt.sess + rl.c.sess + ug.logrt.sess + ug.c.sess + ug.mood.sess +
  plot_layout(guides = "collect", nrow = 1) & theme(legend.position = "none") &
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))  # headroom so brackets are not clipped

# ---- Bottom panels: correlation panels (F-G) --------------------------------

# Panel F: corr_5HT_RL
rl.rt.SE <- source_data %>%
  filter(panel == "corr_5HT_RL") %>%
  mutate(month = factor(month, levels = c("PostStim","W1","M1","M2","M3","M4","M5","M6"),
                        labels = c("Post-Stim","Week 1","Month 1","Month 2","Month 3",
                                   "Month 4","Month 5","Month 6")))

fig.rl.rt.SE.simple <-
  rl.rt.SE %>%
  ggplot(aes(x = nt_change, y = beh_change)) +
  theme_pubr(base_size = 7) +
  facet_wrap(~month, nrow = 1) +
  stat_smooth(method = "lm", alpha = .25, linetype = "solid", linewidth = 0.6,
              color = NT_colors$color[NT_colors$id == "5-HT"],
              fill  = NT_colors$color[NT_colors$id == "5-HT"]) +
  geom_point(size = 1, stroke = 0.5) +
  labs(title = "Response time [RL]",
       x = "5-HT Change\nPost - Pre [Z]",
       y = "Change from Baseline",
       shape = element_blank(), color = element_blank()) +
  guides(color = "none", shape = "none") +
  coord_cartesian(xlim = c(-12, 12), ylim = c(-3.25, 1.5)) +
  theme(panel.spacing = unit(1, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA)) +
  stat_cor(method = "spearman",
           aes(label = paste0(cut(..p..,
                                  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                  labels = c("'****'","'***'","'**'","'*'","'ns'")))),
           label.x = 0, size = 3, label.y = 1)

# Panel G: corr_DA_UG
ug.mood.DA <- source_data %>%
  filter(panel == "corr_DA_UG") %>%
  mutate(month = factor(month, levels = c("PostStim","W1","M1","M2","M3","M4","M5","M6"),
                        labels = c("Post-Stim","Week 1","Month 1","Month 2","Month 3",
                                   "Month 4","Month 5","Month 6")))

fig.ug.m.DA.simple <-
  ug.mood.DA %>%
  ggplot(aes(x = nt_change, y = beh_change)) +
  theme_pubr(base_size = 7) +
  facet_wrap(~month, nrow = 1) +
  stat_smooth(method = "lm", alpha = .25, linetype = "solid", linewidth = 0.6,
              color = NT_colors$color[NT_colors$id == "DA"],
              fill  = NT_colors$color[NT_colors$id == "DA"]) +
  geom_point(size = 1, stroke = 0.5) +
  labs(title = "Mood [UG]",
       x = "DA Change\nPost - Pre [Z]",
       y = "Change from Baseline",
       shape = element_blank(), color = element_blank()) +
  coord_cartesian(xlim = c(-12, 12), ylim = c(-55, 55)) +
  theme(panel.spacing = unit(1, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA)) +
  stat_cor(method = "spearman",
           aes(label = paste0(cut(..p..,
                                  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                  labels = c("'****'","'***'","'**'","'*'","'ns'")))),
           label.x = -5, size = 3, label.y = 43)

figure3_bot <- (fig.rl.rt.SE.simple + fig.ug.m.DA.simple) +
  plot_layout(guides = "collect", nrow = 1) &
  theme(legend.position = 'right')

# Nature double column: 183 mm wide, 5-7 pt text, 8 pt bold panel labels
figure3 <- figure3_top / figure3_bot +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 7))

ggsave(file.path(res_dir, "figure3.pdf"), plot = figure3, device = cairo_pdf,
       width = 183, height = 120, units = "mm")
ggsave(file.path(res_dir, "figure3.png"), plot = figure3,
       width = 183, height = 120, units = "mm", dpi = 300)
