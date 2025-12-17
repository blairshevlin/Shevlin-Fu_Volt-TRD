
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

# Paths
dir = path(here())
nt_dir = dir / "data" / "nt" / "processed"
res_dir = dir / "results" # Updated
clin_dir = dir / "data" / "clinical"

# Subject-level NT Data
load(nt_dir / "UG_RL_NT-Continuous_Raw-Model-Output_11-25-25.RData")

# Generate summary data - focusing only on reward events (RL) and offer events (UG)
rl.EST.Reward = rl.EST %>%
  filter(event == "Reward") %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim"))) %>%
  select(!c(R,P,M,O10,Oz10,R10,Rz10,P10,Pz10,M10,Mz10))

ug.EST.Offer = ug.EST %>%
  filter(event == "Offer") %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim"))) %>%
  select(!c(R,P,M,O10,Oz10,R10,Rz10,P10,Pz10,M10,Mz10))

# IDs
ids_final = unique(ug.EST.Offer$idx)

# Specify colors
NT_colors = data.frame(id = c("DA","5-HT"),
                       color = c("#cb181d","#2171b5"))


ug.trial = ug.EST.Offer %>%
  filter(nt %in% c("DA","SE")) %>%
  pivot_longer(cols = c("O","Oz"),
               names_to = "nt_metric", values_to = "nt_val") %>%
  mutate(offer_bin = ifelse(offer < 4, "$1-3",
                            ifelse(offer < 7, "$4-6","$7-9")),
         offer_bin = factor(offer_bin,
                            levels = c("$1-3","$4-6","$7-9"))) %>%
  group_by(idx,stim,nt,nt_metric,trial,offer,offer_bin) %>%
  reframe(mTrial = mean(nt_val)) %>% ungroup() %>%
  mutate(offerz = scale(offer)[,1]) %>%
  filter(nt_metric == "O") %>%
  group_by(idx,stim,nt) %>%
  summarise(m = mean(mTrial), s = se(mTrial)) %>%
  mutate(nt = factor(nt, levels = c("DA","SE"), labels = c("DA","5-HT")),
         stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         color = ifelse(stim == "Pre-Stim",
                        ifelse(nt == "5-HT","#9ecae1","#fcbba1"),
                        ifelse(nt == "5-HT","#08519c", "#a50f15"))
  ) %>% as.data.frame()

rl.trial = rl.EST.Reward %>%
  filter(nt %in% c("DA","SE")) %>%
  pivot_longer(cols = c("Oz","O"),
               names_to = "nt_metric", values_to = "nt_val") %>%
  mutate(block_type = factor(cond, levels = c("Mixed","Negative","Positive"))
  ) %>%
  group_by(idx,stim,nt,nt_metric,trial,block_type,trial_within_block,prev_rew_raw) %>%
  reframe(mTrial = mean(nt_val),
            prev_rew = ifelse(prev_rew_raw == 0,
                              ifelse(block_type == "Negative" | block_type == "Mixed",-10,prev_rew_raw
                              ),
                              ifelse(prev_rew_raw == 1, 
                                     ifelse(block_type == "Positive" | block_type == "Mixed",10, 0),
                                     prev_rew_raw
                              )
            )
  ) %>% ungroup() %>%
  filter(nt_metric == "O") %>%
    group_by(idx,stim,nt) %>%
  summarise(m = mean(mTrial), s = se(mTrial)) %>%
  mutate(nt = factor(nt, levels = c("DA","SE"), labels = c("DA","5-HT")),
         stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         color = ifelse(stim == "Pre-Stim",
                        ifelse(nt == "5-HT","#9ecae1","#fcbba1"),
                        ifelse(nt == "5-HT","#08519c", "#a50f15"))
  ) %>% as.data.frame()


# Generate raw data summaries
rl.raw = rl.EST.Reward %>%
  filter(nt %in% c("DA","SE")) %>%
  pivot_longer(cols = c("Oz","O"),
               names_to = "nt_metric", values_to = "nt_val") %>%
  mutate(block_type = factor(cond, levels = c("Mixed","Negative","Positive")))

ug.raw = ug.EST.Offer %>%
  filter(nt %in% c("DA","SE")) %>%
  pivot_longer(cols = c("Oz","O"),
               names_to = "nt_metric", values_to = "nt_val") %>%
  mutate(offer_bin = ifelse(offer < 4, "$1-3",
                            ifelse(offer < 7, "$4-6","$7-9")),
         offer_bin = factor(offer_bin, levels = c("$1-3","$4-6","$7-9")))

# Create subject-level averages with error bars (inter-trial variability)
rl.subj_summary = rl.raw %>%
  filter(nt_metric == "O") %>%
  group_by(idx, stim, nt) %>%
  summarise(
    mean_nt = mean(nt_val, na.rm = TRUE),
    se_nt = se(nt_val),
    sd_nt = sd(nt_val, na.rm = TRUE),
    n_trials = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    nt = factor(nt, levels = c("DA","SE"), labels = c("DA","5-HT")),
    stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
    color = ifelse(stim == "Pre-Stim",
                   ifelse(nt == "5-HT","#9ecae1","#fcbba1"),
                   ifelse(nt == "5-HT","#08519c", "#a50f15"))
  )

ug.subj_summary = ug.raw %>%
  filter(nt_metric == "O") %>%
  group_by(idx, stim, nt) %>%
  summarise(
    mean_nt = mean(nt_val, na.rm = TRUE),
    se_nt = se(nt_val),
    sd_nt = sd(nt_val, na.rm = TRUE),
    n_trials = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    nt = factor(nt, levels = c("DA","SE"), labels = c("DA","5-HT")),
    stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
    color = ifelse(stim == "Pre-Stim",
                   ifelse(nt == "5-HT","#9ecae1","#fcbba1"),
                   ifelse(nt == "5-HT","#08519c", "#a50f15"))
  )

# Create offset data for proper positioning
offset_data_RL_raw <- rl.subj_summary %>%
  mutate(
    point_x_pos = ifelse(stim == "Pre-Stim", 1 + 0.1, 2 - 0.1),
    box_x_pos = ifelse(stim == "Pre-Stim", 1 - 0.1, 2 + 0.1)
  )

offset_data_UG_raw <- ug.subj_summary %>%
  mutate(
    point_x_pos = ifelse(stim == "Pre-Stim", 1 + 0.1, 2 - 0.1),
    box_x_pos = ifelse(stim == "Pre-Stim", 1 - 0.1, 2 + 0.1)
  )

# Create the main plots
fig.rl.raw = 
  ggplot() +
  facet_wrap(~nt) +
  #coord_cartesian(ylim=c(45,70))+
  theme_pubr(base_size = 14) +
  geom_errorbar(
    data = offset_data_RL_raw,
    aes(x = point_x_pos, y = mean_nt, 
        ymin = mean_nt - se_nt, ymax = mean_nt + se_nt,
        color = color),
    width = 0.05, alpha = 0.6, linewidth = 0.8
  ) +
  geom_line(
    data = offset_data_RL_raw,
    aes(x = point_x_pos, y = mean_nt, group = idx),
    linewidth = 1, alpha = 0.4, color = "gray50"
  ) +
  geom_point(
    data = offset_data_RL_raw,
    aes(x = point_x_pos, y = mean_nt, color = color),
    size = 2.5, alpha = 0.7
  ) +
  stat_summary(
    data = offset_data_RL_raw,
    aes(x = box_x_pos, y = mean_nt, color = color),
    fun = mean,
    fun.min = function(x) mean(x) - se(x),
    fun.max = function(x) mean(x) + se(x),
    geom = "errorbar",
    linewidth = 2, width = 0.1, alpha = 1
  ) +
  stat_summary(
    data = offset_data_RL_raw,
    aes(x = box_x_pos, y = mean_nt, color = color),
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
    panel.border = element_rect(colour = NA, fill = NA),
    strip.text = element_text(color = "black"),
    panel.spacing = unit(0.2, "cm"),
    axis.text.x = element_text(angle = 15, vjust = 1, hjust = 0.5)
  )

fig.ug.raw = 
ggplot() +
  facet_wrap(~nt) +
  theme_pubr(base_size = 14) +
  #coord_cartesian(ylim=c(45,70))+
  geom_errorbar(
    data = offset_data_UG_raw,
    aes(x = point_x_pos, y = mean_nt, 
        ymin = mean_nt - se_nt, ymax = mean_nt + se_nt,
        color = color),
    width = 0.05, alpha = 0.6, linewidth = 0.8
  ) +
  geom_line(
    data = offset_data_UG_raw,
    aes(x = point_x_pos, y = mean_nt, group = idx),
    linewidth = 1, alpha = 0.4, color = "gray50"
  ) +
  geom_point(
    data = offset_data_UG_raw,
    aes(x = point_x_pos, y = mean_nt, color = color),
    size = 2.5, alpha = 0.7
  ) +
  stat_summary(
    data = offset_data_UG_raw,
    aes(x = box_x_pos, y = mean_nt, color = color),
    fun = mean,
    fun.min = function(x) mean(x) - se(x),
    fun.max = function(x) mean(x) + se(x),
    geom = "errorbar",
    linewidth = 2, width = 0.1, alpha = 1
  ) +
  stat_summary(
    data = offset_data_UG_raw,
    aes(x = box_x_pos, y = mean_nt, color = color),
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
    panel.border = element_rect(colour = NA, fill = NA),
    strip.text = element_text(color = "black"),
    panel.spacing = unit(0.2, "cm"),
    axis.text.x = element_text(angle = 15, vjust = 1, hjust = 0.5)
  )

create_prestim_legend <- function() {
  # Dummy data for Pre-Stim legend
  legend_data_pre <- data.frame(
    nt = factor(c("DA", "5-HT"), levels = c("DA", "5-HT")),
    color = c("#fcbba1", "#9ecae1"),
    x = 1,
    y = 1
  )
  
  ggplot(legend_data_pre, aes(x = x, y = y, color = nt)) +
    geom_point(size = 0) +  # Invisible points
    geom_line(size = 2) +   # This creates the line legend
    scale_color_manual(
      name = "Pre-Stim",
      values = c("DA" = "#fcbba1", "5-HT" = "#9ecae1")
    ) +
    theme_void() +
    theme(
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 10),
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.margin = margin(0, 0, 10, 0)
    ) +
    guides(color = guide_legend(
      override.aes = list(size = 2, alpha = 1),
      title.position = "top",
      title.hjust = 0.5
    ))
}

create_poststim_legend <- function() {
  # Dummy data for Post-Stim legend
  legend_data_post <- data.frame(
    nt = factor(c("DA", "5-HT"), levels = c("DA", "5-HT")),
    color = c("#a50f15", "#08519c"),
    x = 1,
    y = 1
  )
  
  ggplot(legend_data_post, aes(x = x, y = y, color = nt)) +
    geom_point(size = 0) +  # Invisible points
    geom_line(size = 2) +   # This creates the line legend
    scale_color_manual(
      name = "Post-Stim",
      values = c("DA" = "#a50f15", "5-HT" = "#08519c")
    ) +
    theme_void() +
    theme(
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 10),
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.margin = margin(10, 0, 10, 0)
    ) +
    guides(color = guide_legend(
      override.aes = list(size = 2, alpha = 1),
      title.position = "top",
      title.hjust = 0.5
    ))
}

# Extract the legends
prestim_legend <- get_legend(create_prestim_legend())
poststim_legend <- get_legend(create_poststim_legend())

stacked_legends <- 
plot_grid(prestim_legend, poststim_legend, 
                             ncol = 1, align = "v")


# Print summary statistics for the manuscript
cat("Summary Statistics for Raw Data:\n")
cat("Reversal Learning:\n")
print(rl.subj_summary %>% 
      group_by(nt, stim) %>% 
      summarise(
        mean_estimate = mean(mean_nt, na.rm = TRUE),
        sd_estimate = sd(mean_nt, na.rm = TRUE),
        median_trials = median(n_trials, na.rm = TRUE),
        .groups = 'drop'
      ))

cat("\nUltimatum Game:\n")
print(ug.subj_summary %>% 
      group_by(nt, stim) %>% 
      summarise(
        mean_estimate = mean(mean_nt, na.rm = TRUE),
        sd_estimate = sd(mean_nt, na.rm = TRUE),
        median_trials = median(n_trials, na.rm = TRUE),
        .groups = 'drop'
      ))


# Change statistics
rl.wilcox <- rl.subj_summary %>%
  group_by(nt) %>%
  summarise(
    wilcox_test = list(wilcox.test(
      mean_nt[stim == "Pre-Stim"], 
      mean_nt[stim == "Post-Stim"], 
      paired = TRUE
    )),
    .groups = 'drop'
  ) %>%
  mutate(tidy_results = map(wilcox_test, tidy)) %>%
  unnest(tidy_results) %>%
  select(-wilcox_test)

ug.wilcox <- ug.subj_summary %>%
  group_by(nt) %>%
  summarise(
    wilcox_test = list(wilcox.test(
      mean_nt[stim == "Pre-Stim"], 
      mean_nt[stim == "Post-Stim"], 
      paired = TRUE
    )),
    .groups = 'drop'
  ) %>%
  mutate(tidy_results = map(wilcox_test, tidy)) %>%
  unnest(tidy_results) %>%
  select(-wilcox_test)

print(rl.wilcox)
print(ug.wilcox)

# Load clinical data
cl.df = read.csv(clin_dir / "clinical-data_deid_07-10-25.csv" )  %>%
  filter(session == "pre stim") %>%
  select(idx, HDRS) %>%
  mutate(HDRS_Z = scale(HDRS)[,1])

# Merge with NT data for correlation analyses
rl.cl.df = merge(rl.subj_summary, cl.df, by = "idx")
ug.cl.df = merge(ug.subj_summary, cl.df, by = "idx")  



# Run correlations and print results
cat("\nCorrelation Results:\n") 
for (nt_type in c("DA", "5-HT")) {
  for (task in c("Reversal Learning", "Ultimatum Game")) {
    if (task == "Reversal Learning") {
      data_subset = rl.cl.df %>% filter(nt == nt_type, stim == "Pre-Stim")
    } else {
      data_subset = ug.cl.df %>% filter(nt == nt_type, stim == "Pre-Stim")
    }
    
    cor_test = cor.test(data_subset$mean_nt, data_subset$HDRS, method = "spearman")
    
    cat(sprintf("%s - %s:\n", task, nt_type))
    cat(sprintf("  Pearson's r: %.3f\n", cor_test$estimate))
    cat(sprintf("  p-value: %.4f\n", cor_test$p.value))
    cat(sprintf("  n: %d\n\n", nrow(data_subset)))
  }
}

# Generate figures showing correlations
fig.rl.cor = 
ggplot(rl.cl.df %>% filter(stim == "Pre-Stim"), 
                     aes(x = HDRS, y = mean_nt, color = nt)) +
  geom_point(size = 6) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  facet_wrap(~nt) +
  theme_pubr(base_size = 14) +
  #coord_cartesian(ylim=c(-5000,8000))+
  scale_color_manual(values = c("DA" = "#fcbba1", "5-HT" = "#9ecae1")) +
  labs(x = "HDRS", y = "Pre-Stim\nEstimate [Z]", title = "Reversal learning") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = NA, fill = NA),
    strip.text = element_text(color = "black"),
    panel.spacing = unit(0.2, "cm"),
    legend.position = "none"
  ) +
  stat_cor(color = "black",method = "spearman", label.x.npc = "left", label.y = 200, size = 4)

fig.ug.cor = 
ggplot(ug.cl.df %>% filter(stim == "Pre-Stim"), 
                     aes(x = HDRS, y = mean_nt, color = nt)) +
  geom_point(size = 6) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  facet_wrap(~nt) +
  theme_pubr(base_size = 14) +
  #coord_cartesian(ylim=c(-5000,8000))+
  scale_color_manual(values = c("DA" = "#fcbba1", "5-HT" = "#9ecae1")) +
  labs(x = "HDRS", y = "Pre-Stim\nEstimate [Z]", title = "Ultimatum game") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = NA, fill = NA),
    strip.text = element_text(color = "black"),
    panel.spacing = unit(0.2, "cm"),
    legend.position = "none"
  ) +
  stat_cor(color = "black",method = "spearman", label.x.npc = "left", label.y = 70, size = 4)


# Combine plots
main_plot_raw <- (fig.rl.raw | fig.ug.raw ) / (fig.rl.cor | fig.ug.cor) + 
  plot_annotation(tag_levels = "a", 
                  theme = theme(plot.title = element_text(size = 16, face = "bold"),
                                plot.subtitle = element_text(size = 12)))

final_plot_raw <- plot_grid(main_plot_raw, stacked_legends, 
                            ncol = 2, rel_widths = c(4, 1),
                            align = "h", axis = "tb")

# Save the figure
ggsave(res_dir / "figureS1.png", 
       plot = final_plot_raw,
       device = "png",
       width = 14,          # Slightly wider for the distribution plots
       height = 10,         # Taller for the two-row layout
       dpi = 300)
