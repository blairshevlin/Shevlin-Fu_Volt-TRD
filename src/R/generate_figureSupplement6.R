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
# 2025/11/14    Blair Shevlin                           wrote original code

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
dir = path(here())
nt_dir = dir / "data" / "nt" / "processed"
res_dir = dir / "results" # Updated

# Subject-level NT Data
load(nt_dir / "UG_RL_NT-Continuous_9-23-25.RData")

# Generate summary data - focusing only on reward events (RL) and offer events (UG)
rl.EST.Reward = rl.EST %>%
  filter(event == "Reward") %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim"))) %>%
  select(!c(O,R,P,M,O10,Oz10,R10,Rz10,P10,Pz10,M10,Mz10))

ug.EST.Offer = ug.EST %>%
  filter(event == "Offer") %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim"))) %>%
  select(!c(O,R,P,M,O10,Oz10,R10,Rz10,P10,Pz10,M10,Mz10))

# IDs
ids_final = unique(ug.EST.Offer$idx)

# Specify colors
NT_colors = data.frame(id = c("DA","5-HT"),
                       color = c("#cb181d","#2171b5"))

ug.trial = ug.EST.Offer %>%
  filter(nt %in% c("DA","SE")) %>%
  pivot_longer(cols = c("Oz","Rz","Pz","Mz"),
               names_to = "nt_metric", values_to = "nt_val") %>%
  mutate(offer_bin = ifelse(offer < 4, "$1-3",
                            ifelse(offer < 7, "$4-6","$7-9")),
         offer_bin = factor(offer_bin,
                            levels = c("$1-3","$4-6","$7-9"))) %>%
  group_by(idx,stim,nt,nt_metric,trial,offer,offer_bin) %>%
  summarise(mTrial = mean(nt_val)) %>% ungroup() %>%
  mutate(offerz = scale(offer)[,1])

rl.trial = rl.EST.Reward %>%
  filter(nt %in% c("DA","SE")) %>%
  pivot_longer(cols = c("Oz","Rz","Pz","Mz"),
               names_to = "nt_metric", values_to = "nt_val") %>%
  mutate(block_type = factor(cond, levels = c("Mixed","Negative","Positive"))
  ) %>%
  group_by(idx,stim,nt,nt_metric,trial,block_type,trial_within_block,prev_rew_raw) %>%
  summarise(mTrial = mean(nt_val),
            prev_rew = ifelse(prev_rew_raw == 0,
                              ifelse(block_type == "Negative" | block_type == "Mixed",-10,prev_rew_raw
                              ),
                              ifelse(prev_rew_raw == 1, 
                                     ifelse(block_type == "Positive" | block_type == "Mixed",10, 0),
                                     prev_rew_raw
                              )
            )
  ) %>% ungroup()


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
      legend.margin = margin(10, 0, 0, 0)
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

# Stack the legends vertically
stacked_legends <- plot_grid(prestim_legend, poststim_legend, 
                             ncol = 1, align = "v")


fig.ug.trial.byparticipant = 
  ug.trial %>%
  filter(nt_metric == "Oz") %>%
  group_by(idx,stim,nt,trial) %>%
  summarise(m = mean(mTrial)) %>%
  mutate(nt = factor(nt, levels = c("DA","SE"), labels = c("DA","5-HT")),
         stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         color = ifelse(stim == "Pre-Stim",
                        ifelse(nt == "5-HT","#9ecae1", "#fcbba1"),
                        ifelse(nt == "5-HT","#08519c", "#a50f15"))
  ) %>%
  ggplot(aes(x = trial, y = m, fill = color, group = interaction(nt,stim))) +
  theme_pubr(base_size = 10) +
  facet_grid(idx ~ nt) +  # Facet by participant (idx) and neurotransmitter
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "gray") +
  stat_summary(aes(color = color), geom = "line", size = 0.8, linewidth = 1) +
  stat_summary(geom = "ribbon", alpha = 0.35, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_identity() +
  labs(y = "Estimate [z]", x = "Trial", title = "Ultimatum game") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = element_text(color = "black", size = 8),
        panel.spacing = unit(0.1, "cm")) +
  coord_cartesian(ylim = c(-10,10)) +
  guides(color = "none", fill = "none") +
  theme(legend.key = element_rect(color = NA))

# Modified RL trial plot with participant faceting
fig.rl.trial.byparticipant = 
  rl.trial %>%
  filter(nt_metric == "Oz") %>%
  group_by(idx,stim,nt,trial) %>%
  summarise(m = mean(mTrial)) %>%
  mutate(nt = factor(nt, levels = c("DA","SE","NE"), labels = c("DA","5-HT","NE")),
         stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         color = ifelse(stim == "Pre-Stim",
                        ifelse(nt == "5-HT","#9ecae1", "#fcbba1"),
                        ifelse(nt == "5-HT","#08519c","#a50f15"))
  ) %>%
  ggplot(aes(x = trial, y = m, fill = color, group = interaction(nt,stim))) +
  theme_pubr(base_size = 10) +
  facet_grid(idx ~ nt) +  # Facet by participant (idx) and neurotransmitter
  scale_color_identity() +
  scale_fill_identity() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "gray") +
  stat_summary(aes(color = color), geom = "line", size = 0.8, linewidth = 1) +
  stat_summary(geom = "ribbon", alpha = 0.35, show.legend = FALSE) +
  labs(y = "Estimate [z]", x = "Trial", title = "Reversal learning") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = element_text(color = "black", size = 8),
        panel.spacing = unit(0.1, "cm")) +
  coord_cartesian(ylim = c(-10,10)) +
  guides(color = "none", fill = "none") +
  theme(legend.key = element_rect(color = NA))


# Create combined plot with participant-level figures
participant_plot <- (fig.rl.trial.byparticipant | fig.ug.trial.byparticipant) + 
  plot_annotation(tag_levels = "a")

# Combine with legend
final_participant_plot <- plot_grid(participant_plot, stacked_legends, 
                                    ncol = 2, rel_widths = c(5, 1),
                                    align = "h", axis = "tb")
ggsave(res_dir / "figureS6.png", 
       plot = final_participant_plot,
       device = "png",
       width = 16,          # Width in inches
       height = 12,         # Height in inches  
       dpi = 300)           # Resolution (300 DPI for publication quality)
