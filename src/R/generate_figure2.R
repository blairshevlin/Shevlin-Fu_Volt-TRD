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

fig.ug.trial = 
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
  ggplot(   aes(x = trial, y = m, fill = color, group = interaction(nt,stim))) +
  theme_pubr(base_size = 14) +
  facet_wrap(~nt) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  stat_summary( aes(color = color),geom = "line", size = 1, linewidth = 1.5) +
  stat_summary(
    geom = "ribbon", alpha = 0.35,show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_identity() +
  # Common styling
  labs(y = "Estimate [z]", x = "Trial", title = "Ultimatum game") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        strip.text = element_text(color = "black"),
        panel.spacing = unit(0.2, "cm")) +
  coord_cartesian(ylim = c(-7,7)) +
  guides(color = "none",
         fill = "none",
  ) +
  theme(legend.key = element_rect(color = NA)) 

fig.rl.trial = 
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
  ggplot(   aes(x = trial, y = m, fill = color, group = interaction(nt,stim))) +
  theme_pubr(base_size = 14) +
  facet_wrap(~nt) +
  scale_color_identity() +
  scale_fill_identity() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  stat_summary( aes(color = color),geom = "line", size = 1, linewidth = 1.5) +
  stat_summary(
    geom = "ribbon", alpha = 0.35,
    show.legend = FALSE
  ) +
  # Common styling
  labs(y = "Estimate [z]", x = "Trial", title = "Reversal learning") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        strip.text = element_text(color = "black"),
        panel.spacing = unit(0.2, "cm")) +
  coord_cartesian(ylim = c(-7,7)) +
  guides(color = "none",
         fill = "none",
  ) +
  theme(legend.key = element_rect(color = NA)) 

# Create participant avg plots
rl.avg =   rl.trial %>%
  filter(nt_metric == "Oz") %>%
  group_by(idx,stim,nt) %>%
  summarise(m = mean(mTrial)) %>%
  mutate(nt = factor(nt, levels = c("DA","SE"), labels = c("DA","5-HT")),
         stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         color = ifelse(stim == "Pre-Stim",
                        ifelse(nt == "5-HT","#9ecae1","#fcbba1"),
                        ifelse(nt == "5-HT","#08519c", "#a50f15"))
  ) %>% as.data.frame()

ug.avg =   ug.trial %>%
  filter(nt_metric == "Oz") %>%
  group_by(idx,stim,nt) %>%
  summarise(m = mean(mTrial)) %>%
  mutate(nt = factor(nt, levels = c("DA","SE"), labels = c("DA","5-HT")),
         stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         color = ifelse(stim == "Pre-Stim",
                        ifelse(nt == "5-HT","#9ecae1","#fcbba1"),
                        ifelse(nt == "5-HT","#08519c", "#a50f15"))
  ) %>% as.data.frame()

offset_data_RL <- rl.avg %>%
  mutate(
    # For points: Baseline to right, Month 6 to left
    point_x_pos = ifelse(stim == "Pre-Stim", 1 + 0.1, 2 - 0.1),
    # For boxplots: Baseline to left, Month 6 to right (opposite of points)
    box_x_pos = ifelse(stim == "Pre-Stim", 1 - 0.1, 2 + 0.1),
    # Convert sess_fig to numeric for proper plotting
    sess_num = ifelse(stim == "Pre-Stim", 1, 2)
  )
offset_data_UG <- ug.avg %>%
  mutate(
    # For points: Baseline to right, Month 6 to left
    point_x_pos = ifelse(stim == "Pre-Stim", 1 + 0.1, 2 - 0.1),
    # For boxplots: Baseline to left, Month 6 to right (opposite of points)
    box_x_pos = ifelse(stim == "Pre-Stim", 1 - 0.1, 2 + 0.1),
    # Convert sess_fig to numeric for proper plotting
    sess_num = ifelse(stim == "Pre-Stim", 1, 2)
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

fig.rl.avg = 
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  geom_line(data = offset_data_RL,
    aes(x = point_x_pos,y = m, group = idx),
    position = position_dodge2(width = .2),
    size = 1.5, alpha = .15
  ) +
  geom_point(data = offset_data_RL, 
    aes(x = point_x_pos, y = m, group = idx, color = color),
    position = position_dodge2(width = .2), 
    shape = 1,
    stroke = 2,
    alpha = 0.5,
    size = 3
  ) +
  stat_summary(geom = "line", data = offset_data_RL,
    aes(x = box_x_pos, y = m, group = nt),
    linewidth = 1.1)+
  stat_summary(data = offset_data_RL,
    aes(x = box_x_pos,y = m, group = stim, color = color),
    size = 1.1,
    #shape = 15,
    linewidth = 1.1) +
  facet_wrap(~ nt) +
  scale_color_identity() +
  labs(y = "Estimate [z]", x = element_blank(), 
       title = "Reversal learning") +
  theme_pubr(base_size = 14) +
    scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Pre-Stim", "Post-Stim"),
    limits = c(0.5, 2.5)
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 14),
    panel.spacing = unit(0.2, "cm"),
    axis.text.x = element_text(angle = 15, vjust = 0.95, hjust=0.5)
  )

fig.ug.avg = 
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  geom_line(data = offset_data_UG,
    aes(x = point_x_pos,y = m, group = idx),
    position = position_dodge2(width = .2), 
    size = 1.5, 
    alpha = .15) +
  geom_point(data = offset_data_UG, 
    aes(x = point_x_pos, y = m, group = idx, color = color),
    position = position_dodge2(width = .2), 
    shape = 1,
    stroke = 2,
    alpha = 0.5,
    size = 3) +
  stat_summary(geom = "line", data = offset_data_UG,
    aes(x = box_x_pos, y = m, group = nt),
    linewidth = 1.1
  )+
  stat_summary(data = offset_data_UG,aes(x = box_x_pos,y = m, group = stim, color = color), 
  size = 1.1 ,linewidth = 1.1, 
  #shape = 15
  ) +
  facet_wrap(~ nt) +
  scale_color_identity() +
  labs(y = "Estimate [z]", x = element_blank(), 
       title = "Ultimatum game") +
  theme_pubr(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 14),
    panel.spacing = unit(0.2, "cm"),
    axis.text.x = element_text(angle = 15, vjust = 0.95, hjust=0.5)
  )+
    scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Pre-Stim", "Post-Stim"),
    limits = c(0.5, 2.5)
  )

# Task feature data

offset_data_UG_task <- ug.trial %>%
  filter(nt_metric == "Oz") %>%
  group_by(idx,stim,nt,offer_bin) %>%
  summarise(m = mean(mTrial)) %>%
  mutate(nt = factor(nt, levels = c("DA","SE"), labels = c("DA","5-HT")),
         stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         color = ifelse(stim == "Pre-Stim",
                        ifelse(nt == "5-HT","#9ecae1", "#fcbba1"),
                        ifelse(nt == "5-HT","#08519c", "#a50f15")),
         # Create numeric positions for offsetting
         offer_num = as.numeric(offer_bin),
         # Points offset inward, group means offset outward
        # point_x_pos = ifelse(stim == "Pre-Stim", offer_num - 0.1, offer_num + 0.1),
         #group_x_pos = ifelse(stim == "Pre-Stim", offer_num - 0.3, offer_num + 0.3)
         point_x_pos = ifelse(stim == "Pre-Stim", offer_num - 0.15, offer_num + 0.15),
         group_x_pos = ifelse(stim == "Pre-Stim", offer_num - 0.35, offer_num + 0.35)
  )

  offset_data_RL_task <- rl.trial %>%
  filter(nt_metric == "Oz") %>%
  group_by(idx,stim,nt,block_type) %>%
  summarise(m = mean(mTrial)) %>%
  mutate(nt = factor(nt, levels = c("DA","SE"), labels = c("DA","5-HT")),
         stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         block_type = factor(block_type,
                             levels = c("Positive","Mixed","Negative"),
                             labels = c("Reward","Mixed","Punishment")),
         color = ifelse(stim == "Pre-Stim",
                        ifelse(nt == "5-HT","#9ecae1","#fcbba1"),
                        ifelse(nt == "5-HT","#08519c", "#a50f15")),
         # Create numeric positions for offsetting
         block_num = as.numeric(block_type),
         # Points offset inward, group means offset outward
         point_x_pos = ifelse(stim == "Pre-Stim", block_num - 0.15, block_num + 0.15),
         group_x_pos = ifelse(stim == "Pre-Stim", block_num - 0.35, block_num + 0.35) #was 0.3
  )


fig.rl.task = 
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  theme_pubr(base_size = 14) +
  facet_wrap(~nt) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  # Individual points with inward offset
  geom_point(data = offset_data_RL_task,
             aes(x = point_x_pos, y = m, color = color),
             alpha = .75, size =  1.5, stroke = 1.5, shape = 1) + #size was 2
  # Group means (diamonds) with outward offset
  stat_summary(data = offset_data_RL_task,
               aes(x = group_x_pos, y = m, color = color, group = stim),
               size = 1.25, linewidth = 1.25,
               position = position_identity()) +
  scale_color_identity() +
  scale_fill_identity() +
  # Use original block_type labels for x-axis
  scale_x_continuous(
    breaks = c(1, 2, 3),
    labels = c("Reward", "Mixed", "Punishment"),
    limits = c(0.5, 3.5)
  ) +
  labs(y = "Estimate [z]", x = "Block type", title = "Reversal learning") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        strip.text = element_text(color = "black"),
        panel.spacing = unit(0.2, "cm"),
        axis.text.x = element_text(angle = 15, vjust = 0.9, hjust=0.6)) +
  coord_cartesian(ylim = c(-7.5,7.5)) +
  guides(color = "none", fill = "none") +
  theme(legend.key = element_rect(color = NA))

fig.ug.task = 
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  theme_pubr(base_size = 14) +
  facet_wrap(~nt) +
  geom_point(data = offset_data_UG_task, 
    aes(x = point_x_pos, y = m, group = idx, color = color),
    shape = 1,
    stroke = 1.75,
    alpha = 0.75,
    size = 1.5) +
  stat_summary(data = offset_data_UG_task,
  aes(x = group_x_pos,y = m, group = stim, color = color), 
  size = 1.25 ,linewidth = 1.25, 
  #shape = 15
  ) +
  scale_color_identity() +
  scale_fill_identity() +
  # Use original offer_bin labels for x-axis
  scale_x_continuous(
    breaks = c(1, 2, 3),
    labels = c("$1-3", "$4-6", "$7-9"),
    limits = c(0.5, 3.5)
  ) +
  labs(y = "Estimate [z]", x = "Offer", title = "Ultimatum game") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        strip.text = element_text(color = "black"),
        panel.spacing = unit(0.2, "cm")) +
  coord_cartesian(ylim = c(-7.5,7.5)) +
  guides(color = "none", fill = "none") +
  theme(legend.key = element_rect(color = NA))

main_plot <- 
             (fig.rl.trial | fig.ug.trial) / 
             (fig.rl.avg | fig.ug.avg) / 
             (fig.rl.task | fig.ug.task) + 
  plot_annotation(tag_levels = "a")

final_plot <- plot_grid(main_plot, stacked_legends, 
                        ncol = 2, rel_widths = c(4, 1),
                        align = "h", axis = "tb")

ggsave(res_dir / "figure2.png", 
       plot = final_plot,
       device = "png",
       width = 12,          # Width in inches
       height = 10,         # Height in inches  
       dpi = 300)           # Resolution (300 DPI for publication quality)


panel_a <- plot_grid((fig.rl.trial | fig.ug.trial)  + 
  plot_annotation(tag_levels = "a"), stacked_legends, 
                        ncol = 2, rel_widths = c(4, 1),
                        align = "h", axis = "tb")

ggsave(res_dir / "figure2_panelA.png", 
       plot = panel_a,
       device = "png",
       width = 12,          # Width in inches
       height = 4,         # Height in inches  
       dpi = 300)           # Resolution (300 DPI for publication quality)
