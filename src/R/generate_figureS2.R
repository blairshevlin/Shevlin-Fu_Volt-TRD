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

# Paths
dir = path(here())
nt_dir = dir / "data" / "nt" / "processed"
fig_dir = dir / "results" 

# Subject-level NT Data
load(nt_dir / "UG_RL_NT-Continuous_7-14-25.RData")


# IDs
ids_final = unique(ug.EST.Offer$idx)

# Specify colors
id_colors = data.frame(id = ids_final,
                       color = brewer.pal(n=10,name = "Paired"),
                       shape = c(0:9))

NT_colors = data.frame(id = c("DA","5-HT","NE"),
                       color = c("#cb181d","#2171b5","darkgreen"))


# Lets run some trial-level estimate using new measures of trial-level (from cue)
ug.trial = ug.EST.Offer %>%
  pivot_longer(cols = c("Oz","Rz","Pz","Mz","Totz"),
               names_to = "nt_metric", values_to = "nt_val") %>%
  mutate(offer_bin = ifelse(offer < 4, "$1-3",
                            ifelse(offer < 7, "$4-6","$7-9")),
         offer_bin = factor(offer_bin,
                            levels = c("$1-3","$4-6","$7-9"))) %>%
  group_by(idx,stim,nt,nt_metric,trial,offer,offer_bin) %>%
  summarise(mTrial = mean(nt_val)) %>% ungroup() %>%
  mutate(offerz = scale(offer)[,1])

rl.trial = rl.EST.Reward %>%
  pivot_longer(cols = c("Oz","Rz","Pz","Mz","Totz"),
               names_to = "nt_metric", values_to = "nt_val") %>%
  mutate(block_type = factor(cond, levels = c("Mixed","Negative","Positive"))
  ) %>%
  group_by(idx,stim,nt,nt_metric,trial,block_type,trial_within_block,outcome) %>%
  summarise(mTrial = mean(nt_val),
            rew = ifelse(outcome == 0,
                         ifelse(block_type == "Negative" | block_type == "Mixed",-10,outcome
                         ),
                         ifelse(outcome == 1, 
                                ifelse(block_type == "Positive" | block_type == "Mixed",10, 0),
                                outcome
                         )
            ),
            rew_f = ifelse(outcome == 0, "lose", "win"),
  ) %>% ungroup()



fig.ug.trial = 
  ug.trial %>%
  filter(nt_metric == "Oz",
         nt == "NE") %>%
  group_by(idx,stim,nt,trial) %>%
  summarise(m = mean(mTrial)) %>%
  mutate(nt = factor(nt, levels = c("DA","SE","NE"), labels = c("DA","5-HT","NE")),
         stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         color = ifelse(stim == "Pre-Stim", "#c7e9c0","#006d2c")
  ) %>%
  ggplot(   aes(x = trial, y = m, fill = color, group = stim)) +
  theme_pubr(base_size = 14) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  stat_summary( aes(color = color),geom = "line", size = 1, linewidth = 1.5) +
  stat_summary(
    geom = "ribbon", alpha = 0.35,show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_identity() +
  # Common styling
  labs(y = "NE Estimate [z]", x = "Trial", title = "Ultimatum game") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        strip.text = element_text(color = "black"),
        panel.spacing = unit(0.2, "cm")) +
  coord_cartesian(ylim = c(-6, 6)) +
  theme(legend.key = element_rect(color = NA)) 

fig.rl.trial = 
  rl.trial %>%
  filter(nt_metric == "Oz",
         nt == "NE") %>%
  group_by(idx,stim,nt,trial) %>%
  summarise(m = mean(mTrial)) %>%
  mutate(nt = factor(nt, levels = c("DA","SE","NE"), labels = c("DA","5-HT","NE")),
         stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         color = ifelse(stim == "Pre-Stim", "#c7e9c0","#006d2c")
  ) %>%
  ggplot(   aes(x = trial, y = m, fill = color, group = stim)) +
  theme_pubr(base_size = 14) +
  scale_color_identity() +
  scale_fill_identity() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  stat_summary( aes(color = color),geom = "line", size = 1, linewidth = 1.5) +
  stat_summary(
    geom = "ribbon", alpha = 0.35,
    show.legend = FALSE
  ) +
  # Common styling
  labs(y = "NE Estimate [z]", x = "Trial", title = "Reversal learning") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        strip.text = element_text(color = "black"),
        panel.spacing = unit(0.2, "cm")) +
  coord_cartesian(ylim = c(-6, 6)) +
 # guides(color = "none",
#         fill = "none",
 # ) +
  theme(legend.key = element_rect(color = NA)) 

fig.ug.task = ug.trial %>%
  filter(nt_metric == "Oz", nt == "NE") %>%
  group_by(idx,stim,nt,offer_bin) %>%
  summarise(m = mean(mTrial)) %>%
  mutate(nt = factor(nt, levels = c("DA","SE","NE"), labels = c("DA","5-HT","NE")),
         stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         color = ifelse(stim == "Pre-Stim", "#c7e9c0","#006d2c")
  ) %>%
  ggplot(   aes(x = offer_bin, y = m, fill = color, group = interaction(nt,stim))) +
  theme_pubr(base_size = 14) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  geom_point(position = position_dodge2(width=.25),aes(color = color),
             alpha = .25, size=2) +
  stat_summary(size=1,linewidth=1.5,aes(color = color),
               position = position_dodge2(width  = 0.5)) +
  scale_color_identity() +
  scale_fill_identity() +
  # Common styling
  labs(y = "NE Estimate [z]", x = "Offer", title = "Ultimatum game") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        strip.text = element_text(color = "black"),
        panel.spacing = unit(0.2, "cm")) +
  # coord_cartesian(ylim = c(-6, 6)) +
  guides(color = "none",
         fill = "none",
  ) +
  coord_cartesian(ylim = c(-6, 6)) +
  theme(legend.key = element_rect(color = NA)) 

fig.rl.task = rl.trial %>%
  filter(nt_metric == "Oz", nt == "NE") %>%
  group_by(idx,stim,nt,block_type) %>%
  summarise(m = mean(mTrial)) %>%
  mutate(nt = factor(nt, levels = c("DA","SE","NE"), labels = c("DA","5-HT","NE")),
         stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         block_type = factor(block_type,
                             levels = c("Positive","Mixed","Negative"),
                             labels = c("Reward","Mixed","Punishment")
         ),
         color = ifelse(stim == "Pre-Stim", "#c7e9c0","#006d2c")
  ) %>%
  ggplot(   aes(x = block_type, y = m, fill = color, group = stim)) +
  theme_pubr(base_size = 14) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  geom_point(position = position_dodge2(width=.25),aes(color = color),
             alpha = .25, size=2) +
  stat_summary(size=1,linewidth=1.5,aes(color = color),
               position = position_dodge2(width  = 0.5)) +
  scale_color_identity() +
  scale_fill_identity() +
  coord_cartesian(ylim = c(-6, 6)) +
    # Common styling
  labs(y = "NE Estimate [z]", x = "Block type", title = "Reversal learning") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA),
        strip.text = element_text(color = "black"),
        panel.spacing = unit(0.2, "cm"),
       # axis.text.x = element_text(angle = 20, vjust = 1, hjust=0.5)
       ) +
  #coord_cartesian(ylim = c(-6, 6)) +
  guides(color = "none",
         fill = "none",
  ) +
  theme(legend.key = element_rect(color = NA)) 

main_plot <- ((fig.rl.trial / fig.ug.trial)  | (fig.rl.task/fig.ug.task)  ) 

create_legend <- function() {
  # Dummy data for Post-Stim legend
  legend_data <- data.frame(
    stim = factor(c("Pre-Stim", "Post-Stim"), levels = c("Pre-Stim", "Post-Stim")),
    color = c("#c7e9c0", "#006d2c"),
    x = 1,
    y = 1
  )
  
  ggplot(legend_data, aes(x = x, y = y, color = stim)) +
    geom_point(size = 0) +  # Invisible points
    geom_line(size = 2) +   # This creates the line legend
    scale_color_manual(
      name = element_blank(),
      values = c("Pre-Stim" = "#c7e9c0", "Post-Stim" = "#006d2c")
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

ne_legend <- get_legend(create_legend())

final_plot <- plot_grid(main_plot, ne_legend, 
                        ncol = 2, rel_widths = c(4, 0.5),
                        align = "h", axis = "tb")

ggsave(fig_dir / "Extended-Data_Figure2.tiff", 
       plot = final_plot,
       device = "tiff",
       width = 12,          # Width in inches
       height = 5,         # Height in inches  
       dpi = 300)           # Resolution (300 DPI for publication quality)
