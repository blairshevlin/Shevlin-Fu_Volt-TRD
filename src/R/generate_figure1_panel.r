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
dir = path(here())
nt_dir = dir / "data" / "nt" / "processed"
res_dir = dir / "results" # Updated

# Subject-level NT Data
load(nt_dir / "UG_RL_NT-Raw_10-08-25.RData")

# Color scheme
NT_colors = data.frame(id = c("DA","SE"),
                       color = c("#cb181d","#2171b5"))

# Generate summary data - focusing only on reward events (RL) and offer events (UG)
rl.EST.Reward = rl_proc %>%
  filter(event == "Reward") 

sm.rl = rl.EST.Reward %>%
  filter(!is.na(nM)) %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         sample = ifelse(sample == 51,0,(sample-51)*100),
         sm = runmean(nM,3,align = "right")) %>%
  group_by(idx,stim,nt,trial) %>%
  mutate(#zm = scale(nM)[,1],
    szm = runmean(nMz,3,align = "right")) %>%
  filter(sample %in%  c(-400:800))


fig.data  =   sm.rl %>%
  filter(event == "Reward",
         stim == "Post-Stim",
         nt == "SE",
         trial ==  28,
         idx == 5,
         !is.na(prev_rew)) 


window_sum <- fig.data %>%
  filter(sample >= 0 & sample <= 500) %>%
  pull(szm) %>%  # assuming szm is your z-scored variable
  sum()

# Panel 1: Raw nM estimates
p1 <- ggplot(fig.data, aes(x = sample, y = sm)) +  # using sm without -200 adjustment
  theme_pubr(base_size = 14) +
#  geom_hline(yintercept = mean(fig.data$sm), linewidth = 1.1, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linewidth = 1.1, linetype = "dashed", color = "black") +
  geom_line(linewidth = 1.5, color = NT_colors$color[NT_colors$id == "SE"]) +
  labs(x = "time relative to outcome reveal [ms]",
       y = "5-HT [nM]",
       title = "raw model estimate") +
  coord_cartesian(xlim = c(-250, 750)) +
  annotate("text", x = -40, y = max(fig.data$sm), label = "Outcome\nReveal", 
           size = 3.5, hjust = 1, vjust = 1.2)

nM_range <- range(fig.data$sm)
z_range <- range(fig.data$szm)
scale_factor <- diff(nM_range) / diff(z_range)
intercept_adjust <- nM_range[1] - z_range[1] * scale_factor

p2 <- ggplot(fig.data, aes(x = sample)) +
  theme_pubr(base_size = 14) +
  geom_hline(yintercept = 0 * scale_factor + intercept_adjust, linewidth = 1.1, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linewidth = 1.1, linetype = "dashed", color = "black") +
  geom_line(aes(y = szm * scale_factor + intercept_adjust), linewidth = 1.5, color = NT_colors$color[NT_colors$id == "SE"]) +
  labs(x = "time relative to outcome reveal [ms]",
       y = "5-HT [nM]",
       title = "z-scored estimate") +
  coord_cartesian(xlim = c(-250, 750)) +
  # Add second y-axis for z-scores
  scale_y_continuous(
    name = "5-HT [nM]",
    sec.axis = sec_axis(~ (. - intercept_adjust) / scale_factor, name = "5-HT [z]")
  ) +
  theme(axis.title.y.right = element_text(color = NT_colors$color[NT_colors$id == "SE"], angle = 90),
        axis.text.y.right = element_text(color = NT_colors$color[NT_colors$id == "SE"]))

p3 <- ggplot(fig.data, aes(x = sample)) +
  theme_pubr(base_size = 14) +
  geom_hline(yintercept = 0 * scale_factor + intercept_adjust, linewidth = 1.1, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linewidth = 1.1, linetype = "dashed", color = "black") +
  geom_line(aes(y = szm * scale_factor + intercept_adjust), linewidth = 1.5, color = NT_colors$color[NT_colors$id == "SE"]) +
  geom_ribbon(data = subset(fig.data, sample >= 0 & sample <= 500),
              aes(ymax = szm * scale_factor + intercept_adjust, 
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
        axis.text.y.right = element_text(color = NT_colors$color[NT_colors$id == "SE"]))

# Create arrows between panels


# Function to create arrow grob
arrow_grob <- function() {
  arrow(angle = 20, length = unit(0.3, "cm"), type = "closed")
}

# Combine panels with arrows
combined_plot <- grid.arrange(
  p1, 
  rectGrob(gp = gpar(fill = "transparent", col = "transparent")), # spacer for arrow
  p2,
  rectGrob(gp = gpar(fill = "transparent", col = "transparent")), # spacer for arrow  
  p3,
  ncol = 5,
  widths = c(1, 0.1, 1, 0.1, 1)
)
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
ggsave(res_dir / "Figure1_PanelC.png", 
       plot = final_figure,
       device = "png",
       width = 12,          # Width in inches
       height = 4,         # Height in inches  
       dpi = 300)           # Resolution (300 DPI for publication quality)

ggsave(res_dir / "Figure1_PanelC.pdf", 
       plot = final_figure,
       device = "pdf",
       width = 12,          # Width in inches
       height = 4,         # Height in inches  
       dpi = 300)           # Resolution (300 DPI for publication quality)
