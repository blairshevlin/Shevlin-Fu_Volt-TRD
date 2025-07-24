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
dir = path(here())
clin_dir = dir / "data" / "clinical" # Updated
beh_dir = dir / "data" / "behavior"
nt_dir = dir / "data" / "nt" / "processed"
res_dir = dir / "results" # Updated

add_dbs_rows <- function(df) {
  # Get unique IDs
  unique_ids <- unique(df$id)
  
  # Create new dataframe with DBS rows
  dbs_rows <- data.frame(
    id = unique_ids,
    sess = "DBS",
    mChoice = NA,
    mRT = NA,
    mLogRT = NA
  )
  
  # Combine original data with new DBS rows
  result_df <- bind_rows(df, dbs_rows) %>%
    # Create factor for session ordering
    mutate(
      sess = factor(
        sess,
        levels = c("Baseline", "Pre-Stim", "DBS", "Post-Stim", "Week 1",
                   "Month 1", "Month 2", "Month 3", "Month 4", "Month 5", "Month 6")
      )
    ) %>%
    # Sort by ID and session
    arrange(id, sess)
  
  return(result_df)
}


# Subject-level NT Data
load(nt_dir / "UG_RL_NT-Continuous_7-14-25.RData")

# Load clinical data
cl.df = read.csv(clin_dir / "clinical-data_deid_07-10-25.csv" ) 


ug.nt.cl = 
  ug.EST.Offer %>%
  pivot_longer(cols = c("Oz","Rz","Pz","Mz","Totz"),
               names_to = "nt_metric", values_to = "nt_val") %>%
  group_by(idx,stim,nt,nt_metric) %>%
  summarise(mTrial = mean(nt_val)) %>%
  filter(nt_metric == "Totz") %>%
  select(idx,stim,nt,mTrial) %>%
  pivot_wider(values_from = mTrial, names_from = c("stim","nt")) %>%
  mutate(DA_Pre_UG = `Pre-Stim_DA`,
         DA_Post_UG = `Post-Stim_DA`,
         deltaDA_UG = DA_Post_UG - DA_Pre_UG,
         SE_Pre_UG = `Pre-Stim_SE`,
         SE_Post_UG = `Post-Stim_SE`,
         deltaSE_UG = SE_Post_UG - SE_Pre_UG,
         NE_Pre_UG = `Pre-Stim_NE`,
         NE_Post_UG = `Post-Stim_DA`,
         deltaNE_UG = NE_Post_UG - NE_Pre_UG,
  ) %>%
  select(idx,
         DA_Pre_UG,DA_Post_UG,deltaDA_UG,
         SE_Pre_UG,SE_Post_UG,deltaSE_UG,
         NE_Pre_UG,NE_Post_UG,deltaNE_UG
  ) 

rl.nt.cl = 
  rl.EST.Reward %>%
  pivot_longer(cols = c("Oz","Rz","Pz","Mz","Totz"),
               names_to = "nt_metric", values_to = "nt_val") %>%
  group_by(idx,stim,nt,nt_metric) %>%
  summarise(mTrial = mean(nt_val)) %>%
  filter(nt_metric == "Totz") %>%
  select(idx,stim,nt,mTrial) %>%
  pivot_wider(values_from = mTrial, names_from = c("stim","nt")) %>%
  mutate(DA_Pre_RL = `Pre-Stim_DA`,
         DA_Post_RL = `Post-Stim_DA`,
         deltaDA_RL = DA_Post_RL - DA_Pre_RL,
         SE_Pre_RL = `Pre-Stim_SE`,
         SE_Post_RL = `Post-Stim_SE`,
         deltaSE_RL = SE_Post_RL - SE_Pre_RL,
         NE_Pre_RL = `Pre-Stim_NE`,
         NE_Post_RL = `Post-Stim_DA`,
         deltaNE_RL = NE_Post_RL - NE_Pre_RL,
  ) %>%
  select(idx,
         DA_Pre_RL,DA_Post_RL,deltaDA_RL,
         SE_Pre_RL,SE_Post_RL,deltaSE_RL,
         NE_Pre_RL,NE_Post_RL,deltaNE_RL
  ) 

nt.cl = merge(ug.nt.cl,rl.nt.cl,all.x = TRUE)

cl.Oz = merge(cl.df,nt.cl)

cl.Oz.lme = cl.Oz %>%
  filter(!session %in% c("post stim")) %>%
  mutate(sess_f = factor(session, levels  = c("fmri","pre stim", "week 1", "month 1", "month 2", "month 3", "month 4", "month 5", "month 6")) ) %>%
  group_by(idx) %>%
  mutate(baseline_HDRS = HDRS[session == "fmri"],#HDRS[session == "pre stim"],
         baseline_MADRS = MADRS[session == "fmri"],#MADRS[session == "pre stim"],
         m6_HDRS = HDRS[session == "month 6"],
         m6_MADRS = MADRS[session == "month 6"],
         deltaPerDA_UG = (DA_Post_UG - DA_Pre_UG)/DA_Pre_UG,
         deltaPerSE_UG = (SE_Post_UG - SE_Pre_UG)/SE_Pre_UG,
         deltaPerNE_UG = (NE_Post_UG - NE_Pre_UG)/NE_Pre_UG,
         deltaPerDA_RL = (DA_Post_RL - DA_Pre_RL)/DA_Pre_RL,
         deltaPerSE_RL = (SE_Post_RL - SE_Pre_RL)/SE_Pre_RL,
         deltaPerNE_RL = (NE_Post_RL - NE_Pre_RL)/NE_Pre_RL,
         deltaPerHDRS = (HDRS[session == "fmri"] - HDRS[session == "month 6"])/HDRS[session == "fmri"]
  ) %>%
  filter(!session %in% c("pre stim","fmri")) %>%
  ungroup()

data = cl.Oz.lme %>%
  filter(session == "month 6") %>%
  select(idx,HDRS,baseline_HDRS,deltaPerHDRS,
         deltaDA_UG,deltaSE_UG,deltaNE_UG,
         deltaDA_RL,deltaSE_RL,deltaNE_RL) 

# Create synergy score and additional metrics
data_enhanced <- data %>%
  mutate(
    # Synergy when both change in beneficial direction (DA+, SE+)
    synergy_score = case_when(
      deltaDA_UG > 0 & deltaSE_UG > 0 ~ deltaDA_UG * deltaSE_UG,
      deltaDA_UG < 0 & deltaSE_UG < 0 ~ abs(deltaDA_UG * deltaSE_UG),
      TRUE ~ -abs(deltaDA_UG * deltaSE_UG)  # Opposing changes
    ),
    
    # Alternative synergy: normalized product
    synergy_normalized = deltaDA_UG * deltaSE_UG,
    
    # Change magnitude
    change_magnitude =  HDRS - baseline_HDRS,
    #sqrt(deltaDA_UG^2 + deltaSE_UG^2),
    
    # Response categories
    response_category = case_when(
      HDRS <= 8 ~ "Remission",
      deltaPerHDRS > 0.5 ~ "Responder", 
      TRUE ~ "Non-Responder"
    ),
    response_category = factor(response_category,
                               levels = c("Non-Responder",
                                          "Responder",
                                          "Remission")),
    
    # Change pattern categories
    change_pattern = case_when(
      deltaDA_UG > 0 & deltaSE_UG > 0 ~ "Both Increase",
      deltaDA_UG < 0 & deltaSE_UG < 0 ~ "Both Decrease",
      deltaDA_UG > 0 & deltaSE_UG < 0 ~ "DA↑/5-HT↓",
      TRUE ~ "DA↓/5-HT↑"
    )
  )

# Calculate effect sizes
model_full <- lm(HDRS ~ deltaDA_UG * deltaSE_UG, data = data)
effect_sizes <- eta_squared(model_full, partial = FALSE)
cohens_f_overall <- cohens_f(model_full)

# =============================================================================

fig.synergy <- 
  ggplot(data_enhanced, aes(x = synergy_score, y = HDRS)) +
  # Main relationship
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "black", alpha = 0.2) +
  
  # Points colored by change pattern
  geom_point(aes(fill = change_pattern, size = change_magnitude), 
             shape = 21, color = "black", stroke = 1.5, alpha = 0.8) +
  
  # Add correlation info
  geom_text(x = max(data_enhanced$synergy_score), y = max(data_enhanced$HDRS), 
            label = paste0("Pearson's r = ", round(cor(data_enhanced$synergy_score, data_enhanced$HDRS), 3),
                           "\np = ", round(cor.test(data_enhanced$synergy_score, data_enhanced$HDRS)$p.value, 3)),
            hjust = 1.1, vjust = 1.1, size = 5, color = "grey40") +
  
  geom_hline(yintercept = 8,
             linetype = "dashed",
             linewidth = 0.75) +
  # Add remission dat
  geom_text(x = min(data_enhanced$synergy_score),
            y = 9, label = "Clinical Remission",
            size = 5, hjust = 0
  ) +
  
  scale_fill_manual(values = c("Both Increase" = "#2166ac", "Both Decrease" = "#762a83",
                               "DA↑/5-HT↓" = "#d73027", "DA↓/5-HT↑" = "#f4a582"),
                    name = "Change Pattern") +
  scale_size_continuous(range = c(12, 3), name = "Δ HDRS-17",
                        breaks = c(-10, -15, -20),
                        guide = guide_legend(title.position = "top", 
                                             title.hjust = 0.5,
                                             override.aes = list(fill = "black", alpha = 0.5))
  ) +
  
  labs(
    x = "ΔDA × ΔSE",#"Joint Neurochemical Change\nΔDA × ΔSE",
    y = "Month 6 HDRS-17",
    title = "Neurochemical Profiles Predict Treatment Response",
    #subtitle = "Positive synergy (both DA and SE increase) associated with lower depression"
  ) +
  
  theme_pubr(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    legend.position = "right"
  )


fig.profiles <- 
  ggplot(data_enhanced, aes(x = deltaDA_UG, y = deltaSE_UG)) +
  # Quadrant backgrounds
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, 
           alpha = 0.5, fill = "#2166ac") +  # Both positive
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, 
           alpha = 0.5, fill = "#762a83") +  # Both negative
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, 
           alpha = 0.5, fill = "#d73027") +     # DA+/SE-
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, 
           alpha = 0.5, fill = "#f4a582") +     # DA-/SE+
  
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  
  # Points sized by HDRS and colored by response
  geom_point(aes(size = change_magnitude, fill = change_pattern), #fill = response_category), 
             shape = 21, color = "black", 
             stroke = 1.5, alpha = 0.8) +
  
  # Quadrant labels
  annotate("text", x = 60, y = 50, 
           label = "Both ↑",
           size = 4, fontface = "bold") +
  annotate("text", x = -60, y = -50, label = "Both ↓", 
           size = 4, fontface = "bold") +
  annotate("text", x = 60, y = -50, 
           label = "DA↑/5-HT↓", 
           size = 4, fontface = "bold") +
  annotate("text", x = -60, y = 50, 
           label = "DA↓/5-HT↑",
           size = 4, fontface = "bold") +
  
  scale_fill_manual(values = c("Both Increase" = "#2166ac", "Both Decrease" = "#762a83",
                               "DA↑/5-HT↓" = "#d73027", "DA↓/5-HT↑" = "#f4a582"),
                    name = "Change Pattern") +
  scale_size_continuous(range = c(12, 3), name = "Δ HDRS-17",
                        breaks = c(-10, -15, -20),
                        guide = guide_legend(title.position = "top", 
                                             title.hjust = 0.5,
                                             override.aes = list(fill = "black", alpha = 0.5))
  ) +
  labs(
    x = "Change in DA (ΔDA)",
    y = "Change in 5-HT (Δ5-HT)", 
    title = "Individual Patient Neurochemical Profiles",
  ) +
  
  theme_pubr(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    legend.position = "right"
  )

# Create offset data with boxplot positions opposite to point positions

cl.O.fig = cl.Oz %>%
  filter(!session %in% c('fmri')) %>%
  mutate(idx = factor(idx),
         sess_fig = factor(session,levels  = c("pre stim", "post stim","week 1", "month 1", "month 2", "month 3", "month 4", "month 5", "month 6"),
                           labels = c("Baseline","DBS", "Week 1", "Month 1", "Month 2", "Month 3", "Month 4", "Month 5", "Month 6")) ) %>%
  group_by(idx) %>%
  mutate(
    HDRS_CP = (HDRS - HDRS[session == "pre stim"] )/HDRS[session == "pre stim"],
    remission = ifelse(HDRS < 8, 1, 0),
    responder = factor(ifelse(abs(HDRS_CP) > .5,1,0),
                       levels = c(1,0),
                       labels = c("Responder","Nonresponder")),
    M6_responder = responder[sess_fig == "Month 6"],
    M6_remission = remission[sess_fig == "Month 6"],
    cohort = ifelse(idx %in% c(1:5),"RC+S","PC")
  )%>%
  ungroup() %>% as.data.frame()

offset_data <- cl.O.fig %>%
  filter(sess_fig %in% c("Baseline", "Month 6")) %>%
  mutate(
    # For points: Baseline to right, Month 6 to left
    point_x_pos = ifelse(sess_fig == "Baseline", 1 + 0.1, 2 - 0.1),
    # For boxplots: Baseline to left, Month 6 to right (opposite of points)
    box_x_pos = ifelse(sess_fig == "Baseline", 1 - 0.1, 2 + 0.1),
    # Convert sess_fig to numeric for proper plotting
    sess_num = ifelse(sess_fig == "Baseline", 1, 2)
  )

# Plot with everything positioned properly
fig.hdrs.change =
  ggplot() +
  theme_pubr(base_size = 14) +
  geom_hline(yintercept = 8,
             linetype = "dashed",
             linewidth = 0.75) +
  # Add remission data
  annotate("text",
           x = 0.5,
           y = 9.5, 
           label = "Clinical Remission",
           size = 5, hjust = 0
  ) +
  # Boxplots with offset to opposite side of points
  geom_boxplot(
    data = offset_data,
    aes(x = box_x_pos, y = HDRS, group = sess_fig),
    linewidth = 1.1, 
    outlier.alpha = 0,
    width = 0.15,  # Slightly narrower boxplots
    show.legend = F
  ) +
  # Lines using point positions
  geom_line(
    data = offset_data,
    aes(x = point_x_pos, y = HDRS, group = idx),
    linewidth = 1.5, 
    alpha = .25, 
    color = "black"
  ) +
  # Points
  geom_point(
    data = offset_data,
    aes(x = point_x_pos, y = HDRS, group = idx),
    size = 3, 
    stroke = 1.5, 
    alpha = .35, 
    show.legend = F
  ) +
  scale_color_manual(values = id_colors$color) +
  coord_cartesian(ylim=c(0,32))+
  # Custom x-axis labels at the original positions
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Baseline", "Month 6"),
    limits = c(0.5, 2.5)
  ) +
  labs(x = element_blank(), y = "HDRS-17",
       title = "Month 6 Symptom Change")


fig.4 = (fig.hdrs.change + fig.profiles + fig.synergy) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 16, face = "bold"),
                                plot.subtitle = element_text(size = 12, color = "grey40"))
  ) +
  plot_layout(guides = "collect")

ggsave(res_dir / "figure4.png", 
       plot = fig.4,
       device = "png",
       width = 16,          # Width in inches
       height = 5,         # Height in inches  
       dpi = 300)           # Resolution (300 DPI for publication quality)
