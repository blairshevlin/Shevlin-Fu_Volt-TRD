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


# Subject-level Data
nt.df = read.csv(file = nt_dir / "NT_Oz_Means_06-27-2024.csv") %>% select(!X)

# Load clinical data
cl.df = read.csv(clin_dir / "clinical-data_deid_07-10-25.csv" ) 
cl.Oz = merge(cl.df,nt.df)

cl.Oz.lme = cl.Oz %>%
  filter(!session %in% c("post stim")) %>%
  mutate(sess_f = factor(session, levels  = c("fmri","pre stim", "week 1", "month 1", "month 2", "month 3", "month 4", "month 5", "month 6")) ) %>%
  group_by(id) %>%
  mutate(baseline_HDRS = HDRS[session == "fmri"],#HDRS[session == "pre stim"],
         baseline_MADRS = MADRS[session == "fmri"],#MADRS[session == "pre stim"],
         m6_HDRS = HDRS[session == "month 6"],
         m6_MADRS = MADRS[session == "month 6"],
         deltaPerDA_UG = (DA_Post_UG - DA_Pre_UG)/DA_Pre_UG,
         deltaPerSE_UG = (SE_Post_UG - SE_Pre_UG)/SE_Pre_UG,
         deltaPerDA_RL = (DA_Post_RL - DA_Pre_RL)/DA_Pre_RL,
         deltaPerSE_RL = (SE_Post_RL - SE_Pre_RL)/SE_Pre_RL,
         deltaPerHDRS = (HDRS[session == "fmri"] - HDRS[session == "month 6"])/HDRS[session == "fmri"]
  ) %>%
  filter(!session %in% c("pre stim","fmri")) %>%
  ungroup()

cl.change = cl.Oz.lme %>% 
  select(id,sess_f, baseline_HDRS, HDRS, deltaSE_UG, deltaSE_RL,  deltaDA_UG, deltaDA_RL ) %>%
  distinct() %>%
  mutate(HDRS_C = baseline_HDRS - HDRS,
         month = recode(sess_f,
                        "week 1" = "W1",
                        "month 1" =  "M1",
                        "month 2" = "M2",
                        "month 3" = "M3",
                        "month 4" = "M4",
                        "month 5" = "M5",
                        "month 6" = "M6"),
         HDRS_CP = (HDRS_C)/baseline_HDRS) %>%
  select(!sess_f)

# Ultimatum Game behavior
ug.beh = read.csv(file = beh_dir / "UG_Beh_AllSess_06_26_2024.csv") %>% 
  mutate(sess = factor(session, 
                       levels = c("fmri","pre_stim","post_stim","1w","m1","m2","m3","m4","m5","m6"),
                       labels = c("fMRI", "Pre-Stim", "Post-Stim", "Week 1", "Month 1", 
                                  "Month 2","Month 3", "Month 4", "Month 5", "Month 6"))
  ) %>% as.data.frame()

ug.beh.means = ug.beh %>% group_by(id,sess) %>%
  summarise(mChoice = mean(rej==0),
            mRT = mean(rt),
            mLogRT = mean(log(rt))) %>%   
  mutate(id = factor(id),
         sess = recode(sess,"fMRI" = "Baseline")) %>%
  as.data.frame()

ug.beh.means = add_dbs_rows(ug.beh.means)

ug.nt = nt.df  %>%
  filter(!is.na(DA_Pre_UG)) %>%
  mutate(DA = DA_Post_UG - DA_Pre_UG,
         SE = SE_Post_UG - SE_Pre_UG) %>%
  select(c(id, DA,SE))

ug.beh.change = 
  ug.beh.means %>%
  mutate(mRT = mLogRT) %>%
  select(!mLogRT) %>%
  pivot_wider(names_from = sess, values_from = c(mChoice, mRT)) %>%
  mutate(# Use fMRI as baseline
    mChoice_PostStim_Change = `mChoice_Post-Stim` - mChoice_Baseline,
    mChoice_W1_Change = `mChoice_Week 1` - mChoice_Baseline,
    mChoice_M1_Change = `mChoice_Month 1` - mChoice_Baseline,
    mChoice_M2_Change = `mChoice_Month 2` - mChoice_Baseline,
    mChoice_M3_Change = `mChoice_Month 3` - mChoice_Baseline,
    mChoice_M4_Change = `mChoice_Month 4` - mChoice_Baseline,
    mChoice_M5_Change = `mChoice_Month 5` - mChoice_Baseline,
    mChoice_M6_Change = `mChoice_Month 6` - mChoice_Baseline,
    mRT_PostStim_Change = `mRT_Post-Stim` - mRT_Baseline,
    mRT_W1_Change = `mRT_Week 1` - mRT_Baseline,
    mRT_M1_Change = `mRT_Month 1` - mRT_Baseline,
    mRT_M2_Change = `mRT_Month 2` - mRT_Baseline,
    mRT_M3_Change = `mRT_Month 3` - mRT_Baseline,
    mRT_M4_Change = `mRT_Month 4` - mRT_Baseline,
    mRT_M5_Change = `mRT_Month 5` - mRT_Baseline,
    mRT_M6_Change = `mRT_Month 6` - mRT_Baseline
  ) %>%
  select(c(id,mChoice_PostStim_Change,mChoice_W1_Change,mChoice_M1_Change,mChoice_M2_Change,
           mChoice_M3_Change,mChoice_M4_Change,mChoice_M5_Change,mChoice_M6_Change,
           mRT_PostStim_Change,mRT_W1_Change,mRT_M1_Change,mRT_M2_Change,
           mRT_M3_Change,mRT_M4_Change,mRT_M5_Change,mRT_M6_Change ))

ug.beh.nt = merge(ug.beh.change, ug.nt) 

ug.beh.nt.long = ug.beh.nt %>% 
  pivot_longer(cols = c(DA,SE), names_to = "NT", values_to = "ests") %>%
  pivot_longer(cols = c(mChoice_PostStim_Change,mChoice_W1_Change,mChoice_M1_Change,mChoice_M2_Change,
                        mChoice_M3_Change,mChoice_M4_Change,mChoice_M5_Change,mChoice_M6_Change,
                        mRT_PostStim_Change,mRT_W1_Change,mRT_M1_Change,mRT_M2_Change,
                        mRT_M3_Change,mRT_M4_Change,mRT_M5_Change,mRT_M6_Change), 
               names_to = "change") %>%
  mutate(month = str_extract(change, "(?<=_)[^_]+(?=_)"),
         beh = ifelse(grepl("mChoice", change), "choice", "rt"),
         baseline = ifelse(grepl("Change2",change),"w1","fmri"))

ug.beh.change.w =
  ug.beh.nt.long %>%
  select(!c(change,NT,ests,baseline)) %>%
  filter(month!="PostStim") 

ug.beh.change.cl = 
  merge(cl.change, ug.beh.change.w, by = c("id","month")) %>%
  distinct() %>%
  pivot_wider(names_from=beh, values_from = value)

ug.beh.change.cl %>%
  select(id,deltaDA_UG,deltaSE_UG) %>% 
  distinct() %>%
  pivot_longer(cols = c(deltaDA_UG,deltaSE_UG)) %>%
  group_by(name) %>%
  summarise(m = mean(value),
            sd = sd(value),
            hsf = sd(value)/2,
            `m+0.5` = m + hsf,
            `m-0.5` = m - hsf,
            `m+1` = m + sd,
            `m-1` = m -sd)

data = cl.Oz.lme %>%
  filter(session == "month 6") %>%
  select(id,HDRS,baseline_HDRS,deltaPerHDRS,deltaDA_UG,deltaSE_UG) 

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
# FIGURE 1: Synergy Score vs Depression Outcome
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
  scale_size_continuous(range = c(12, 4), name = "Δ HDRS-17",
                        breaks = c(-5, -10, -15, -20),
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

# =============================================================================
# FIGURE 2: Effect Size Visualization
# =============================================================================

# Prepare effect size data
effect_data <- data.frame(
  term = c("ΔDA", "ΔSE", "ΔDA × ΔSE"),
  eta_squared = c(effect_sizes$Eta2[1], effect_sizes$Eta2[2], effect_sizes$Eta2[3]),
  ci_low = c(effect_sizes$CI_low[1], effect_sizes$CI_low[2], effect_sizes$CI_low[3]),
  ci_high = c(effect_sizes$CI_high[1], effect_sizes$CI_high[2], effect_sizes$CI_high[3]),
  interpretation = case_when(
    c(effect_sizes$Eta2[1], effect_sizes$Eta2[2], effect_sizes$Eta2[3]) < 0.01 ~ "Negligible",
    c(effect_sizes$Eta2[1], effect_sizes$Eta2[2], effect_sizes$Eta2[3]) < 0.06 ~ "Small",
    c(effect_sizes$Eta2[1], effect_sizes$Eta2[2], effect_sizes$Eta2[3]) < 0.14 ~ "Medium",
    TRUE ~ "Large"
  )
) %>%
  mutate(
    term = factor(term, levels = c("ΔDA × ΔSE", "ΔDA", "ΔSE")),  # Reorder for importance
    interpretation = factor(interpretation, levels = c("Negligible", "Small", "Medium", "Large"))
  )

fig.effects <- 
  ggplot(effect_data, aes(x = eta_squared, y = term)) +
  # Effect size bars
  geom_col(alpha = 0.8, width = 0.6) +
  
  # Confidence intervals
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2, alpha = 0.7) +
  
  # Add effect size values
  geom_text(aes(label = round(eta_squared, 3)), 
            hjust = -0.1, vjust = -0.1,
            size = 3.5, fontface = "bold") +
  
  # Reference lines for effect size interpretation
  geom_vline(xintercept = c(0.01, 0.06, 0.14), linetype = "dashed", alpha = 0.5) +
  annotate("text", x = c(-0.01, 0.035, 0.11), y = 0.5, 
           label = c("Small", "Medium", "Large"), angle = 90, 
           size = 3.5, hjust = 0) +
  
  labs(
    x = "η² (Proportion of Variance Explained)",
    y = "",
    title = "Effect Sizes: Interaction Dominates Individual Effects", 
    subtitle = paste0("Overall model Cohen's f = ", round(cohens_f_overall$Cohens_f, 3), " (large effect)")
  ) +
  
  theme_pubr(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    axis.text.y = element_text(size = 11, face = "bold")
  )

# =============================================================================
# FIGURE 3: Change Pattern Breakdown
# =============================================================================

pattern_summary <- data_enhanced %>%
  group_by(change_pattern) %>%
  summarise(
    n = n(),
    mean_HDRS = mean(HDRS),
    sd_HDRS = sd(HDRS),
    mean_synergy = mean(synergy_score),
    .groups = 'drop'
  ) %>%
  mutate(
    se_HDRS = sd_HDRS / sqrt(n),
    change_pattern = factor(change_pattern, 
                            levels = c("Both Decrease", "DA-/SE+", "DA+/SE-", "Both Increase"))
  )

fig.patterns <- 
  ggplot(pattern_summary, aes(x = change_pattern, y = mean_HDRS, fill = change_pattern)) +
  # Bars with error bars
  geom_col(alpha = 0.8, width = 0.7) +
  geom_errorbar(aes(ymin = mean_HDRS - se_HDRS, ymax = mean_HDRS + se_HDRS),
                width = 0.2, alpha = 0.7) +
  
  # Add sample sizes
  geom_text(aes(label = paste0("n=", n)), y = 1, color = "white", fontface = "bold") +
  
  # Add mean values
  geom_text(aes(label = round(mean_HDRS, 1)), 
            y = pattern_summary$mean_HDRS + pattern_summary$se_HDRS + 1,
            fontface = "bold", size = 4) +
  
  scale_fill_manual(values = c("Both Increase" = "#2166ac", "DA+/SE-" = "#d73027",
                               "DA-/SE+" = "#f4a582", "Both Decrease" = "#762a83")) +
  
  labs(
    x = "Neurochemical Change Pattern",
    y = "Month 6 HDRS-17",
    title = "Treatment Response by Change Pattern",
    subtitle = "Both DA and SE increasing shows best outcomes"
  ) +
  
  theme_pubr(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    legend.position = "none"#,
   # axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  coord_cartesian(ylim=c(0,25))

# =============================================================================
# FIGURE 4: Individual Patient Profiles
# =============================================================================

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
  geom_point(aes(size = change_magnitude), #fill = response_category), 
             shape = 21, color = "black", 
             fill = "black",
             stroke = 1.5, alpha = 0.5) +
  
  # Quadrant labels
  annotate("text", x = 5, y = 4.5, 
           #label = "Joint\n(Both ↑)", 
           label = "Both ↑",
           size = 5, fontface = "bold") +
  annotate("text", x = -5, y = -4.5, label = "Both ↓", 
           size = 5, fontface = "bold") +
  annotate("text", x = 5, y = -4.5, 
           #label = "Asymmetric\n(DA↑/5-HT↓)", 
           label = "DA↑/5-HT↓", 
           size = 5, fontface = "bold") +
  annotate("text", x = -5, y = 4.5, 
           #label = "Asymmetric\n(DA↓/5-HT↑)", 
           label = "DA↓/5-HT↑",
           size = 5, fontface = "bold") +
  
  scale_size_continuous(range = c(12, 4), name = "Δ HDRS-17",
                        breaks = c(-10,-15,-20),
                        guide = guide_legend(title.position = "top", 
                                             title.hjust = 0.5)) +
 # scale_fill_manual(values = c("Remission" = "#2166ac", 
  #                             "Responder" = "darkgreen",
   #                            "Non-Responder" = "#d73027"),
    #                name = "Clinical Response") +
  
  labs(
    x = "Change in DA (ΔDA)",
    y = "Change in 5-HT (Δ5-HT)", 
    title = "Individual Patient Neurochemical Profiles",
   # subtitle = "Green quadrant shows synergistic changes with best outcomes"
  ) +
  
  theme_pubr(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    legend.position = "right"
  )

# =============================================================================
# COMBINED LAYOUT
# =============================================================================

# Combine plots
combined_plot <- (fig.patterns + fig.profiles) / (fig.synergy + fig.effects) +
  plot_annotation(
    title = "Neurochemical Synergy Predicts SCC DBS Response",
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 12, color = "grey40"))
  )

ggsave("neurochemical_synergy_analysis.tiff", 
       plot = combined_plot,
       device = "tiff",
       width = 16,          # Width in inches
       height = 12,         # Height in inches  
       dpi = 300,           # Resolution (300 DPI for publication quality)
       compression = "lzw") # Compression method

# Display individual plots
print("Figure 1: Synergy Score vs Depression")
print(p1_synergy)

print("\nFigure 2: Effect Sizes")  
print(p2_effects)

print("\nFigure 3: Change Patterns")
print(p3_patterns)

print("\nFigure 4: Individual Profiles")
print(p4_profiles)

print("\nCombined Figure:")
print(combined_plot)

# Print summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat("Synergy Score Correlation with HDRS:\n")
cor_result <- cor.test(data_enhanced$synergy_score, data_enhanced$HDRS)
cat("r =", round(cor_result$estimate, 3), ", p =", round(cor_result$p.value, 3), "\n\n")

cat("Effect Sizes (η²):\n")
print(effect_sizes)

cat("\nChange Pattern Summary:\n")
print(pattern_summary)



offset_data <-  cl.Oz %>%
  filter(!session %in% c('fmri')) %>%
  mutate(id = factor(id),
         sess_fig = factor(session,levels  = c("pre stim", "post stim","week 1", "month 1", "month 2", "month 3", "month 4", "month 5", "month 6"),
                           labels = c("Baseline","DBS", "Week 1", "Month 1", "Month 2", "Month 3", "Month 4", "Month 5", "Month 6")) ) %>%
  group_by(id) %>%
  mutate(
    HDRS_CP = (HDRS - HDRS[session == "pre stim"] )/HDRS[session == "pre stim"],
    remission = ifelse(HDRS < 8, 1, 0),
    responder = factor(ifelse(abs(HDRS_CP) > .5,1,0),
                       levels = c(1,0),
                       labels = c("Responder","Nonresponder")),
    M6_responder = responder[sess_fig == "Month 6"],
    M6_remission = remission[sess_fig == "Month 6"],
    cohort = ifelse(id %in% c(809,811,812,814,815),"RC+S","PC")
  )%>%
  ungroup() %>% as.data.frame() %>%
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
            y = 9, label = "Clinical Remission",
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
    aes(x = point_x_pos, y = HDRS, group = id),
    linewidth = 1.5, 
    alpha = .25, 
    color = "black"
  ) +
  # Points
  geom_point(
    data = offset_data,
    aes(x = point_x_pos, y = HDRS, group = id),
    size = 3, 
    stroke = 1.5, 
    alpha = .35, 
    show.legend = F
  ) +
  coord_cartesian(ylim=c(0,32))+
  # Custom x-axis labels at the original positions
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Baseline", "Month 6"),
    limits = c(0.5, 2.5)
  ) +
  labs(x = element_blank(), y = "HDRS-17",
       title = "Month 6 Symptom Change") +
  theme(    plot.title = element_text(size = 14, face = "bold")
            )


combined_three = (fig.hdrs.change + fig.profiles + fig.synergy) +
  plot_annotation(tag_levels = "a",
    #title = "Neurochemical Patterns Predicts SCC DBS Response",
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 12, color = "grey40"))
  )

ggsave("neurochemical_synergy_analysis_3plot.tiff", 
       plot = combined_three,
       device = "tiff",
       width = 18,          # Width in inches
       height = 6,         # Height in inches  
       dpi = 300)           # Resolution (300 DPI for publication quality)
      
