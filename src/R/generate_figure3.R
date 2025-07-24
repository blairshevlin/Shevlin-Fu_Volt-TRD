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
dir = path(here())
nt_dir = dir / "data" / "nt" / "processed"
beh_dir = dir / "data" / "behavior"
res_dir = dir / "results" # Updated

# Subject-level NT Data
load(nt_dir / "UG_RL_NT-Continuous_7-14-25.RData")

# IDs
ids_final = unique(ug.EST.Offer$idx)

# Specify colors
NT_colors = data.frame(id = c("DA","5-HT"),
                       color = c("#cb181d","#2171b5"))

# Load behavioral data
rl.beh = read.csv(file = beh_dir / "rl-data_deid_07-10-25.csv")

ug.beh = read.csv(file = beh_dir / "ug-data_deid_07-10-25.csv")

# Get session-level averages
ug.beh.means = ug.beh %>% group_by(idx,sess) %>%
  summarise(mChoice = mean(rej==0),
            mRT = mean(rt),
            mLogRT = mean(log(rt))) %>%   
  mutate(idx = factor(idx),
         sess = recode(sess,"fMRI" = "Baseline")) %>%
  as.data.frame()

ug.mood.means = ug.beh %>% group_by(idx,sess) %>%
  filter(!is.na(rt_mood)) %>%
  summarise(mMood = mean(mood),
            mRT = mean(rt_mood ),
            mLogRT = mean(log(rt))) %>% 
  mutate(idx = factor(idx),
         sess = recode(sess,"fMRI" = "Baseline")) %>%
  as.data.frame()

rl.beh.means = rl.beh %>%
  group_by(idx,sess) %>%
  summarise(mChoice = mean(opt),
            mRT = mean(rt),
            mLogRT = mean(log(rt))) %>% 
  mutate(idx = factor(idx),
         sess = recode(sess,"fMRI" = "Baseline")) %>%
  as.data.frame()

# Add DBS line
add_dbs_rows <- function(df) {
  # Get unique IDs
  unique_ids <- unique(df$idx)
  
  # Create new dataframe with DBS rows
  dbs_rows <- data.frame(
    idx = unique_ids,
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
    arrange(idx, sess)
  
  return(result_df)
}

add_dbs_rows_mood <- function(df) {
  # Get unique IDs
  unique_ids <- unique(df$idx)
  
  # Create new dataframe with DBS rows
  dbs_rows <- data.frame(
    idx = unique_ids,
    sess = "DBS",
    mMood = NA,
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
    arrange(idx, sess)
  
  return(result_df)
}

rl.beh.means = add_dbs_rows(rl.beh.means)
ug.beh.means = add_dbs_rows(ug.beh.means)
ug.mood.means = add_dbs_rows_mood(ug.mood.means)

# Function to calculate significance and position stars

calculate_significance <- function(data, value_col, time_col, baseline_level, 
                                   comparison_level, subject_col = "idx", paired = TRUE) {
  # Extract baseline and comparison data
  baseline_data <- data[data[[time_col]] == baseline_level, ]
  comparison_data <- data[data[[time_col]] == comparison_level, ]
  
  # Ensure data is properly paired by subject
  merged_data <- merge(
    baseline_data[, c(subject_col, value_col)],
    comparison_data[, c(subject_col, value_col)],
    by = subject_col
  )
  
  # Perform paired t-test
  t_test <- wilcox.test(merged_data[[paste0(value_col, ".x")]], 
                   merged_data[[paste0(value_col, ".y")]], 
                   paired = paired)
  
  # Convert p-value to significance symbols
  sig <- if(t_test$p.value <= 0.001) "***"
  else if(t_test$p.value <= 0.01) "**"
  else if(t_test$p.value <= 0.05) "*"
  else "ns"
  
  return(sig)
}

calculate_star_positions <- function(data, value_col, bracket_height_offset = 3) {
  max_val <- max(data[[value_col]], na.rm = TRUE)
  
  # Position stars just above where brackets would be
  # Adjust bracket_height_offset based on your bracket positions
  positions <- max_val + bracket_height_offset + c(2, 7, 12)
  return(positions)
}

significance_labels_ug_c <- sapply(c("Month 1", "Month 3", "Month 6"), function(comp_level) {
  calculate_significance(ug.beh.means, "mChoice","sess" , "Baseline", 
                         comp_level, "idx", paired = TRUE)
}
)

significance_labels_rl_c <- sapply(c("Month 1", "Month 3", "Month 6"), function(comp_level) {
  calculate_significance(rl.beh.means, "mChoice","sess" , "Baseline", 
                         comp_level, "idx", paired = TRUE)
}
)

significance_labels_ug_logrt <- sapply(c("Month 1", "Month 3", "Month 6"), function(comp_level) {
  calculate_significance(ug.beh.means, "mLogRT","sess" , "Baseline", 
                         comp_level, "idx", paired = TRUE)
}
)

significance_labels_rl_logrt <- sapply(c("Month 1", "Month 3", "Month 6"), function(comp_level) {
  calculate_significance(rl.beh.means, "mLogRT","sess" , "Baseline", 
                         comp_level, "idx", paired = TRUE)
}
)

significance_labels_ug_mood <- sapply(c("Month 1", "Month 3", "Month 6"), function(comp_level) {
  calculate_significance(ug.mood.means, "mMood","sess" , "Baseline", 
                         comp_level, "idx", paired = TRUE)
}
)

rl.c.sess = 
  rl.beh.means %>%
  filter(sess %in% c("Baseline","DBS","Month 1","Month 3","Month 6")) %>%
  ggplot(aes(x = sess, y = mChoice))+
  theme_pubr(base_size = 14) +
  geom_rect(xmin="DBS",xmax="Month 1", color="gray",
            ymin=-Inf, ymax=0.45, 
            fill="gray", alpha=0.15) +
  geom_rect(xmin="Month 1",xmax="Month 6", color="black",
            ymin=-Inf, ymax=0.45, 
            fill="black", alpha=0.15) +
  annotate("text", x = 2.5, y = 0.425, 
           label = "OFF",fontface = "bold",
           color = "white",
           size = 3.5)+ 
  annotate("text", x = "Month 3", y = 0.425, 
           label = "ON",fontface = "bold",
           color = "white",size = 3.5) + 
  geom_vline(xintercept = "DBS", linewidth = 2) +
  coord_cartesian(ylim= c(0.425,1) ) +
  geom_boxplot(data = rl.beh.means[rl.beh.means$sess %in% c("Baseline","DBS","Month 1", "Month 3", "Month 6"),],
               linewidth = 1.1,outlier.alpha = 0,
               show.legend = F) +
  geom_point(data = rl.beh.means[rl.beh.means$sess %in% c("Baseline","Month 1", "Month 3", "Month 6"),],
             size=2,color = "purple",alpha=.5,
             position = position_dodge2(width = .3),
             stroke = 1.75,
             aes( group = idx
             ),
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x =  element_blank(), y = "p(optimal)",
       title = "Choice [RL]",
       color = element_blank(),
       shape = element_blank()) +
  scale_shape_manual(values = c(16,4)) +
  stat_compare_means( 
    label="p.signif",
    comparisons = list(
      c("Baseline","Month 1"),
      c("Baseline","Month 3"),
      c("Baseline","Month 6")),
    method = "t.test",
    bracket.size = 0.5,
    tip.length = 0.02,
    step.increase = 0.08,label.x = c(2.5,2.5,2.5)
  )  

ug.c.sess = 
  ug.beh.means %>%
  filter(sess %in% c("Baseline","DBS","Month 1","Month 3","Month 6")) %>%
  ggplot(aes(x = sess, y = mChoice))+
  theme_pubr(base_size = 14) +
  geom_rect(xmin="DBS",xmax="Month 1", color="gray",
            ymin=-Inf, ymax=0.05, 
            fill="gray", alpha=0.15) +
  geom_rect(xmin="Month 1",xmax="Month 6", color="black",
            ymin=-Inf, ymax=0.05, 
            fill="black", alpha=0.15) +
  annotate("text", x = 2.5, y = 0.005, 
           label = "OFF",fontface = "bold",
           color = "white",
           size = 3.5)+ 
  annotate("text", x = "Month 3", y = 0.005, 
           label = "ON",fontface = "bold",
           color = "white",size = 3.5) + 
  geom_vline(xintercept = "DBS", linewidth = 2) +
  #coord_cartesian(ylim= c(0,1.2) ) +
  geom_boxplot(data = ug.beh.means[ug.beh.means$sess %in% c("Baseline","DBS","Month 1", "Month 3", "Month 6"),],
               linewidth = 1.1,outlier.alpha = 0,
               show.legend = F) +
  geom_point(data = ug.beh.means[ug.beh.means$sess %in% c("Baseline","Month 1", "Month 3", "Month 6"),],
             size=2,color = "purple",alpha=.5,
             position = position_dodge2(width = .3),
             stroke = 1.75,aes( group = idx
             ),
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x =  element_blank(), y = "p(accept)",
       title = "Choice [UG]",
       color = element_blank(),
       shape = element_blank()) +
  scale_shape_manual(values = c(16,4)) +
  stat_compare_means( 
    label="p.signif",
    comparisons = list(
      c("Baseline","Month 1"),
      c("Baseline","Month 3"),
      c("Baseline","Month 6")),
    method = "t.test",
    bracket.size = 0.5,
    tip.length = 0.02,
    step.increase = 0.08,label.x = c(2.5,2.5,2.5)
  )  

ug.logrt.sess = 
  ug.beh.means %>%
  filter(sess %in% c("Baseline","DBS","Month 1","Month 3","Month 6")) %>%
  ggplot(aes(x = sess, y = mLogRT))+
  theme_pubr(base_size = 14) +
  geom_rect(xmin="DBS",xmax="Month 1", color="gray",
            ymin=-Inf, ymax=-0.4, 
            fill="gray", alpha=0.15) +
  geom_rect(xmin="Month 1",xmax="Month 6", color="black",
            ymin=-Inf, ymax= -0.4, 
            fill="black", alpha=0.15) +
  annotate("text", x = 2.5, y = -0.5, 
           label = "OFF",fontface = "bold",
           color = "white",
           size = 3.5)+ 
  annotate("text", x = "Month 3", y = -0.5, 
           label = "ON",fontface = "bold",
           color = "white",size = 3.5) + 
  geom_vline(xintercept = "DBS", linewidth = 2) +
  geom_boxplot(data = ug.beh.means[ug.beh.means$sess %in% c("Baseline","DBS","Month 1", "Month 3", "Month 6"),],
               linewidth = 1.1,outlier.alpha = 0,
               show.legend = F) +
  geom_point(data = ug.beh.means[ug.beh.means$sess %in% c("Baseline","Month 1", "Month 3", "Month 6"),],
             size=2,color = "purple",alpha=.5,
             position = position_dodge2(width = .3),
             stroke = 1.75,aes( group = idx
             ),
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x =  element_blank(), y = "log sec",
       color = element_blank(),title = "Response time [UG]",
       shape = element_blank()) +
  scale_shape_manual(values = c(16,4)) +
  stat_compare_means( 
    label="p.signif",
    comparisons = list(
      c("Baseline","Month 1"),
      c("Baseline","Month 3"),
      c("Baseline","Month 6")),
    method = "t.test",
    bracket.size = 0.5,
    tip.length = 0.02,
    step.increase = 0.08,label.x = c(2.5,2.5,2.5)
  )  


rl.logrt.sess = 
  rl.beh.means %>%
  filter(sess %in% c("Baseline","DBS","Month 1","Month 3","Month 6")) %>%
  ggplot(aes(x = sess, y = mLogRT))+
  theme_pubr(base_size = 14) +
  geom_rect(xmin="DBS",xmax="Month 1", color="gray",
            ymin=-Inf, ymax=-3.75, 
            fill="gray", alpha=0.15) +
  geom_rect(xmin="Month 1",xmax="Month 6", color="black",
            ymin=-Inf, ymax=-3.75, 
            fill="black", alpha=0.15) +
  annotate("text", x = 2.5, y =-4, 
           label = "OFF",fontface = "bold",
           color = "white",
           size = 3.5)+ 
  annotate("text", x = "Month 3", y = -4, 
           label = "ON",fontface = "bold",
           color = "white",size = 3.5) + 
  geom_vline(xintercept = "DBS", linewidth = 2) +
  #coord_cartesian(ylim= c(-0.6,2) ) +
  geom_boxplot(data = rl.beh.means[rl.beh.means$sess %in% c("Baseline","DBS","Month 1", "Month 3", "Month 6"),],
               linewidth = 1.1,outlier.alpha = 0,
               show.legend = F) +
  geom_point(data = rl.beh.means[rl.beh.means$sess %in% c("Baseline","Month 1", "Month 3", "Month 6"),],
             size=2,color = "purple",alpha=.5,
             position = position_dodge2(width = .3),
             stroke = 1.75,aes( group = idx
             ),
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x =  element_blank(), y = "log sec",
       title = "Response time [RL]",
       color = element_blank(),
       shape = element_blank()) +
  scale_shape_manual(values = c(16,4)) +
  stat_compare_means( 
    label="p.signif",
    comparisons = list(
      c("Baseline","Month 1"),
      c("Baseline","Month 3"),
      c("Baseline","Month 6")),
    method = "t.test",
    bracket.size = 0.5,
    tip.length = 0.02,
    step.increase = 0.08,label.x = c(2.5,2.5,2.5)
  )  


ug.mood.sess = 
  ug.mood.means %>%
  filter(sess %in% c("Baseline","DBS","Month 1","Month 3","Month 6")) %>%
  ggplot(aes(x = sess, y = mMood))+
  theme_pubr(base_size = 14) +
  geom_rect(xmin="DBS",xmax="Month 1", color="gray",
            ymin=-Inf, ymax=-1, 
            fill="gray", alpha=0.15) +
  geom_rect(xmin="Month 1",xmax="Month 6", color="black",
            ymin=-Inf, ymax=-1, 
            fill="black", alpha=0.15) +
  annotate("text", x = 2.5, y = -5, 
           label = "OFF",fontface = "bold",
           color = "white",
           size = 3.5)+ 
  annotate("text", x = "Month 3", y = -5, 
           label = "ON",fontface = "bold",
           color = "white",size = 3.5) + 
  geom_vline(xintercept = "DBS", linewidth = 2) +
  #coord_cartesian(ylim= c(0,1.2) ) +
  geom_boxplot(data = ug.mood.means[ug.mood.means$sess %in% c("Baseline","DBS","Month 1", "Month 3", "Month 6"),],
               linewidth = 1.1,outlier.alpha = 0,
               show.legend = F) +
  geom_point(data = ug.mood.means[ug.mood.means$sess %in% c("Baseline","Month 1", "Month 3", "Month 6"),],
             size=2,color = "purple",alpha=.5,
             position = position_dodge2(width = .3),
             stroke = 1.75,aes( group = idx
             ),
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x =  element_blank(), y = "rating",
       title = "Mood [UG]",
       color = element_blank(),
       shape = element_blank()) +
  scale_shape_manual(values = c(16,4)) +
  stat_compare_means( 
    label="p.signif",
      comparisons = list(
        c("Baseline","Month 1"),
        c("Baseline","Month 3"),
        c("Baseline","Month 6")),
      method = "t.test",
      bracket.size = 0.5,
      tip.length = 0.02,
      step.increase = 0.08,label.x = c(2.5,2.5,2.5)
    )  
    

figure3_top = rl.logrt.sess+rl.c.sess +ug.logrt.sess+ ug.c.sess + ug.mood.sess  +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect",nrow = 1) & theme(legend.position = "none")


ggsave(res_dir / "figure3_top.png", 
       plot = figure3_top,
       device = "png",
       width = 15,          # Width in inches
       height = 4,         # Height in inches  
       dpi = 300)           # Resolution (300 DPI for publication quality)


############
# NT ~ Beh #
############

# Simplify changes in NT across sessions and task (by NT type)
rl.nt = rl.EST.Reward %>%
  filter(nt %in% c("DA","SE")) %>%
  group_by(idx,nt,stim) %>%
  summarise(est = mean(Oz)) %>%
  group_by(idx,nt) %>%
  summarise(dff = est[stim == "Post-Stim"] - est[stim == "Pre-Stim"]) %>%
  group_by(idx) %>%
  reframe(DA = dff[nt == "DA"],
         SE = dff[nt == "SE"])

ug.nt = ug.EST.Offer %>%
  filter(nt %in% c("DA","SE")) %>%
  group_by(idx,nt,stim) %>%
  summarise(est = mean(Oz)) %>%
  group_by(idx,nt) %>%
  summarise(dff = est[stim == "Post-Stim"] - est[stim == "Pre-Stim"]) %>%
  group_by(idx) %>%
  reframe(DA = dff[nt == "DA"],
          SE = dff[nt == "SE"])


rl.beh.change = 
  rl.beh.means %>%
  mutate(mRT = mLogRT) %>%
  select(!mLogRT) %>%
  pivot_wider(names_from = sess, values_from = c(mChoice, mRT)) %>%
  mutate(# Use Baseline as baseline
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
  select(c(idx,mChoice_PostStim_Change,mChoice_W1_Change,mChoice_M1_Change,mChoice_M2_Change,
           mChoice_M3_Change,mChoice_M4_Change,mChoice_M5_Change,mChoice_M6_Change,
           mRT_PostStim_Change,mRT_W1_Change,mRT_M1_Change,mRT_M2_Change,
           mRT_M3_Change,mRT_M4_Change,mRT_M5_Change,mRT_M6_Change,
  ))

rl.beh.nt = merge(rl.beh.change, rl.nt) 

rl.beh.nt.long = rl.beh.nt %>% 
  pivot_longer(cols = c(DA,SE), names_to = "NT", values_to = "ests") %>%
  pivot_longer(cols = c(mChoice_PostStim_Change,mChoice_W1_Change,mChoice_M1_Change,mChoice_M2_Change,
                        mChoice_M3_Change,mChoice_M4_Change,mChoice_M5_Change,mChoice_M6_Change,
                        mRT_PostStim_Change,mRT_W1_Change,mRT_M1_Change,mRT_M2_Change,
                        mRT_M3_Change,mRT_M4_Change,mRT_M5_Change,mRT_M6_Change), names_to = "change") %>%
  mutate(month = str_extract(change, "(?<=_)[^_]+(?=_)"),
         beh = ifelse(grepl("mChoice", change), "choice", "rt"))

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
  select(c(idx,mChoice_PostStim_Change,mChoice_W1_Change,mChoice_M1_Change,mChoice_M2_Change,
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


# UG Mood
ug.mood.change = 
  ug.mood.means %>%
  mutate(mRT = mLogRT) %>%
  select(!mLogRT) %>%
  pivot_wider(names_from = sess, values_from = c(mMood, mRT)) %>%
  mutate(# Use fMRI as baseline
    mMood_PostStim_Change = `mMood_Post-Stim` - mMood_Baseline,
    mMood_W1_Change = `mMood_Week 1` - mMood_Baseline,
    mMood_M1_Change = `mMood_Month 1` - mMood_Baseline,
    mMood_M2_Change = `mMood_Month 2` - mMood_Baseline,
    mMood_M3_Change = `mMood_Month 3` - mMood_Baseline,
    mMood_M4_Change = `mMood_Month 4` - mMood_Baseline,
    mMood_M5_Change = `mMood_Month 5` - mMood_Baseline,
    mMood_M6_Change = `mMood_Month 6` - mMood_Baseline,
    mRT_PostStim_Change = `mRT_Post-Stim` - mRT_Baseline,
    mRT_W1_Change = `mRT_Week 1` - mRT_Baseline,
    mRT_M1_Change = `mRT_Month 1` - mRT_Baseline,
    mRT_M2_Change = `mRT_Month 2` - mRT_Baseline,
    mRT_M3_Change = `mRT_Month 3` - mRT_Baseline,
    mRT_M4_Change = `mRT_Month 4` - mRT_Baseline,
    mRT_M5_Change = `mRT_Month 5` - mRT_Baseline,
    mRT_M6_Change = `mRT_Month 6` - mRT_Baseline
  ) %>%
  select(c(idx,mMood_PostStim_Change,mMood_W1_Change,mMood_M1_Change,mMood_M2_Change,
           mMood_M3_Change,mMood_M4_Change,mMood_M5_Change,mMood_M6_Change,
           mRT_PostStim_Change,mRT_W1_Change,mRT_M1_Change,mRT_M2_Change,
           mRT_M3_Change,mRT_M4_Change,mRT_M5_Change,mRT_M6_Change ))

ug.mood.nt = merge(ug.mood.change, ug.nt) 

ug.mood.nt.long = ug.mood.nt %>% 
  pivot_longer(cols = c(DA,SE), names_to = "NT", values_to = "ests") %>%
  pivot_longer(cols = c(mMood_PostStim_Change,mMood_W1_Change,mMood_M1_Change,mMood_M2_Change,
                        mMood_M3_Change,mMood_M4_Change,mMood_M5_Change,mMood_M6_Change,
                        mRT_PostStim_Change,mRT_W1_Change,mRT_M1_Change,mRT_M2_Change,
                        mRT_M3_Change,mRT_M4_Change,mRT_M5_Change,mRT_M6_Change), 
               names_to = "change") %>%
  mutate(month = str_extract(change, "(?<=_)[^_]+(?=_)"),
         beh = ifelse(grepl("mMood", change), "mood", "rt"),
         baseline = ifelse(grepl("Change2",change),"w1","fmri"))


# Run correlations and tidy the results

rl.beh.nt.long %>%
  group_by(NT,beh,month) %>%
  do(tidy(cor.test(.$value, .$ests, method = "spearman"))) %>%
  # Filter for significant results (p < 0.05)
  filter(p.value < .05)

# All SE~RT for reporting
rl.beh.nt.long %>%
  group_by(NT,beh,month) %>%
  do(tidy(cor.test(.$value, .$ests, method = "spearman"))) %>%
  # Filter for significant results (p < 0.05)
  filter(NT == "SE", beh == "rt")

ug.beh.nt.long %>%
  group_by(NT,beh,month) %>%
  do(tidy(cor.test(.$value, .$ests, method = "spearman"))) %>%
  # Filter for significant results (p < 0.05)
  filter(p.value < .05)

ug.mood.nt.long %>%
  group_by(NT,beh,month) %>%
  do(tidy(cor.test(.$value, .$ests, method = "spearman"))) %>%
  # Filter for significant results (p < 0.05)
  filter(p.value < .05)

# All DA ~ mood
ug.mood.nt.long %>%
  group_by(NT,beh,month) %>%
  do(tidy(cor.test(.$value, .$ests, method = "spearman"))) %>%
  # Filter for significant results (p < 0.05)
  filter(NT == "DA", beh == "mood")


rl.rt.SE =   rl.beh.nt.long %>%
  mutate(month = factor(month, levels=c("PostStim","W1","M1","M2","M3","M4","M5","M6"),
                        labels = c("Post-Stim","Week 1", "Month 1", "Month 2", "Month 3",
                                   "Month 4","Month 5","Month 6") ) ) %>%
  filter(NT == "SE",
         beh == "rt",
         month %in% c("Month 1", "Month 3", "Month 6"))

ug.mood.DA =   ug.mood.nt.long %>%
  mutate(month = factor(month, levels=c("PostStim","W1","M1","M2","M3","M4","M5","M6"),
                        labels = c("Post-Stim","Week 1", "Month 1", "Month 2", "Month 3",
                                   "Month 4","Month 5","Month 6") ) ) %>%
  filter(NT == "DA",
         beh == "mood",
         month %in% c("Month 1", "Month 3", "Month 6")) 


fig.rl.rt.SE.simple = 
  rl.rt.SE %>%
  ggplot(aes(x = ests, y = value))+
  theme_pubr(base_size = 16) +
  facet_wrap(~month, nrow = 1) +
  stat_smooth(method = "lm", 
              alpha = .25,
              linetype = "solid",
              linewidth = 2.1,
              color = NT_colors$color[NT_colors$id == "5-HT"], 
              fill = NT_colors$color[NT_colors$id == "5-HT"] ) +
  geom_point(size = 3,stroke = 1.5) +
  labs(title = "Response time [RL]",
       x = "5-HT Change\nPost - Pre [Z]",
       # x = element_blank(),
       y = "Change from Baseline",
       # y = "Change from Baseline\nResponse Times (log-sec)",
       shape = element_blank(),
       color = element_blank()) +
  guides(color = "none",shape = "none") +
  coord_cartesian(xlim=c(-10,10),ylim = c(-3.25, 1.5)) +
  theme(
    panel.spacing = unit(1, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = NA, fill = NA)) +
  stat_cor(method="spearman",aes(label =paste0(cut(..p.., 
                                                   breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                   labels = c("'****'", "'***'", "'**'", "'*'","'ns'"))
  )),
  label.x = 0, size = c(8,12,12), label.y = 1)


fig.ug.m.DA.simple =
  ug.mood.DA %>%
  ggplot(aes(x = ests, y = value))+
  theme_pubr(base_size = 16) +
  facet_wrap(~month, nrow = 1) +
  stat_smooth(method = "lm", 
              alpha = .25,
              linetype = "solid",
              linewidth = 2.1,
              color = NT_colors$color[NT_colors$id == "DA"], 
              fill = NT_colors$color[NT_colors$id == "DA"], ) +
  geom_point(size = 3,stroke = 1.5) +
  labs(title = "Mood [UG]",
       x = "DA Change\nPost - Pre [Z]",
       y = "Change from Baseline",
       shape = element_blank(),
       color = element_blank()) +
  coord_cartesian(xlim=c(-10,10), ylim = c(-55,55)) +
  theme(
    panel.spacing = unit(1, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = NA, fill = NA)) +
  stat_cor(method="spearman",aes(label =paste0(cut(..p.., 
                                                   breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                   labels = c("'****'", "'***'", "'**'", "'*'", "'ns'"))
  )),label.x = -5, size = c(12,12,8), label.y = 43)


figure3_bot = (fig.rl.rt.SE.simple + fig.ug.m.DA.simple) +
  plot_annotation(tag_levels = list(c('f', 'g'), '1')) +
  plot_layout(guides = "collect",nrow = 1) &
  theme(legend.position = 'right') 

ggsave(res_dir / "figure3_bot.png", 
       plot = figure3_bot,
       device = "png",
       width = 15,          # Width in inches
       height = 4,         # Height in inches  
       dpi = 300)           # Resolution (300 DPI for publication quality)

