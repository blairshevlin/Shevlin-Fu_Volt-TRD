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
# 2025/10/08    Blair Shevlin                           include updated NT data, OR session
# 2025/11/04    Blair Shevlin                           minor edits to finalize
# 2025/11/14    Blair Shevlin                           split into seperate figures

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
clin_dir = dir / "data" / "clinical"

# Load clinical data
cl.df = read.csv(clin_dir / "clinical-data_deid_07-10-25.csv" ) 

# Subject-level NT Data
load(nt_dir / "UG_RL_NT-Continuous_9-23-25.RData")

# Extract relevant dataframes
rl.EST.Reward = rl.EST %>%
  filter(event == "Reward") %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))

ug.EST.Offer = ug.EST %>%
  filter(event == "Offer") %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim"))) 

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
  
  return(t_test)
}

calculate_star_positions <- function(data, value_col, bracket_height_offset = 3) {
  max_val <- max(data[[value_col]], na.rm = TRUE)
  
  # Position stars just above where brackets would be
  # Adjust bracket_height_offset based on your bracket positions
  positions <- max_val + bracket_height_offset + c(2, 7, 12)
  return(positions)
}

significance_labels_ug_c_alt <- sapply(c("Pre-Stim","Post-Stim"), function(comp_level) {
  calculate_significance(ug.beh.means, "mChoice","sess" , "Baseline", 
                         comp_level, "idx", paired = TRUE)
}
)

significance_labels_rl_c_alt <- sapply(c("Pre-Stim","Post-Stim"), function(comp_level) {
  calculate_significance(rl.beh.means, "mChoice","sess" , "Baseline", 
                         comp_level, "idx", paired = TRUE)
}
)

significance_labels_ug_logrt_alt <- sapply(c("Pre-Stim","Post-Stim"), function(comp_level) {
  calculate_significance(ug.beh.means, "mLogRT","sess" , "Baseline", 
                         comp_level, "idx", paired = TRUE)
}
)

significance_labels_rl_logrt_alt <- sapply(c("Pre-Stim","Post-Stim"), function(comp_level) {
  calculate_significance(rl.beh.means, "mLogRT","sess" , "Baseline", 
                         comp_level, "idx", paired = TRUE)
}
)

significance_labels_ug_mood_alt <- sapply(c("Pre-Stim","Post-Stim"), function(comp_level) {
  calculate_significance(ug.mood.means, "mMood","sess" , "Baseline", 
                         comp_level, "idx", paired = TRUE)
}
)

sessions = c("Baseline", "Pre-Stim", "Post-Stim", "DBS","Month 1","Month 3","Month 6")

rl.beh.means.plot = rl.beh.means %>%
    filter(sess %in% sessions) %>%
    mutate(sess= factor(sess, levels = sessions))
ug.beh.means.plot = ug.beh.means %>%
    filter(sess %in% sessions) %>%
    mutate(sess= factor(sess, levels = sessions))
ug.mood.means.plot = ug.mood.means %>%
    filter(sess %in% sessions) %>%
    mutate(sess= factor(sess, levels = sessions))

rl.c.sess = 
  rl.beh.means.plot %>%
  mutate(sess= factor(sess, levels = sessions)) %>%
  ggplot(aes(x = sess, y = mChoice))+
  theme_pubr(base_size = 14) +
  geom_vline(xintercept = "DBS", linewidth = 2) +
    geom_boxplot(data = rl.beh.means.plot,
               linewidth = 1.1,outlier.alpha = 0,
               show.legend = F) +
  geom_point(data = rl.beh.means.plot,
             size=2,color = "purple",alpha=.5,
             position = position_dodge2(width = .3),
             stroke = 1.75,
             aes( group = idx
             ),shape = 18
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
      c("Baseline", "Post-Stim"),
      c("Pre-Stim", "Post-Stim") ),  
    method = "wilcox",
    paired = TRUE,
    bracket.size = 0.5,
    tip.length = 0.02,
    step.increase = 0.08,label.x = 2.5
  ) 

ug.c.sess = 
  ug.beh.means.plot %>%
  ggplot(aes(x = sess, y = mChoice))+
  theme_pubr(base_size = 14) +
  geom_vline(xintercept = "DBS", linewidth = 2) +
  #coord_cartesian(ylim= c(0,1.2) ) +
  geom_boxplot(data = ug.beh.means.plot,
               linewidth = 1.1,outlier.alpha = 0,
               show.legend = F) +
  geom_point(data = ug.beh.means.plot,
             size=2,color = "purple",alpha=.5,
             position = position_dodge2(width = .3),
             stroke = 1.75,aes( group = idx
             ),shape = 18
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
      c("Baseline", "Post-Stim"),
      c("Pre-Stim", "Post-Stim") ), 
    method = "wilcox",    paired = TRUE,
    bracket.size = 0.5,
    tip.length = 0.02,
    step.increase = 0.08,label.x = 2.5
  ) 

ug.logrt.sess = 
  ug.beh.means.plot %>%
  ggplot(aes(x = sess, y = mLogRT))+
  theme_pubr(base_size = 14) +
  geom_vline(xintercept = "DBS", linewidth = 2) +
  geom_boxplot(data = ug.beh.means.plot,
               linewidth = 1.1,outlier.alpha = 0,
               show.legend = F) +
  geom_point(data = ug.beh.means.plot,
             size=2,color = "purple",alpha=.5,
             position = position_dodge2(width = .3),
             stroke = 1.75,aes( group = idx
             ),shape = 18
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x =  element_blank(), y = "log sec",
       color = element_blank(),title = "Response time [UG]",
       shape = element_blank()) +
  scale_shape_manual(values = c(16,4)) +
  stat_compare_means( 
    label="p.signif",
    comparisons = list(
      c("Baseline", "Post-Stim"),
      c("Pre-Stim", "Post-Stim") ), 
    method = "wilcox",    paired = TRUE,

    bracket.size = 0.5,
    tip.length = 0.02,
    step.increase = 0.08,label.x =2.5
  )  

rl.logrt.sess = 
  rl.beh.means.plot %>%
  ggplot(aes(x = sess, y = mLogRT))+
  theme_pubr(base_size = 14) +
  geom_vline(xintercept = "DBS", linewidth = 2) +
  #coord_cartesian(ylim= c(-0.6,2) ) +
  geom_boxplot(data = rl.beh.means.plot,
               linewidth = 1.1,outlier.alpha = 0,
               show.legend = F) +
  geom_point(data = rl.beh.means.plot,
             size=2,color = "purple",alpha=.5,
             position = position_dodge2(width = .3),
             stroke = 1.75,aes( group = idx
             ),shape = 18
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
      c("Baseline", "Post-Stim"),
      c("Pre-Stim", "Post-Stim") ), 
    method = "wilcox",    paired = TRUE,

    bracket.size = 0.5,
    tip.length = 0.02,
    step.increase = 0.08,label.x = 2.5
  )  

ug.mood.sess = 
  ug.mood.means.plot %>%
  ggplot(aes(x = sess, y = mMood))+
  theme_pubr(base_size = 14) +
  geom_vline(xintercept = "DBS", linewidth = 2) +
  #coord_cartesian(ylim= c(0,1.2) ) +
  geom_boxplot(data = ug.mood.means.plot,
               linewidth = 1.1,outlier.alpha = 0,
               show.legend = F) +
  geom_point(data = ug.mood.means.plot,
             size=2,color = "purple",alpha=.5,
             position = position_dodge2(width = .3),
             stroke = 1.75,aes( group = idx
             ),shape = 18
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x =  element_blank(), y = "rating",
       title = "Mood [UG]",
       color = element_blank(),
       shape = element_blank()) +
  scale_shape_manual(values = c(16,4)) +
  stat_compare_means( 
    label= "p.signif",
    comparisons = list(
      c("Baseline", "Post-Stim"),
      c("Pre-Stim", "Post-Stim") ), 
      method = "wilcox",    paired = TRUE,

      bracket.size = 0.5,
      tip.length = 0.02,
      step.increase = 0.08 #,label.x = c(2.5,2.5,2.5)
    )  
    

figure3 = rl.logrt.sess+rl.c.sess+ ug.logrt.sess+ug.c.sess + ug.mood.sess  +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect",nrow = 1) & theme(legend.position = "none")


ggsave(res_dir / "figureS3_raw.png", 
       plot = figure3,
       device = "png",
       width = 15,          # Width in inches
       height = 4,         # Height in inches  
       dpi = 300)           # Resolution (300 DPI for publication quality)