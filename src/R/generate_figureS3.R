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
library(ggpubr)
library(here)
library(fs)
library(here)
library(rstatix)

# Paths
dir = path(here())
clin_dir = dir / "data" / "clinical" # Updated
res_dir = dir / "results" # Updated

# Load clinical data
cl.df = read.csv(clin_dir / "clinical-data_deid_07-10-25.csv" ) 


cl.O.fig = cl.df %>%
  filter(!session %in% c('fmri')) %>%
  mutate(id = factor(idx),
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
    cohort = ifelse(id %in% c(1:5),"RC+S","PC")
  )%>%
  ungroup() %>% 
  filter(!sess_fig %in% c("Month 2", "Month 4", "Month 5")) %>%
  as.data.frame()

test.cohort <-     
  cl.O.fig %>%
  filter(sess_fig!="DBS")%>%
  group_by(sess_fig) %>%
  t_test(HDRS ~ cohort) %>%
  add_significance("p")  %>%
  add_xy_position(x = "sess_fig", dodge = 1) %>%
  mutate(y.position = 30,
         xmin = c(0.75, 2.75, 3.75, 4.75, 5.75),
         xmax = c(1.25, 3.25, 4.25, 5.25, 6.25))

fig.hdrs.timeline.cohort =
cl.O.fig %>%
  ggplot(aes(x = sess_fig, y = HDRS)) +
  theme_pubr(base_size = 14) +
  geom_vline(xintercept = "DBS", linewidth = 2) +
  stat_summary(geom = "line",
               aes(color = cohort, group = cohort),
               position = position_dodge2(width = .75),
               linewidth = 1.5 )+ 
  stat_summary(aes( shape = cohort, color = cohort),
               position = position_dodge2(width = .75),
               size = 1.5, linewidth = 1.5 )+  
  geom_point(data = cl.O.fig[cl.O.fig$sess_fig!="DBS",],
             size=2, alpha = .5,#color = "purple"
             position = position_dodge2(width = .75),
             stroke = 1.75,aes( shape = cohort, color = cohort),
  ) +
  labs(x = element_blank(),
       y = "HDRS-17",
       shape = "Cohort",
       color = "Cohort") +
  scale_shape_manual(values = c(16,15)) +
  scale_color_brewer(type = "qual",palette = 2) +
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=0.75),
        legend.position = c(0.1, 0.2))+
  stat_pvalue_manual(
    test.cohort,  tip.length = 0) 



ggsave(res_dir / "Extended-Data_Figure3.png", 
       plot = fig.hdrs.timeline.cohort,
       device = "png",
       width = 6.5,          # Width in inches
       height = 3.5,         # Height in inches  
       dpi = 300)           # Resolution (300 DPI for publication quality)
