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
# 2025/09/19    Blair Shevlin                           generating alternative NT metrics
# 2025/10/08    Blair Shevlin                           outputing raw NT
# 2025/11/25    Blair Shevlin                           adapted for raw model outputs

# Goal of this code is to process the NT data for the UG and RL tasks into usable formats while combining with behavioral data.

rm(list = ls())

library(tidyverse)
library(fs)
library(here)
library(RColorBrewer)

# Primary directories
dir = path(here())

beh_dir = dir / "data" /"behavior"
in_nt = dir / "data" / "nt" / "raw"
out_nt = dir / "data" / "nt" / "processed"

# Behavior
ug_beh = read.csv(file = beh_dir / "ug-data_deid_07-10-25.csv")
rl_beh = read.csv(file = beh_dir / "rl-data_deid_07-10-25.csv")

# Colors
id_colors = data.frame(idx = ug_beh$idx %>% unique(),
                       color = brewer.pal(n=10,name = "Paired"))

# Restrict to OR session
ug_OR = ug_beh %>%
  filter(sess %in% c("Pre-Stim","Post-Stim"))
rl_OR = rl_beh %>%
  filter(sess %in% c("Pre-Stim","Post-Stim"))

# Load NT data
ug_files = list.files(path = in_nt, pattern = "UG")
rl_files = list.files(path = in_nt, pattern = "RL")

# Keep only files with "raw-model-output"
ug_files <- ug_files[grepl("raw-model-output", ug_files)]
rl_files <- rl_files[grepl("raw-model-output", rl_files)]

# Remove RL files from participants 1&2 who were misalligned
rl_files <- rl_files[!grepl("deid-1_", rl_files) & !grepl("deid-2_", rl_files)]

# Empty container
ug_nt = NULL
rl_nt = NULL

# Combine UG NT data and merge with behavior
for (f in 1:length(ug_files) ) {
    print(paste0("Processing file: ", ug_files[f]))
    tmp_ug = read.csv(paste0(in_nt,"/",ug_files[f])) %>%
    mutate(sess = stim)

     tmp_beh = ug_OR %>%
      filter(idx == unique(tmp_ug$idx),
            sess == unique(tmp_ug$stim)) %>%
      mutate(prev_offer_raw = ifelse(trial==1,NA,lag(offer)),
             prev_offer = ifelse(is.na(prev_offer_raw),NA,
                                 ifelse(prev_offer_raw > offer, "-PE","+PE"))
      )

    tmp_ug_m = merge(tmp_ug, tmp_beh, by = c("idx","trial","sess"), all.x = TRUE)

    ug_nt = rbind(ug_nt, tmp_ug_m)
}

# Combine RL NT data and merge with behavior
for (f in 1:length(rl_files) ) {
    print(paste0("Processing file: ", rl_files[f]))

    tmp_rl = read.csv(paste0(in_nt,"/",rl_files[f])) %>%
    mutate(sess = stim)

     tmp_beh = rl_OR %>%
        filter(idx == unique(tmp_rl$idx),
            sess == unique(tmp_rl$stim)) %>%
        mutate(prev_rew_raw = ifelse(trial_within_block==1,NA,lag(outcome)),
                    prev_rew = ifelse(is.na(prev_rew_raw),NA,
                                    ifelse(prev_rew_raw == 1,
                                            ifelse(cond %in% c("Mixed", "Positive"),10,0),
                                            ifelse(cond %in% c("Mixed", "Negative"),-10,0)
                                    )
                    )
            )

    tmp_rl_m = merge(tmp_rl, tmp_beh, by = c("idx","trial","sess"), all.x = TRUE)


    rl_nt = rbind(rl_nt, tmp_rl_m)

#    print(    nrow(rl_nt[is.na(rl_nt$rt),]) ) # Check for missing RTs

}

# Z-score and re-number samples to be consistent with other format
ug_proc = ug_nt %>%
  pivot_longer(cols = c("SE","DA","NE"), values_to = "nM",
               names_to = "nt") %>%
  filter(rt < 10) %>%
  group_by(idx,event,nt) %>%
  mutate(nMz = scale(nM)[,1],
        sample = sample + 31) %>%
  ungroup()
rl_proc = rl_nt %>%
  pivot_longer(cols = c("SE","DA","NE"), values_to = "nM",
               names_to = "nt") %>%
  filter(rt < 10) %>%
  group_by(idx,event,nt) %>%
  mutate(nMz = scale(nM)[,1],
        sample = sample + 51 ) %>%
  ungroup()

save(ug_proc,rl_proc,
     file = out_nt / "UG_RL_NT-Raw_Model_Output_11-25-25.RData")

# Construct trial-level estimates
# 9.16.2025 - Extending samples
# 9.23.2025 - Adding metrics with additional samples
# 11.25.2025 - Adapting for raw model outputs

rl.EST = rl_proc %>%
  filter(!is.na(nM),
         ) %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         nt = factor(nt, levels = c("SE","DA","NE")),
         event = factor(event,levels=c("Cue","Choice","Reward"))) %>%
  group_by(idx,stim,trial,nt,event,cond,block,opt,rt,prev_rew_raw,trial_within_block,outcome,reversal) %>%
  reframe(# Overall estimate (window)
            O = sum(nM[(sample<=56 & sample>=51)]),
            Oz = sum(nMz[(sample<=56 & sample>=51)]),
            # Relative estimate (window)
            R = sum(nM[(sample<=56 & sample>51)]) - (nM[sample==51] * length(c(52:56))),
            Rz = sum(nMz[(sample<=56 & sample>51)]) - (nMz[sample==51] * length(c(52:56))),
            # Peak estimate over trial
            P = max(nM[(sample<=56 & sample>=51)]), 
            Pz = max(nMz[(sample<=56 & sample>=51)]),
            # Avg estimate over trial
            M = mean(nM[(sample<=56 & sample>=51)]),
            Mz = mean(nMz[(sample<=56 & sample>=51)]),
            # Estimate over 10 samples
            O10 = sum(nM[(sample<=61 & sample>=51)]),
            Oz10 = sum(nMz[(sample<=61 & sample>=51)]),
            R10 = sum(nM[(sample<=61 & sample>51)]) - (nM[sample==51] * length(c(52:61))),
            Rz10 = sum(nMz[(sample<=61 & sample>51)]) - (nMz[sample==51] * length(c(52:61))),
            P10 = max(nM[(sample<=61 & sample>=51)]),
            Pz10 = max(nMz[(sample<=61 & sample>=51)]),
            M10 = mean(nM[(sample<=61 & sample>=51)]),
            Mz10 = mean(nMz[(sample<=61 & sample>=51)])
            )

ug.EST = ug_proc %>%
  filter(!is.na(nM)
         ) %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         nt = factor(nt, levels = c("SE","DA","NE")),
         event = factor(event,levels=c("Offer","Choice")),
         offer_bin = ifelse(offer <4,"Low",
                            ifelse(offer>6,"High","Middle")),
         offer_bin = factor(offer_bin,
                            levels = c("Low","Middle","High"))) %>%
  group_by(idx,stim,trial,nt,event,rej,rt,offer,prev_offer_raw,offer_bin) %>%
  reframe(O = sum(nM[(sample<=36 & sample>=31)]),
            Oz = sum(nMz[(sample<=36 & sample>=31)]),
            R = sum(nM[(sample<=36 & sample>31)]) - (nM[sample==31] * length(c(32:36))),
            Rz = sum(nMz[(sample<=36 & sample>31)]) - (nMz[sample==31] * length(c(32:36))),
            # Peak estimate over trial
            P = max(nM[(sample<=36 & sample>=31)]), 
            Pz = max(nMz[(sample<=36 & sample>=31)]),
            # Avg estimate over trial
            M = mean(nM[(sample<=36 & sample>=31)]),
            Mz = mean(nMz[(sample<=36 & sample>=31)]),
            # Total estimate over 10 samples
            O10 = sum(nM[(sample<=41 & sample>=31)]),
            Oz10 = sum(nMz[(sample<=41 & sample>=31)]),
            R10 = sum(nM[(sample<=41 & sample>31)]) - (nM[sample==31] * length(c(32:41))),
            Rz10 = sum(nMz[(sample<=41 & sample>31)]) - (nMz[sample==31] * length(c(32:41))),
            P10 = max(nM[(sample<=41 & sample>=31)]),
            Pz10 = max(nMz[(sample<=41 & sample>=31)]),
            M10 = mean(nM[(sample<=41 & sample>=31)]),
            Mz10 = mean(nMz[(sample<=41 & sample>=31)])
  )

# Add colors to everything
rl.EST = merge(rl.EST,id_colors)
ug.EST = merge(ug.EST, id_colors)

save(rl.EST,ug.EST,
      file = out_nt / "UG_RL_NT-Continuous_Raw-Model-Output_11-25-25.RData")

