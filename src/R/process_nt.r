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

# Z-score
ug_proc = ug_nt %>%
  group_by(idx,event,nt) %>%
  mutate(nMz = scale(nM)[,1]) %>%
  ungroup()
rl_proc = rl_nt %>%
  group_by(idx,event,nt) %>%
  mutate(nMz = scale(nM)[,1]) %>%
  ungroup()

# Construct trial-level estimates
rl.EST = rl_proc %>%
  filter(!is.na(nM),
         rt < 10,
         ) %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         nt = factor(nt, levels = c("SE","DA","NE")),
         event = factor(event,levels=c("Cue","Choice","Reward"))) %>%
  group_by(idx,stim,trial,nt,event,cond,block,opt,rt,prev_rew_raw,trial_within_block,outcome,reversal) %>%
  summarise(# Overall estimate (window)
            O = sum(nM[(sample<56 & sample>51)]),
            Oz = sum(nMz[(sample<56 & sample>51)]),
            # Relative estimate (window)
            R = sum(nM[(sample<56 & sample>51)]) - (nM[sample==50] * length(c(51:56))),
            Rz = sum(nMz[(sample<56 & sample>51)]) - (nMz[sample==50] * length(c(51:56))),
            # Peak estimate over trial
            P = max(nM), 
            Pz = max(nMz),
            # Avg estimate over trial
            M = mean(nM),
            Mz = mean(nMz),
            # Total estimate over trial
            Tot = sum(nM),
            Totz = sum(nMz)
            )

ug.EST = ug_proc %>%
  filter(!is.na(nM),
         rt < 10,
         ) %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")),
         nt = factor(nt, levels = c("SE","DA","NE")),
         event = factor(event,levels=c("Offer","Choice")),
         offer_bin = ifelse(offer <4,"Low",
                            ifelse(offer>6,"High","Middle")),
         offer_bin = factor(offer_bin,
                            levels = c("Low","Middle","High"))) %>%
  group_by(idx,stim,trial,nt,event,rej,rt,offer,prev_offer_raw,offer_bin) %>%
  summarise(O = sum(nM[(sample<36 & sample>31)]),
            Oz = sum(nMz[(sample<36 & sample>31)]),
            R = sum(nM[(sample<36 & sample>31)]) - (nM[sample==30] * length(c(31:36))),
            Rz = sum(nMz[(sample<36 & sample>31)]) - (nMz[sample==30] * length(c(31:36))),
            # Peak estimate over trial
            P = max(nM), 
            Pz = max(nMz),
            # Avg estimate over trial
            M = mean(nM),
            Mz = mean(nMz),
            # Total estimate over trial
            Tot = sum(nM),
            Totz = sum(nMz)
            )

# Generate summary data - focusing only on reward events (RL) and offer events (UG)
rl.EST.Reward = rl.EST %>%
  filter(event == "Reward") %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))

rl.EST.Reward = merge(rl.EST.Reward,id_colors)

ug.EST.Offer = ug.EST %>%
  filter(event == "Offer") %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim"))) 

ug.EST.Offer = merge(id_colors, ug.EST.Offer)

save(rl.EST.Reward,rl.EST.Reward, file = out_nt / "UG_RL_NT-Continuous_7-14-25.RData")
