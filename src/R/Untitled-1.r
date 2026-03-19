###########################
# Trial-level regressions #
###########################

rl.EST.pr %>% head()

rl_trial_data <- merge(rl.EST.pr, cl.df.baseline)  %>%
  mutate(
    # Create meaningful predictors
    block_type_factor = block_type,
    log_rt = log(rt),
    rt_z = scale(rt)[,1],
    trial_z = scale(trial)[,1],
    
    # Create reversal indicators
    reversal_trial = reversal,
    post_reversal = ifelse(reversal == 1 | 
                          lag(reversal, 1, default = 0) == 1 | 
                          lag(reversal, 2, default = 0) == 1 | 
                          lag(reversal, 3, default = 0) == 1 | 
                          lag(reversal, 4, default = 0) == 1, 1, 0), # reversal trial + 4 trials after
    pre_reversal = ifelse(lead(reversal, 1, default = 0) == 1 | 
                         lead(reversal, 2, default = 0) == 1 | 
                         lead(reversal, 3, default = 0) == 1 | 
                         lead(reversal, 4, default = 0) == 1 | 
                         lead(reversal, 5, default = 0) == 1, 1, 0), # 5 trials before reversal
    
    # Error trials (assuming optimal choice exists)
    error_trial = ifelse(opt == 0, 1, 0), # Did not choose optimal option
    
    # RT bins for analysis
    rt_slow = ifelse(rt > median(rt, na.rm = TRUE), 1, 0),
    
    # Depression severity
    HDRS_z = scale(HDRS)[,1] ,

    # Prediction error: current reward - previous reward (simple version)
    pred_error = rew - prev_rew,
    # Categorize prediction errors
    pred_error_type = case_when(
      pred_error > 0 ~ "Positive_PE",   # Better than expected
      pred_error < 0 ~ "Negative_PE",   # Worse than expected  
      pred_error == 0 ~ "No_PE",        # As expected
      TRUE ~ "Missing"
    ),
    pred_error_type = factor(pred_error_type, levels = c("No_PE", "Negative_PE", "Positive_PE"))
  )


# Model 1: Basic trial-level effects
rl.DA.trial.lm1 <- lmer(Oz ~ outcome_f * stim + 
                       block_type_factor * stim + 
                       log_rt * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = rl_trial_data %>% filter(nt == "DA"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
rl.SE.trial.lm1 <- lmer(Oz ~ outcome_f * stim + 
                       block_type_factor * stim + 
                       log_rt * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = rl_trial_data %>% filter(nt == "SE"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
rl.NE.trial.lm1 <- lmer(Oz ~ outcome_f * stim + 
                       block_type_factor * stim + 
                       log_rt * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = rl_trial_data %>% filter(nt == "NE"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))

#Model 2: Error and reversal effects
rl.DA.trial.lm2 <- lmer(Oz ~ error_trial * stim + 
                       reversal_trial * stim + 
                       post_reversal * stim +
                       rt_slow * stim +
                       trial_z * stim +
                       (1 + stim | idx),
                     data = rl_trial_data %>% filter(nt == "DA"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
rl.SE.trial.lm2 <- lmer(Oz ~ error_trial * stim + 
                       reversal_trial * stim + 
                       post_reversal * stim +
                       rt_slow * stim +
                       trial_z * stim +
                       (1 + stim | idx),
                     data = rl_trial_data %>% filter(nt == "SE"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
rl.NE.trial.lm2 <- lmer(Oz ~ error_trial * stim + 
                       reversal_trial * stim + 
                       post_reversal * stim +
                       rt_slow * stim +
                       trial_z * stim +
                       (1 + stim | idx),
                     data = rl_trial_data %>% filter(nt == "NE"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))

#Model 3: Current and previous reward effects
rl.DA.trial.lm3 <- lmer(Oz ~ rew_f * stim + 
                    prev_rew_f * stim + 
                    rew_f * prev_rew_f * stim +  # interaction between current and previous reward
                    trial_z * stim +
                    (1 + stim | idx),
                  data = rl_trial_data %>% filter(nt == "DA", !is.na(prev_rew)),
                  REML = FALSE,
                  control = lmerControl(optimizer = "bobyqa"))
rl.SE.trial.lm3 <- lmer(Oz ~ rew_f * stim + 
                    prev_rew_f * stim + 
                    rew_f * prev_rew * stim +  # interaction between current and previous reward
                    trial_z * stim +
                    (1 + stim | idx),
                  data = rl_trial_data %>% filter(nt == "SE", !is.na(prev_rew)),
                  REML = FALSE,
                  control = lmerControl(optimizer = "bobyqa"))
rl.NE.trial.lm3 <- lmer(Oz ~ rew_f * stim + 
                    prev_rew_f * stim + 
                    rew_f * prev_rew_f * stim +  # interaction between current and previous reward
                    trial_z * stim +
                    (1 + stim | idx),
                  data = rl_trial_data %>% filter(nt == "NE", !is.na(prev_rew)),
                  REML = FALSE,
                  control = lmerControl(optimizer = "bobyqa"))

# Model 4: Prediction error effects (explicitly modeled)
rl.DA.trial.lm4 <- lmer(Oz ~ pred_error_type * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = rl_trial_data %>% filter(nt == "DA", !is.na(prev_rew)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
rl.SE.trial.lm4 <- lmer(Oz ~ pred_error_type * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = rl_trial_data %>% filter(nt == "SE", !is.na(prev_rew)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
rl.NE.trial.lm4 <- lmer(Oz ~ pred_error_type * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = rl_trial_data %>% filter(nt == "NE", !is.na(prev_rew)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))

## UG ##
########
ug.EST.pr %>% head()

ug_trial_data <- merge(ug.EST.pr, cl.df.baseline)  %>%
  mutate(
    # Create meaningful predictors
    outcome_binary = 1*(rej==0),  # Accept = 1, Reject = 0
    offer_bin_f = factor(offer_bin_f, levels = c("Medium","Low","High")),
    prev_offer_bin_f = ifelse(prev_offer_raw < 4, "Low", ifelse(prev_offer_raw > 7, "High", "Medium")),
    prev_offer_bin_f = factor(prev_offer_bin_f, levels = c("Medium","Low","High")),
    log_rt = log(rt),
    rt_z = scale(rt)[,1],
    trial_z = scale(trial)[,1],
    
    # Accept trials
    accept_trial = ifelse(outcome_binary == 1, 1, 0), # Accepted
    
    # RT bins for analysis
    rt_slow = ifelse(rt > median(rt, na.rm = TRUE), 1, 0),

    # Previous offer and outcome
    offer_change = offer - prev_offer_raw,

    # Categorize prediction errors
    offer_change_type = case_when(
      offer_change > 0 ~ "Offer_Better",   # Better than expected
      offer_change < 0 ~ "Offer_Worse",   # Worse than expected  
      offer_change == 0 ~ "Same_Offer",        # As expected
      TRUE ~ "Missing"
    ),
    offer_change_type = factor(offer_change_type, levels = c("Same_Offer", "Offer_Worse", "Offer_Better")),

    # Depression severity
    HDRS_z = scale(HDRS)[,1] 
  )

# Model 1: Basic trial-level effects
ug.DA.trial.lm1 <- lmer(Oz ~ offer_bin_f * stim + 
                       log_rt * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "DA"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
ug.SE.trial.lm1 <- lmer(Oz ~ offer_bin_f * stim + 
                       log_rt * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "SE"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
ug.NE.trial.lm1 <- lmer(Oz ~ offer_bin_f * stim + 
                       log_rt * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "NE"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))

#Model 2: Error and reversal effects
ug.DA.trial.lm2 <- lmer(Oz ~ accept_trial * stim + 
                       rt_slow * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "DA"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
ug.SE.trial.lm2 <- lmer(Oz ~ accept_trial * stim + 
                       rt_slow * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "SE"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
ug.NE.trial.lm2 <- lmer(Oz ~ accept_trial * stim + 
                       rt_slow * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "NE"),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))

#Model 3: Current and previous reward effects
ug.DA.trial.lm3 <- lmer(Oz ~ offer_bin_f * stim + 
                    prev_offer_bin_f * stim + 
                    offer_bin_f * prev_offer_bin_f * stim +  # interaction between current and previous reward
                    trial_z * stim +
                    (1 + stim | idx),
                  data = ug_trial_data %>% filter(nt == "DA", !is.na(prev_offer_raw)),
                  REML = FALSE,
                  control = lmerControl(optimizer = "bobyqa"))
ug.SE.trial.lm3 <- lmer(Oz ~ offer_bin_f * stim + 
                    prev_offer_bin_f * stim + 
                    offer_bin_f * prev_offer_bin_f * stim +  # interaction between current and previous reward
                    trial_z * stim +
                    (1 + stim | idx),
                  data = ug_trial_data %>% filter(nt == "SE", !is.na(prev_offer_raw)),
                  REML = FALSE,
                  control = lmerControl(optimizer = "bobyqa"))
ug.NE.trial.lm3 <- lmer(Oz ~ offer_bin_f * stim + 
                    prev_offer_bin_f * stim + 
                    offer_bin_f * prev_offer_bin_f * stim +  # interaction between current and previous reward
                    trial_z * stim +
                    (1 + stim | idx),
                  data = ug_trial_data %>% filter(nt == "NE", !is.na(prev_offer_raw)),
                  REML = FALSE,
                  control = lmerControl(optimizer = "bobyqa"))

# Model 4: Prediction error/offer change effects (explicitly modeled)
ug.DA.trial.lm4 <- lmer(Oz ~ offer_change_type * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "DA", !is.na(prev_offer_raw)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
ug.SE.trial.lm4 <- lmer(Oz ~ offer_change_type * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "SE", !is.na(prev_offer_raw)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
ug.NE.trial.lm4 <- lmer(Oz ~ offer_change_type * stim + 
                       trial_z * stim +
                       (1 + stim | idx),
                     data = ug_trial_data %>% filter(nt == "NE", !is.na(prev_offer_raw)),
                     REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))


###################################
# Summary Behavioral Analyses     #
###################################

# Reversal Learning Task Summary Statistics
rl_behavioral_summary <- rl_trial_data %>%
  group_by(idx, stim) %>%
  summarise(
    mean_rt = mean(rt, na.rm = TRUE),
    median_rt = median(rt, na.rm = TRUE),
    sd_rt = sd(rt, na.rm = TRUE),
    error_rate = mean(error_trial, na.rm = TRUE),
    optimal_choice_rate = mean(opt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(stim = ifelse(stim == "Pre-Stim", "pre", "post")) %>%
  pivot_wider(
    names_from = stim,
    values_from = c(mean_rt, median_rt, sd_rt, error_rate, optimal_choice_rate),
    names_sep = "_"
  )

# Ultimatum Game Summary Statistics  
ug_behavioral_summary <- ug_trial_data %>%
  group_by(idx, stim) %>%
  summarise(
    mean_rt = mean(rt, na.rm = TRUE),
    median_rt = median(rt, na.rm = TRUE),
    sd_rt = sd(rt, na.rm = TRUE),
    accept_rate = mean(accept_trial, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(stim = ifelse(stim == "Pre-Stim", "pre", "post")) %>%
  pivot_wider(
    names_from = stim,
    values_from = c(mean_rt, median_rt, sd_rt, accept_rate),
    names_sep = "_"
  )

# Summary tables for manuscript
create_summary_table <- function(data, task_name) {

  data =rl_behavioral_summary
  if(task_name == "Reversal Learning") {
    summary_stats <- data %>%
      summarise(
        across(c(mean_rt_pre, mean_rt_post, error_rate_pre, error_rate_post, 
                optimal_choice_rate_pre, optimal_choice_rate_post), 
               list(mean = ~mean(.x, na.rm = TRUE), 
                    sd = ~sd(.x, na.rm = TRUE)), 
               .names = "{.col}_{.fn}")
      ) %>%
      pivot_longer(everything(), names_to = "measure", values_to = "value") %>%
      separate(measure, into = c("metric", "time", "stat"), sep = "_") %>%
      pivot_wider(names_from = c(time, stat), values_from = value) %>%
      mutate(
        pre_mean_sd = paste0(round(pre_mean, 2), " (", round(pre_sd, 2), ")"),
        post_mean_sd = paste0(round(post_mean, 2), " (", round(post_sd, 2), ")")
      ) %>%
      select(metric, pre_mean_sd, post_mean_sd)
    
  } else {  # Ultimatum Game
    summary_stats <- data %>%
      summarise(
        across(c(mean_rt_pre, mean_rt_post, accept_rate_pre, accept_rate_post), 
               list(mean = ~mean(.x, na.rm = TRUE), 
                    sd = ~sd(.x, na.rm = TRUE)), 
               .names = "{.col}_{.fn}")
      ) %>%
      pivot_longer(everything(), names_to = "measure", values_to = "value") %>%
      separate(measure, into = c("metric", "time", "stat"), sep = "_") %>%
      pivot_wider(names_from = c(time, stat), values_from = value) %>%
      mutate(
        pre_mean_sd = paste0(round(pre_mean, 2), " (", round(pre_sd, 2), ")"),
        post_mean_sd = paste0(round(post_mean, 2), " (", round(post_sd, 2), ")")
      ) %>%
      select(metric, pre_mean_sd, post_mean_sd)
  }
  
  return(summary_stats)
}

# Generate summary tables
rl_summary_table <- create_summary_table(rl_behavioral_summary, "Reversal Learning")
ug_summary_table <- create_summary_table(ug_behavioral_summary, "Ultimatum Game")

# Statistical tests for behavioral changes
rl_rt_test <- t.test(rl_behavioral_summary$mean_rt_post, 
                     rl_behavioral_summary$mean_rt_pre, paired = TRUE)
rl_error_test <- t.test(rl_behavioral_summary$error_rate_post, 
                        rl_behavioral_summary$error_rate_pre, paired = TRUE)
rl_optimal_test <- t.test(rl_behavioral_summary$optimal_choice_rate_post, 
                          rl_behavioral_summary$optimal_choice_rate_pre, paired = TRUE)

ug_rt_test <- t.test(ug_behavioral_summary$mean_rt_post, 
                     ug_behavioral_summary$mean_rt_pre, paired = TRUE)
ug_accept_test <- t.test(ug_behavioral_summary$accept_rate_post, 
                         ug_behavioral_summary$accept_rate_pre, paired = TRUE)

# Create results summary
behavioral_test_results <- data.frame(
  Task = c("Reversal Learning", "Reversal Learning", "Reversal Learning",
           "Ultimatum Game", "Ultimatum Game"),
  Measure = c("Response Time", "Error Rate", "Optimal Choice Rate",
              "Response Time", "Accept Rate"),
  t_statistic = c(rl_rt_test$statistic, rl_error_test$statistic, rl_optimal_test$statistic,
                  ug_rt_test$statistic, ug_accept_test$statistic),
  p_value = c(rl_rt_test$p.value, rl_error_test$p.value, rl_optimal_test$p.value,
              ug_rt_test$p.value, ug_accept_test$p.value),
  effect_direction = c(
    ifelse(rl_rt_test$statistic > 0, "Increase", "Decrease"),
    ifelse(rl_error_test$statistic > 0, "Increase", "Decrease"),
    ifelse(rl_optimal_test$statistic > 0, "Increase", "Decrease"),
    ifelse(ug_rt_test$statistic > 0, "Increase", "Decrease"),
    ifelse(ug_accept_test$statistic > 0, "Increase", "Decrease")
  )
) %>%
  mutate(
    p_value_formatted = case_when(
      p_value < 0.001 ~ "< 0.001",
      p_value < 0.01 ~ paste0("= ", round(p_value, 3)),
      TRUE ~ paste0("= ", round(p_value, 2))
    ),
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# Print summary tables
cat("\n=== REVERSAL LEARNING BEHAVIORAL SUMMARY ===\n")
print(kable(rl_summary_table, 
            col.names = c("Measure", "Pre-Stimulation", "Post-Stimulation"),
            caption = "Reversal Learning Task: Behavioral Measures (Mean ± SD)"))

cat("\n=== ULTIMATUM GAME BEHAVIORAL SUMMARY ===\n")
print(kable(ug_summary_table,
            col.names = c("Measure", "Pre-Stimulation", "Post-Stimulation"),
            caption = "Ultimatum Game: Behavioral Measures (Mean ± SD)"))

cat("\n=== STATISTICAL TESTS FOR BEHAVIORAL CHANGES ===\n")
print(kable(behavioral_test_results %>% 
            select(Task, Measure, t_statistic, p_value_formatted, effect_direction, significance),
            col.names = c("Task", "Measure", "t", "p", "Direction", "Sig"),
            caption = "Paired t-tests for Pre vs Post-Stimulation Behavioral Changes"))

# Simple correlation with depression severity (if desired)
if("HDRS" %in% names(rl_behavioral_summary)) {
  cat("\n=== CORRELATIONS WITH DEPRESSION SEVERITY ===\n")
  
  rl_hdrs_cors <- rl_behavioral_summary %>%
    summarise(
      rt_hdrs_cor = cor(rt_change, HDRS, use = "complete.obs"),
      error_hdrs_cor = cor(error_rate_change, HDRS, use = "complete.obs"),
      optimal_hdrs_cor = cor(optimal_rate_change, HDRS, use = "complete.obs")
    )
  
  ug_hdrs_cors <- ug_behavioral_summary %>%
    summarise(
      rt_hdrs_cor = cor(rt_change, HDRS, use = "complete.obs"),
      accept_hdrs_cor = cor(accept_rate_change, HDRS, use = "complete.obs")
    )
  
  cat("Reversal Learning - HDRS correlations:\n")
  cat(sprintf("RT change: r = %.3f\n", rl_hdrs_cors$rt_hdrs_cor))
  cat(sprintf("Error rate change: r = %.3f\n", rl_hdrs_cors$error_hdrs_cor))
  cat(sprintf("Optimal choice change: r = %.3f\n", rl_hdrs_cors$optimal_hdrs_cor))
  
  cat("\nUltimatum Game - HDRS correlations:\n")
  cat(sprintf("RT change: r = %.3f\n", ug_hdrs_cors$rt_hdrs_cor))
  cat(sprintf("Accept rate change: r = %.3f\n", ug_hdrs_cors$accept_hdrs_cor))
}

# Save results for manuscript
write.csv(behavioral_test_results, "behavioral_summary_statistics.csv", row.names = FALSE)                     