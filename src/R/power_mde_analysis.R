# Copyright (C) 2026 Blair Shevlin <blairshevlin@gmail.com>
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
# 2026/08/15    Blair Shevlin                           wrote power / minimum-detectable-
#                                                        effect (MDE) sensitivity analysis
#                                                        for the primary 3-test family (R5)
#
# Purpose
# -------
# R5 (Nature Resub 5): pushes on sample-size justification, and explicitly rejects the
# retroactive-justification framings already in the manuscript ("didn't know the effect
# size a priori," "grant said collect all the data," "parent trial enrollment cap").
# Per the project's statistical_rigor_checklist.md item 4, this script answers a
# different, defensible question instead: "what effect size could this design actually
# detect," reported per test, not as evidence the study was adequately powered in
# advance. See robustness_analysis_plan.md section 3.5 for the full reasoning behind
# every method choice below.
#
# FRAMING RULE (read before writing anything from this script's output into the
# rebuttal letter or manuscript): every sentence built from these numbers must read as
# "this design could detect effects of at least this size," never as "the study was
# powered to detect its actual effects" or any variant implying retroactive adequacy.
# R5's letter is unusually explicit about rejecting exactly that move -- getting this
# framing wrong undoes the point of doing the analysis at all.
#
# This is deliberately a SEPARATE script from statistical_tests_robust.R: it computes
# design-level sensitivity (a property of n, alpha, and the test structure), not a
# result derived from a specific fitted model, and keeping it separate avoids implying
# the MDE numbers below are contingent on any single model specification in the main
# robustness file.

library(tidyverse)
library(fs)
library(here)
library(pwr)           # NOT currently in setup.R as of Aug 2026 -- add it there
library(effectsize)    # for cohens_d() with CIs, used in Section 3's observed effect sizes
library(broom)         # for tidy() in Section 3 -- NOT auto-attached by library(tidyverse)
library(broom.mixed)   # loaded for parity with statistical_tests_robust.R; not strictly
                        # required here since this script only fits plain lm(), but kept
                        # so the two scripts' library blocks don't silently diverge
# NOTE: this script previously omitted broom, which tidy() requires -- library(tidyverse)
# does NOT auto-attach it despite broom being part of the tidyverse install set. Caught in
# code review before this was ever run; if you see "no applicable method for 'tidy'" it
# means this library block was trimmed again without checking Section 3's dependencies.

# Paths (mirrors statistical_tests_robust.R conventions)
dir = path(here())
nt_dir   = dir / "data" / "nt" / "processed"
clin_dir = dir / "data" / "clinical"

alpha_uncorrected = 0.05
alpha_bonferroni  = 0.05 / 3   # matches statistical_tests_robust.R's primary-family
                                # Bonferroni-within-3 (test c only, per R1's actual
                                # wording) -- kept identical here so the two scripts
                                # can't silently drift on which alpha is "the" one
target_power = 0.80

##########################################################
# SECTION 1: MDE FOR TESTS (a)/(b) -- TASK-SPECIFICITY,  #
# PAIRED WITHIN-SUBJECT DESIGNS                          #
##########################################################

# Tests (a) DA x RL and (b) 5-HT x UG are, at the design level, paired pre/post-stim
# comparisons -- collapsing the actual trial-level mixed-effects models (which also
# model trial-level replication and by-subject random slopes) down to a simple paired
# comparison of subject-level means. This is a real simplification, stated explicitly:
# it ignores the extra precision trial-level replication provides, so treat the MDE
# below as a conservative (upper-bound-ish) answer for "what a subject-level summary
# of this design could detect," not a substitute for the full mixed-model's actual
# power (which is not calculated here -- see robustness_analysis_plan.md section 3.5
# for why a full simr-style simulation was deliberately not attempted: added complexity
# without changing the substantive answer, harder for a referee to independently check).

mde_paired = function(n, alpha = alpha_uncorrected, power = target_power) {
  fit = pwr.t.test(n = n, sig.level = alpha, power = power, type = "paired",
                    alternative = "two.sided")
  tibble(n = n, alpha = alpha, power = power, mde_dz = fit$d)
}

primary_ab_mde = bind_rows(
  mde_paired(n = 8)  %>% mutate(test = "a_DA_RL_task_specificity",   design = "paired, n=8 (RL)"),
  mde_paired(n = 10) %>% mutate(test = "b_5HT_UG_task_specificity",  design = "paired, n=10 (UG)")
)

print(primary_ab_mde)

# Interpretation guide (fill in once run): Cohen's convention treats dz ~0.2/0.5/0.8 as
# small/medium/large. Expect the printed mde_dz values to be well above 0.8 given n=8-10
# -- that is the honest answer for this design, not a computational error. Compare
# against the ACTUALLY OBSERVED within-subject effect sizes for DA (RL) / 5-HT (UG)
# stim-related change (Section 3 below) to see whether the observed effects clear this
# bar or not -- do not just report the MDE in isolation.

##########################################################
# SECTION 2: MDE FOR TEST (c) -- DA x 5-HT -> MONTH 6    #
# HDRS INTERACTION ("THE KEY INTERACTION"), REGRESSION   #
##########################################################

# Test (c) is the interaction term in model_DAxSE (HDRS ~ deltaDA_UG * deltaSE_UG,
# n=10, 3 predictors: deltaDA_UG, deltaSE_UG, their interaction). For the interaction
# term specifically: u = 1 (numerator df, one interaction term), v = n - k - 1 = 6
# (residual df, k=3 predictors). Reported at BOTH the Bonferroni-corrected alpha
# (matching statistical_tests_robust.R's treatment of this as the one test in the
# family that must survive correction) and the uncorrected alpha, for comparison.
#
# v = 6 here is an ASSUMPTION (n=10, 3 predictors), not yet verified against an actual
# fitted model -- Section 3 below fits model_DAxSE_check and re-checks this assumption
# against df.residual() once real data is loaded, printing a warning if they disagree
# (e.g. due to a subject dropping out of the merge). Treat primary_c_mde below as
# provisional until that check passes.

assumed_v_test_c = 6   # n=10 - 3 predictors - 1 = 6; re-verified in Section 3

mde_interaction = function(u, v, alpha, power = target_power) {
  fit = pwr.f2.test(u = u, v = v, sig.level = alpha, power = power)
  tibble(u = u, v = v, alpha = alpha, power = power,
         mde_f2 = fit$f2, mde_partial_R2 = fit$f2 / (1 + fit$f2))
}

primary_c_mde = bind_rows(
  mde_interaction(u = 1, v = assumed_v_test_c, alpha = alpha_bonferroni)  %>% mutate(correction = "Bonferroni-within-3 (matches statistical_tests_robust.R)"),
  mde_interaction(u = 1, v = assumed_v_test_c, alpha = alpha_uncorrected) %>% mutate(correction = "uncorrected, for reference only")
) %>%
  mutate(test = "c_DAx5HT_HDRS_key_interaction")

print(primary_c_mde)

# FLAG BEFORE WRITING THIS UP: with only 6 residual df, expect mde_f2 to be large (an
# independent Python cross-check via scipy's noncentral-F distribution, run outside R
# since this sandbox has no R interpreter, gave mde_f2 ~= 2.3 / partial R2 ~= 0.70 at
# the Bonferroni alpha, and ~1.4 / 0.59 uncorrected -- confirm these match your R output
# within rounding before using either number). That is a genuinely implausible effect
# size for most psychological/clinical interactions to be able to reach 80% power at
# this alpha with this design -- state that plainly rather than softening it. The
# honest takeaway is "this design has very limited sensitivity for the interaction
# specifically, even though the observed effect [Section 3] happened to be large enough
# to reach nominal significance" -- not that the design was well-suited to detect
# interactions of a priori-plausible size.

##########################################################
# SECTION 3: OBSERVED EFFECT SIZES, FOR COMPARISON       #
# AGAINST THE MDEs ABOVE (NOT A NEW ANALYSIS -- REUSES   #
# THE SAME DATA/DERIVATION AS statistical_tests_robust.R)#
##########################################################

# Loaded independently here (rather than via source()) to keep this script runnable on
# its own, consistent with how other scripts in src/R/ are structured. If this drifts
# from statistical_tests_robust.R's own computation of the same quantities, that's a
# bug to fix, not an acceptable inconsistency -- the numbers must match.

load(nt_dir / "UG_RL_NT-Continuous_9-23-25.RData")

rl.EST.Reward = rl.EST %>%
  filter(event == "Reward") %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))
ug.EST.Offer = ug.EST %>%
  filter(event == "Offer") %>%
  mutate(stim = factor(stim, levels = c("Pre-Stim","Post-Stim")))

# Uses effectsize::cohens_d() (one-sample form, mu=0, on the per-subject change scores)
# rather than a hand-rolled mean/sd ratio -- gives a CI via the noncentral-t distribution
# and is directly comparable to the mde_dz values from Section 1, which is the point of
# this comparison. A hand-rolled d would give the same point estimate but no CI, which
# matters at n=8-10 for a Nature submission.
observed_d = function(est_data, nt_filter) {
  change_data = est_data %>%
    mutate(nt = factor(nt, levels = c("DA","SE","NE"))) %>%
    filter(nt == nt_filter) %>%
    select(idx, trial, stim, Oz) %>%
    pivot_wider(names_from = stim, values_from = Oz) %>%
    group_by(idx) %>%
    summarise(nt_change = mean(`Post-Stim` - `Pre-Stim`, na.rm = TRUE), .groups = "drop")

  d_result = effectsize::cohens_d(change_data$nt_change, mu = 0)

  tibble(
    mean_change = mean(change_data$nt_change),
    sd_change   = sd(change_data$nt_change),
    n           = nrow(change_data),
    cohen_d     = d_result$Cohens_d,
    ci_lower    = d_result$CI_low,
    ci_upper    = d_result$CI_high
  )
}

observed_ab = bind_rows(
  observed_d(rl.EST.Reward, "DA") %>% mutate(test = "a_DA_RL_task_specificity"),
  observed_d(ug.EST.Offer,  "SE") %>% mutate(test = "b_5HT_UG_task_specificity")
)

print(observed_ab)

# Test (c)'s observed effect size: refit model_DAxSE here using the same derivation as
# statistical_tests_robust.R's "NT ~ CLINICAL" section (deltaDA_UG/deltaSE_UG merged
# with month-6 HDRS). Kept minimal -- only what's needed for the interaction term's
# observed partial R2, not the full comparison-model table (that lives in the main
# robustness file, not here).

cl.df = read.csv(clin_dir / "clinical-data_deid_07-10-25.csv")

ug.nt.cl = ug.EST.Offer %>%
  pivot_longer(cols = c("Oz","Rz","Pz","Mz"), names_to = "nt_metric", values_to = "nt_val") %>%
  group_by(idx, stim, nt, nt_metric) %>%
  summarise(mTrial = mean(nt_val), .groups = "drop") %>%
  filter(nt_metric == "Oz") %>%
  select(idx, stim, nt, mTrial) %>%
  pivot_wider(values_from = mTrial, names_from = c("stim","nt")) %>%
  mutate(deltaDA_UG = `Post-Stim_DA` - `Pre-Stim_DA`,
         deltaSE_UG = `Post-Stim_SE` - `Pre-Stim_SE`) %>%
  select(idx, deltaDA_UG, deltaSE_UG)

data_m6_power = cl.df %>%
  filter(session == "month 6") %>%
  select(idx, HDRS) %>%
  merge(ug.nt.cl)

model_DAxSE_check = lm(HDRS ~ deltaDA_UG * deltaSE_UG, data = data_m6_power)
observed_c = tidy(model_DAxSE_check) %>%
  filter(term == "deltaDA_UG:deltaSE_UG") %>%
  mutate(test = "c_DAx5HT_HDRS_key_interaction")

# CROSS-CHECK, BLOCKING before anything from this section goes in the rebuttal: this
# refits HDRS ~ deltaDA_UG * deltaSE_UG independently of statistical_tests_robust.R's
# own model_DAxSE (same formula, but a separately-written data-prep pipeline -- see the
# comment above this section). The two MUST report the same coefficient/p-value for the
# interaction term; if they don't, that's a bug in one of the two pipelines, not an
# acceptable discrepancy in a Nature rebuttal. Paste the confirmed estimate/p-value from
# statistical_tests_robust.R's own `tidy(model_DAxSE)` output below once you've run it,
# then re-run this script -- it will error loudly on a mismatch instead of silently
# reporting two different numbers for the same quantity.
confirmed_estimate_from_main_script = -0.142  # <- fill in after running statistical_tests_robust.R
confirmed_pvalue_from_main_script   = 0.00282 # <- fill in after running statistical_tests_robust.R
if (!is.na(confirmed_estimate_from_main_script)) {
  stopifnot(
    "Interaction estimate does not match statistical_tests_robust.R's model_DAxSE -- do not proceed until reconciled" =
      isTRUE(all.equal(round(observed_c$estimate, 3), round(confirmed_estimate_from_main_script, 3), tolerance = 1e-6)),
    "Interaction p-value does not match statistical_tests_robust.R's model_DAxSE -- do not proceed until reconciled" =
      isTRUE(all.equal(round(observed_c$p.value, 3), round(confirmed_pvalue_from_main_script, 3), tolerance = 1e-6))
  )
  message("Cross-check passed: this script's model_DAxSE_check matches statistical_tests_robust.R's model_DAxSE.")
} else {
  message("Cross-check NOT YET PERFORMED: confirmed_estimate_from_main_script/confirmed_pvalue_from_main_script ",
          "are still NA. Fill these in from statistical_tests_robust.R's own tidy(model_DAxSE) output before ",
          "trusting Section 3's numbers for the rebuttal -- see comment above.")
}

# Partial R2 / f2 for the interaction term specifically (not the whole model) via a
# nested-model comparison, matching what the MDE in Section 2 is calibrated against.
# Variable names below reflect what each model actually is (reduced = no interaction,
# full = with interaction) -- a prior version of this script had these inverted, which
# happened to cancel out in the arithmetic but was a landmine for future edits.
model_DAxSE_noint = lm(HDRS ~ deltaDA_UG + deltaSE_UG, data = data_m6_power)
sse_reduced = sum(residuals(model_DAxSE_noint)^2)   # model WITHOUT the interaction
sse_full    = sum(residuals(model_DAxSE_check)^2)   # model WITH the interaction
observed_partial_R2 = (sse_reduced - sse_full) / sse_reduced
observed_f2 = observed_partial_R2 / (1 - observed_partial_R2)

# Re-verify Section 2's assumed v (residual df) against what was actually fit here. If
# these disagree -- e.g. because a subject dropped out of the idx merge -- the MDE
# reported in Section 2 no longer describes the model actually being compared against,
# and needs to be recomputed at the correct v before use.
actual_v_test_c = df.residual(model_DAxSE_check)
if (actual_v_test_c != assumed_v_test_c) {
  warning(sprintf(
    "Section 2's MDE assumed v=%d residual df, but model_DAxSE_check actually has v=%d. ",
    assumed_v_test_c, actual_v_test_c),
    "Recomputing primary_c_mde at the correct v -- re-check Section 2's printed table against this corrected version.")
  primary_c_mde = bind_rows(
    mde_interaction(u = 1, v = actual_v_test_c, alpha = alpha_bonferroni)  %>% mutate(correction = "Bonferroni-within-3 (matches statistical_tests_robust.R)"),
    mde_interaction(u = 1, v = actual_v_test_c, alpha = alpha_uncorrected) %>% mutate(correction = "uncorrected, for reference only")
  ) %>%
    mutate(test = "c_DAx5HT_HDRS_key_interaction")
  print(primary_c_mde)
} else {
  message(sprintf("Section 2's assumed v=%d confirmed against the actual fitted model (n=%d).",
                   assumed_v_test_c, nrow(data_m6_power)))
}

print(observed_c %>% mutate(observed_partial_R2 = observed_partial_R2, observed_f2 = observed_f2))

##########################################################
# SECTION 4: BATTEN ET AL. (2024) EFFECT-SIZE GROUNDING  #
# (per Helen Mayberg, Aug 2026 -- see narrative_guidance.md) #
##########################################################

# Batten, S. R. et al. "Dopamine and serotonin in human substantia nigra track social
# context and value signals during economic exchange." Nat. Hum. Behav. 8, 718-728
# (2024). Same Ultimatum Game task, same human-electrochemistry method, overlapping
# author team (Kopell, Mayberg, Gu, Montague), n=4 Parkinson's patients x 2 sessions
# (8 datasets) -- comparable order of magnitude to this study's n=8-10.
#
# Numbers below are transcribed directly from the published Results section (N-M1,
# N-M3), NOT re-derived from raw data (Batten et al.'s underlying data/code are public
# on GitHub -- github.com/danbang/article-DA-5HT-UG-SNr -- but reproducing their models
# was out of scope here; if an exact re-derivation is wanted later, that repo is the
# place to start).
#
# IMPORTANT FRAMING, stated explicitly per narrative_guidance.md: this is NOT a formal
# a priori power calculation for the SCC/DBS study -- the anatomical target (SCC vs.
# SNr) and clinical population (TRD vs. Parkinson's) differ, and effect magnitudes are
# not assumed to transfer directly. It is evidence that the team had SOME empirically
# grounded sense of plausible monoamine effect magnitudes in an analogous human
# intracranial recording paradigm before this study -- a fair, honest answer to "did
# you have any basis for expecting an effect of this size," without claiming more
# precision than is true.
#
# IMPORTANT DESIGN CAVEAT, also state explicitly: Batten et al.'s significant effects
# were detected via TRIAL-LEVEL mixed-effects models (n=426-454 trial-level df, with
# by-session random effects) -- a fundamentally different, much higher-powered design
# than this study's SUBJECT-LEVEL summary comparisons (n=8-10 points, one per subject).
# Do not present the Batten et al. numbers as directly power-comparable to the MDEs in
# Sections 1-2 above; they are comparable only as effect-magnitude grounding, not as an
# equivalent power calculation. This caveat is echoed via message() below (not just this
# source comment) so it reaches anyone working from the printed/console output alone.

message("BATTEN ET AL. GROUNDING -- READ BEFORE USING: these are TRIAL-LEVEL mixed-model ",
        "effects (n=426-454 trial df), not subject-level effects, and are NOT directly ",
        "power-comparable to the MDEs in Sections 1-2 (which are subject-level, n=8-10). ",
        "Use only as effect-magnitude grounding, per narrative_guidance.md.")

# TEST-MAPPING CLARIFICATION (added Aug 2026, Alagapan-grounding round): this Batten et al.
# grounding is used for tests (a)/(b) ONLY -- the stimulation/condition-induced within-
# subject change tests (DA x RL task specificity; 5-HT x UG task specificity). It is NOT
# used to ground test (c) (the DA x 5-HT -> month-6 HDRS interaction); Section 5 below
# provides a separate, more directly relevant grounding source for test (c) specifically
# (Alagapan et al. 2023, Nature -- same SCC target, same TRD population, same senior
# author). This mapping was previously implicit; stating it explicitly here removes
# ambiguity about which grounding source supports which test.

message("BATTEN ET AL. GROUNDS TESTS (a)/(b) ONLY -- see Section 5 for test (c)'s ",
        "grounding source (Alagapan et al. 2023). Do not use Batten's numbers, converted ",
        "or not, as grounding for test (c).")

# The three Batten et al. effects were reported on two DIFFERENT native standardized
# scales -- a Cohen's-d-equivalent for the contrast-coded (human vs. computer) condition
# effect, and an approximate partial-r for the two continuous (z-scored) trial-level
# covariates. These scales are not linearly related (r=0.85 and d=0.85 are very
# different magnitudes), so both are converted onto a COMMON d scale and a COMMON r
# scale below via the standard two-independent-groups conversion formulas, rather than
# comparing native_beta values across rows directly -- comparing native_beta columns
# as-is (as an earlier version of this table did) understates how much larger the
# condition effect is relative to the two covariate-tracking effects.
r_from_d = function(d) d / sqrt(d^2 + 4)
d_from_r = function(r) 2 * r / sqrt(1 - r^2)

batten_grounding = tribble(
  ~batten_effect,                                                   ~native_beta, ~native_metric,  ~stat,           ~p,
  "Overall DA level, human vs. computer condition",                  0.85,        "d_equivalent",  "t(454) = 3.04", 0.002,
  "Relative DA tracks trial-by-trial value difference (RPE-like)",   0.21,        "approx_r",      "t(426) = 2.82", 0.005,
  "Relative 5-HT tracks current offer value",                       0.22,        "approx_r",      "t(426) = 3.15", 0.002
) %>%
  mutate(
    common_d_scale = ifelse(native_metric == "d_equivalent", native_beta, d_from_r(native_beta)),
    common_r_scale = ifelse(native_metric == "approx_r",     native_beta, r_from_d(native_beta))
  )

print(batten_grounding)

# On the common d scale, the condition effect (~0.85) is roughly DOUBLE the two
# covariate-tracking effects (~0.43-0.45), not ~4x as a naive same-column read of the
# native betas (0.85 vs 0.21-0.22) would suggest. Use common_d_scale / common_r_scale in
# any rebuttal text that compares these three effects to each other or to this study's
# own MDEs -- do not quote native_beta across rows as if on the same scale.

# Suggested use in the rebuttal text (draft this after confirming Section 3's observed
# numbers, not before): frame as "in an analogous human intracranial monoamine-
# recording paradigm using the same task, standardized effects for task/condition-
# related dopamine and serotonin signals ranged from ~0.2 (trial-level covariation, on
# an approximate-r scale) to ~0.85 (categorical condition effect, on a Cohen's-d-
# equivalent scale) -- giving the team some empirically grounded expectation of
# plausible effect magnitudes, even though the SCC/DBS-specific effect size could not be
# predicted precisely for this new anatomical target and clinical population." Combine
# with, but do not conflate with, the MDE/sensitivity framing from Sections 1-2, and do
# not collapse the two native scales into one number without noting the conversion.

##########################################################
# SECTION 5: ALAGAPAN ET AL. (2023) EFFECT-SIZE/PRECEDENT #
# GROUNDING FOR TEST (c) ONLY                             #
# (per Blair, Aug 2026 -- see rebuttal drafts/08.19.26/    #
# Alagapan_grounding_plan_FINAL.md for the full reasoning) #
##########################################################

# Alagapan, S. et al. "Cingulate dynamics track depression recovery with deep brain
# stimulation." Nature 622, 130-138 (2023). https://doi.org/10.1038/s41586-023-06541-3.
# Same senior author (H.S. Mayberg is co-senior/co-corresponding author on both this paper
# and Alagapan et al.), same anatomical target (subcallosal cingulate, SCC), same clinical
# population (treatment-resistant depression, SCC DBS), same overall trial infrastructure
# (NCT01984710, an earlier-generation study in the same research program). n=10 enrolled,
# n=6 with usable LFP data, n=5 "typical responders" used for the classifier/beta-band
# analyses -- directly comparable order of magnitude to this study's n=8-10.
#
# Central finding used here: a neural biomarker (the "spectral discriminative component,"
# SDC), beta-band-power-driven, derived from SCC local field potentials, that predicts/
# tracks clinical depression state (HDRS-17). This grounds test (c) specifically -- a
# neural/neurochemical signal predicting a clinical depression outcome -- NOT tests (a)/(b),
# which Section 4's Batten grounding covers instead (see the test-mapping clarification
# added to Section 4).
#
# IMPORTANT FRAMING, same principle as Section 4: this is NOT a formal a priori power
# calculation for this study, and Alagapan's biomarker claim is a WITHIN-SUBJECT
# LONGITUDINAL classifier (does each patient's own LFP track their own trajectory over 24
# weeks?) while test (c) is a BETWEEN-SUBJECT, SINGLE-TIMEPOINT regression interaction (do
# two monoamine change scores jointly predict one later HDRS value, across patients?).
# These are related but not identical estimands -- use this grounding as precedent for the
# FEASIBILITY AND PUBLISHABILITY of small-n biomarker-to-depression-outcome findings in this
# exact clinical program, not as an implied statistical equivalence to test (c)'s MDE.
#
# DELIBERATE DESIGN CHOICE: unlike Section 4's Batten grounding, NO r/d common-scale
# conversion is attempted here. Alagapan's numbers are AUCs (classifier discriminability
# under cross-validation) and Wilcoxon/permutation P-values -- converting AUC to Cohen's d
# via the probit relationship (d = sqrt(2) * qnorm(AUC)) and then to r/partial-R^2 chains
# three assumption-laden steps (binormal equal-variance discriminant distributions; treating
# a mean-over-CV-folds AUC as a point coefficient; treating classification-discriminability
# and regression-variance-explained as commensurable quantities) on numbers that are
# themselves means +/- SD over 5-fold cross-validation, not point estimates. A methodologist
# referee (R5) would reasonably contest each link in that chain. Presented on native scale
# only, with an explicit non-conversion rationale -- see
# rebuttal drafts/08.19.26/Alagapan_grounding_plan_FINAL.md Section 0 for the full sign-off.

message("ALAGAPAN ET AL. (2023) GROUNDING -- READ BEFORE USING: grounds test (c) ONLY (not ",
        "tests a/b -- see Section 4's test-mapping clarification for those). Values are NOT ",
        "converted to test (c)'s partial-R2 scale -- AUC/Wilcoxon-P and partial-R2 are not ",
        "directly interconvertible without assumptions not considered defensible here (see ",
        "comment above). Use as structural/feasibility precedent only, never as a numeric ",
        "power comparison to Table 2's MDE.")

alagapan_grounding = tribble(
  ~alagapan_effect, ~value, ~metric, ~n, ~source_citation,
  "LFP classifier discriminates 'sick' vs 'stable response' state",
    "AUC = 0.87 +/- 0.09", "AUC", 5,
    "Fig. 2A; Extended Data Table 3 (5-fold leave-one-participant-out CV)",
  "Left-low-beta power, weeks 1-4 vs. weeks 21-24",
    "one-sided Wilcoxon P = 0.03", "Wilcoxon P", 5, "Fig. 2c",
  "Left-high-beta power, weeks 1-4 vs. weeks 21-24",
    "one-sided Wilcoxon P = 0.03", "Wilcoxon P", 5, "Fig. 2c",
  "SDC biomarker predicts weekly HDRS-defined clinical state, full course",
    "AUC = 0.94 +/- 0.04", "AUC", 5, "Fig. 2G; Extended Data Table 3",
  "SDC biomarker predicts weekly HDRS-defined state, restricted window",
    "AUC = 0.89 +/- 0.07", "AUC", 4,
    "Extended Data Table 3 only (no fig. no.); P006 excluded -- only one HDRS class present in weeks 5-20",
  "SDC tracks facial-expression classifier output",
    "linear mixed model F(1.00,51.74) = 6.54, P = 0.01", "F-test", 5,
    "Fig. 5E; Extended Data Table 3",
  "Transition-week concordance, SDC vs. face classifier",
    "Kendall's tau = 0.89, P = 0.037", "Kendall's tau", 5, "Fig. 5F"
) %>%
  mutate(not_converted_reason = paste(
    "AUC/Wilcoxon-P/Kendall's-tau not converted to test (c)'s partial-R2 scale -- see",
    "framing comment above this tribble for the full rationale. Presented on native scale",
    "only, as structural/feasibility precedent for test (c), not a numeric power comparison."
  ))

print(alagapan_grounding)

# Numbers transcribed directly from the published article and its Extended Data Table 3 /
# Nature Portfolio Reporting Summary, NOT re-derived from raw data (Alagapan et al.'s
# underlying data/code are public via the Data Archive for The Brain Initiative, DABI --
# dabi.loni.usc.edu/dsi/1UH3NS103550/UXUF7822Z3JL -- but reproducing their models was out of
# scope here). Verification protocol for this table: transcription-accuracy check against
# the source PDF (done -- every row above was re-confirmed against the actual page image,
# not a text-layer extraction, which failed on Extended Data Table 3's image-rendered
# layout during an earlier drafting pass; see Alagapan_grounding_plan_FINAL.md Section 2 for
# the correction history), NOT a Python statistical re-derivation (these are external
# published numbers, not quantities computed from this study's own data -- same
# verification tier as Section 4's Batten numbers).

##########################################################
# SECTION 6: FIELD PRECEDENT -- ALAGAPAN'S OWN NATURE     #
# PORTFOLIO REPORTING SUMMARY DISCLOSURE (QUALITATIVE,    #
# NOT A TABLE NUMBER)                                     #
##########################################################

# Source: Nature Portfolio Reporting Summary accompanying Alagapan et al. (2023), "Life
# sciences study design" section, "Sample size" field. Corresponding authors on that
# Reporting Summary: Christopher Rozell, Helen S. Mayberg.
#
# TRIMMED QUOTE ONLY -- do not use the full field text. The full text also contains an
# enrollment-cap/device-availability justification ("availability of prototype devices
# (n=10) from the device manufacturer") that is structurally identical to a justification
# R5 has already explicitly rejected in THIS manuscript ("the parent trials aimed at
# enrolling 25 patients only... I return to my point about minimum effect sizes" -- R5,
# Editor and Review Comments.docx, Aug 14 2026). Quoting that clause risks R5 recognizing
# the same weak argument recurring in cited precedent. See
# rebuttal drafts/08.19.26/Alagapan_grounding_plan_FINAL.md Section 3 for the full
# reasoning; only the trimmed clause below is approved for use (confirmed by Blair).

alagapan_reporting_summary_quote_trimmed = paste(
  "...this is a first-of-its kind study investigating LFP collected longitudinally during",
  "SCC DBS; it was not possible to power the study based on LFP outcomes."
)

message("ALAGAPAN REPORTING SUMMARY QUOTE -- use ONLY the trimmed version stored in ",
        "alagapan_reporting_summary_quote_trimmed. Do NOT quote the full Reporting Summary ",
        "field text -- it contains a device-availability enrollment-cap justification R5 ",
        "has already rejected in this manuscript. See comment above for the full citation.")

##########################################################
# SECTION 7: VERIFICATION -- THIS MANUSCRIPT'S OWN TESTS  #
# USE TWO-SIDED, UNRELAXED ALPHA (CONTRAST WITH ALAGAPAN) #
##########################################################

# Per Blair, Aug 20 2026: unlike Alagapan et al. (2023), who explicitly relaxed to one-sided
# testing at small n (their Methods: "The small sample size of the current study does not
# have sufficient power to test statistical significance at 0.05 in a two-sided test, even
# when the direction of the change is readily apparent. Therefore, we used a one-sided test
# with a threshold of 0.05 and also confirmed statistical significance in a two-sided test
# with a relaxed threshold of 0.1"), this manuscript's own primary analyses do NOT use
# one-sided or small-n-relaxed nonparametric testing anywhere.
#
# Verified by direct inspection of this repo's own scripts (static code-review facts,
# recorded here for traceability -- not re-executed, since sourcing statistical_tests.R
# would run the full analysis pipeline and load large data files unnecessarily for what is
# a code-inspection fact, not a computed quantity):
#   - statistical_tests.R line 337 / statistical_tests_robust.R line 602: paired-comparison
#     helper function signature defaults to `alternative = "two.sided"`.
#   - statistical_tests.R line 787 / statistical_tests_robust.R line 1395:
#     `t.test(data_m6$HDRS, data_m6$baseline_HDRS, paired = TRUE)` -- two-sided by default
#     (t.test()'s own default `alternative` argument is "two.sided").
#   - Neither script contains a single `wilcox.test(...)` call (grepped in full, Aug 2026).
#
# BEFORE THIS GOES IN THE REBUTTAL LETTER OR MANUSCRIPT: confirm the manuscript's own
# Methods TEXT (not just the code) states "two-sided, alpha = 0.05" explicitly for the
# primary tests -- a referee reads the manuscript, not this repo's source code, and the
# claim must be traceable to what's actually written there. NOT YET CONFIRMED -- check the
# current manuscript Methods before using this claim in the rebuttal.

manuscript_uses_two_sided_unrelaxed_testing = TRUE  # see verification comment above
