
# =============================================================================
# MATCHED SUBSAMPLING FOR NEUROIMAGING STUDY
# VBM Sensitivity Analysis: 1:1 age- and sex-matched subsample (n=32 per group)
#
# Purpose:
#   Select a balanced subsample of n=32 healthy older adults (HO) matched
#   to the 32 MCI patients on age and sex. This matched subsample will be
#   used for the VBM/SBM sensitivity analysis (Supplementary Table S4).
#
# Three matching strategies are compared:
#   (A) Scaled Euclidean distance
#   (B) Standard Mahalanobis distance
#   (C) Robust Mahalanobis distance
#
# ----------------------------------------------------------------------------
# WHAT IS EACH DISTANCE METRIC?
# ----------------------------------------------------------------------------
#
# (A) Scaled Euclidean distance
#   The standard straight-line distance between two points in covariate space,
#   after each covariate has been standardized (z-scored: subtract mean, divide
#   by SD). Standardization is essential here because age (range ~56–88) and
#   sex (0/1) are on very different scales — without it, age would completely
#   dominate the distance and sex would be effectively ignored.
#   Limitation: treats all covariates as independent and equally variable,
#   which is fine here since age and sex are unlikely to be strongly correlated.
#
# (B) Standard Mahalanobis distance
#   An extension of Euclidean distance that additionally accounts for the
#   correlation structure between covariates, using the sample covariance
#   matrix. Two subjects who differ on a covariate that is highly variable
#   in the sample will be penalized less than two subjects who differ on a
#   covariate that is tightly distributed. Result: a scale-invariant,
#   correlation-adjusted distance.
#   Limitation: the sample covariance matrix can be distorted by outliers,
#   which matters here because the HO age range (56–88) is wider than MCI
#   (59–83) and includes some younger/older extremes.
#
# (C) Robust Mahalanobis distance
#   Same principle as (B), but the covariance matrix is estimated using the
#   Minimum Covariance Determinant (MCD) algorithm, which finds the subset
#   of observations with the most compact covariance structure and uses that
#   for estimation. This downweights extreme observations (outliers in age),
#   producing a more stable distance estimate in small samples.
#   Recommended when covariates have unequal spread across groups or mild
#   outliers — both of which apply here.
#
# ----------------------------------------------------------------------------
# BALANCE METRIC — Standardized Mean Difference (SMD)
# ----------------------------------------------------------------------------
#   SMD = (mean_treated - mean_control) / pooled_SD
#   This is equivalent to Cohen's d. It measures the magnitude of the
#   remaining imbalance on each covariate after matching, independent of
#   sample size (unlike p-values, which are sensitive to N).
#
#   Thresholds (Rubin 2001; Austin 2011):
#     SMD < 0.10  → excellent balance  ✅
#     SMD < 0.25  → acceptable balance ⚠
#     SMD ≥ 0.25  → poor balance       ❌ (re-match or adjust caliper)
#
# ----------------------------------------------------------------------------
# PACKAGE: MatchIt (Ho, Imai, King & Stuart, J. Stat. Softw., 2011)
# ----------------------------------------------------------------------------
#   MatchIt implements non-parametric preprocessing for causal inference.
#   It reduces covariate imbalance between groups before statistical analysis,
#   without imposing a parametric model on the outcome variable.
#   Reference: doi:10.18637/jss.v042.i08
#
# Author: Konstantin Toussas
# Date:   May 2026
# =============================================================================


# ── 0. Load packages ──────────────────────────────────────────────────────────

library(MatchIt)   # Matching algorithms and balance diagnostics
library(dplyr)     # Data manipulation
library(readxl)


# ── 1. Load and prepare data ──────────────────────────────────────────────────

df <- read_excel(file.choose())   # demographics table: see data_templates/demographics_template.csv
out_dir <- dirname(file.choose())  # choose any file in the folder where the matched-subsample CSVs should be written
# Columns: CODE (subject ID), AGE (years), SEX (F/M), GROUP (HO/MCI)

# Encode GROUP as binary: MCI = 1 (patients), HO = 0 (control)
# MatchIt convention: the group of interest (here MCI) is coded as 1
df$group_bin <- ifelse(df$GROUP == "MCI", 1, 0)

# Encode SEX as a factor (required for exact matching in MatchIt)
df$SEX <- as.factor(df$SEX)

cat("=== Full sample demographics (before matching) ===\n")
df %>%
  group_by(GROUP) %>%
  summarise(
    n          = n(),
    age_mean   = round(mean(AGE), 2),
    age_sd     = round(sd(AGE), 2),
    n_female   = sum(SEX == "F"),
    n_male     = sum(SEX == "M")
  ) %>%
  print()


# ── 2. Helper function: run matching + extract SMD ────────────────────────────

run_matching <- function(distance_method, label, seed = 42) {
  set.seed(seed)
  m <- matchit(
    group_bin ~ AGE,         # Match on age; sex handled via exact =
    data     = df,
    method   = "nearest",    # Greedy 1:1 nearest-neighbor
    distance = distance_method,
    exact    = ~ SEX,        # Hard constraint: match within same sex only
    ratio    = 1,            # One HO per MCI patient
    replace  = FALSE         # No reuse of control participants
  )
  s <- summary(m, standardize = TRUE)

  # Extract post-matching SMD table
  smd <- as.data.frame(s$sum.matched)
  smd$Covariate <- rownames(smd)
  smd$Method    <- label
  smd <- smd[, c("Method", "Covariate",
                 "Std. Mean Diff.", "Var. Ratio",
                 "eCDF Mean", "eCDF Max")]

  list(matchit_obj = m, summary_obj = s, smd = smd, label = label)
}


# ── 3. Run all three matching strategies ──────────────────────────────────────

results_A <- run_matching("scaled_euclidean", "A: Scaled Euclidean")
results_B <- run_matching("mahalanobis",       "B: Standard Mahalanobis")
results_C <- run_matching("robust_mahalanobis","C: Robust Mahalanobis")


# ── 4. Print individual summaries ─────────────────────────────────────────────

for (res in list(results_A, results_B, results_C)) {
  cat(paste0("\n=== ", res$label, " ===\n"))
  print(res$summary_obj)
}


# ── 5. Side-by-side SMD comparison across all three methods ───────────────────

cat("\n")
cat("================================================================\n")
cat(" SMD COMPARISON — All three distance methods (post-matching)\n")
cat("================================================================\n")
cat(" Target: SMD < 0.10 (excellent) | SMD < 0.25 (acceptable)\n")
cat("----------------------------------------------------------------\n")

all_smd <- rbind(results_A$smd, results_B$smd, results_C$smd)
print(all_smd, row.names = FALSE)

# Flag any poor balance
poor <- all_smd[abs(all_smd$`Std. Mean Diff.`) > 0.25, ]
if (nrow(poor) > 0) {
  cat("\n⚠️  WARNING: Poor balance (SMD > 0.25) detected:\n")
  print(poor[, c("Method", "Covariate", "Std. Mean Diff.")], row.names = FALSE)
} else {
  cat("\n✅  All covariates: SMD < 0.25 for all three methods.\n")
}


# ── 6. Formal balance tests for each matched sample ───────────────────────────

for (res in list(results_A, results_B, results_C)) {
  matched <- match.data(res$matchit_obj)
  cat(paste0("\n--- Formal balance tests: ", res$label, " ---\n"))
  cat("Age t-test:  ")
  tt <- t.test(AGE ~ group_bin, data = matched)
  cat(sprintf("t = %.3f, df = %.1f, p = %.4f\n", tt$statistic, tt$parameter, tt$p.value))
  cat("Sex chi-sq: ")
  cs <- chisq.test(table(matched$GROUP, matched$SEX))
  cat(sprintf("X² = %.3f, df = %d, p = %.4f\n", cs$statistic, cs$parameter, cs$p.value))
}


# ── 7. Matched sample demographics for all three methods ──────────────────────

for (res in list(results_A, results_B, results_C)) {
  matched <- match.data(res$matchit_obj)
  cat(paste0("\n=== Matched demographics: ", res$label, " ===\n"))
  matched %>%
    group_by(GROUP) %>%
    summarise(
      n          = n(),
      age_mean   = round(mean(AGE), 2),
      age_sd     = round(sd(AGE), 2),
      n_female   = sum(SEX == "F"),
      n_male     = sum(SEX == "M"),
      pct_female = round(mean(SEX == "F") * 100, 1)
    ) %>% print()
}


# ── 8. Save matched subject lists for all three methods ───────────────────────

filenames <- c(
  "matched_subsample_euclidean.csv",
  "matched_subsample_standard_mahal.csv",
  "matched_subsample_robust_mahal.csv"
)

for (i in seq_along(list(results_A, results_B, results_C))) {
  res     <- list(results_A, results_B, results_C)[[i]]
  matched <- match.data(res$matchit_obj)
  write.csv(
    matched[, c("CODE", "GROUP", "AGE", "SEX")],
    file.path(out_dir, filenames[i]),
    row.names = FALSE
  )
  cat(sprintf("✅  Saved: %s\n", filenames[i]))
}

# ============================================================
# Post-hoc power analysis — Toussas et al. (BRAINCOM-2026-214)
# Two-sample t-test, two-tailed, alpha = 0.05
# n1 = 32 (aMCI), n2 = 58 (HO)
# ============================================================

library(pwr)

n1    <- 58
n2    <- 32
alpha <- 0.05

# Harmonic mean
n_harm <- (2 * n1 * n2) / (n1 + n2)

# Minimum d for >=80% and >=95% power
find_d <- function(target_power) {
  uniroot(function(d) {
    pwr.t2n.test(n1 = n1, n2 = n2, d = d,
                 sig.level = alpha, alternative = "two.sided")$power - target_power
  }, lower = 0.3, upper = 1.5)$root
}

d_80 <- find_d(0.80)
d_95 <- find_d(0.95)

cat(sprintf("Harmonic mean n:           %.0f\n",   round(n_harm)))
cat(sprintf("Min d for >=80%% power:    %.2f\n",   ceiling(d_80 * 100) / 100))
cat(sprintf("Min d for >=95%% power:    %.2f\n",   ceiling(d_95 * 100) / 100))

