# Multimodal morphometry in amnestic MCI — analysis code

Code for: Toussas K, Marie D, ..., Bréchet L. *Multimodal morphometry reveals hippocampal atrophy and temporoparietal thinning in mild cognitive impairment.* Brain Communications (2026). Repository: https://github.com/MemStim-Lab/toussas-2026-morphometry-amci

Participant data (MRI, cognitive scores) are not included; see the Data availability statement in the paper. Full acquisition, preprocessing and statistical parameters are described in the Methods section of the paper; this README lists the software versions, the files, and the order in which they are run.

## Software

MATLAB R2023a with SPM12 and CAT12 v12.8.2 including its TFCE toolbox; R with MatchIt, dplyr, readxl, pwr, ggplot2, ggpubr, rstatix; Python 3.10 with numpy, pandas, scipy, statsmodels, pingouin.

## Files and run order

1. **Preprocessing** — CAT12 cross-sectional segmentation (SANLM denoising, DARTEL normalisation to MNI at 1.5 mm, modulation, 8 mm FWHM smoothing) and surface pipeline (cortical thickness, resampled to the 32k mesh, 12 mm FWHM smoothing). Run in the CAT12 GUI with default settings; no batch file was saved.
2. `02_group_statistics/2_2_vbm_two_sample_ttest_TIV_age_sex.m` — SPM12 two-sample t-test on the smoothed modulated GM maps (HO vs aMCI, covariates TIV, age, sex). Participants and covariates are read from `demographics.csv` (template in `data_templates/`). Inference: cluster-level FWE p < 0.05, k ≥ 30 voxels. The SBM model (two-sample t-test on thickness maps, covariate age; vertex-wise permutation-based FWE p < 0.05 from the CAT12 TFCE toolbox, ≥ 30 vertices) was run with the same design in the CAT12 GUI.
3. `03_cluster_volumes/get_totals.m` — sums GM volume (ml) within each FWE-significant cluster mask for every participant. `2_4_group_diff_GMV_HO_vs_MCI.R` compares and plots the cluster volumes between groups.
4. `04_matched_subsample/matching_subsample.R` — selects the 32 HO matched 1:1 on age and sex to the aMCI group (MatchIt; the robust Mahalanobis-distance solution was used). The VBM and SBM models of step 2 were then re-run on these 64 participants (Supplementary Figure S2, Table S4).
5. `05_brain_behavior/3_1_correlation_matrix_FDR.ipynb` — Pearson / Kendall correlations between hippocampal cluster volume and the 17 cognitive measures, Benjamini–Hochberg FDR. `3_2_correlation_plots.R` — scatter plots of Figure 3.

Group-level result maps for Supplementary Figure S2 C–D are in `results_maps/`, and the CAT12 display settings used for those panels in `figures_display_settings/`.

## Contact

Lucie Bréchet, University of Geneva — lucie.brechet@unige.ch
