# AGP_confounders

Presented are scripts performing core tasks in the study "Host variables confound gut microbiota studies of human disease" by Vujkovic-Cvijin, I. et al. 2020. Below is a summary of script functions.

1_unmatch_nperms.R:
Case-control cohorts are constructed with controls chosen only by proximity to cases by location metadata (termed 'unmatched') for all input conditions. 25 permutations of unmatched cohorts are created for each condition.

2_fishwilx_unmatched_cohorts.R:
Significance of differences between cases and controls are tested. Benjamini-Hochberg false discovery rate-corrected Q values are calculated to determine which input confounder variables differ in their distributions between cases and controls.

3_unmatch_match_nperms_adonis.R:
25 permutations of cohorts are created for 'unmatched' case-control cohorts and case-control cohorts matched for all confounder variables found to differ between cases in controls in the preceding script. PERMANOVA P and F statistic values are calculated and graphical visualizations are produced which display differences in matched vs. unmatched cohorts.
