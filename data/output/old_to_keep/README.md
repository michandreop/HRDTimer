# old_to_keep

Files from data/output/ carrying information NOT present in the new supplemental tables
(figures/supplemental_tables_v2). Original relative paths preserved.

- `HRDTimer_filters.csv` — Per-sample early/late/NA SBS1&SBS3 exposure bounds (QC-filter tracks, Fig S13). NOT in supp tables.
- `SuppFig7_exposure_summary.csv` — Per-sample early/late/NA exposures for ALL signatures + CIs. Supp tables tabulate only SBS1/SBS3.
- `PCAWG_SCANB_INFORM_SBS1_Age_plot_inc_sub.csv` — SBS1 burden INCLUDING subclonal/cell-private (differs from SupTable10 = MRCA/non-sub).
- `Breast_cohort_sample_removal_reasons.csv` — Free-text QC removal reasons. Supp tables carry only a boolean flag.
- `May27_TimingRun/PCAWG_Ovary_WGD_HRD_TimingResults_timing_nboot200_0001_only_prob_change_boot_pSub.csv` — Ovary pi1/pi2, c, Nt per sample. SupTable9 (ovary) has only HRDTime/WGDTime.
- `Metastasis_SBS1_Age_plot.csv` — Primary-relapse: extra cols not in SupTable11 (purity, raw CpG counts, shared, Relapse_time).
- `NormalTissue_SBS1_Age_plot.csv` — Normal tissue: extra cols not in SupTable10 (raw SBS1, SNV/indel counts, case).
