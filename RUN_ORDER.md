# Run order

The steps of the September 2026 recalculation in the order they ran, with the script each step called (its path in this repository) and the environment switches it set. The two switches select the analysis population: `T90_FULL_NIGHTS_ONLY=1` runs the step on the 15,551 full diagnostic nights (the ranking family and every analysis that uses a sleep amount), `T90_COLUMN=spo2_pct_below_90_sleep` runs it with the sleep-period T90 as the exposure. Steps that ran on EC2 (`wait_and_fetch.py`) collected the outputs of `launch_v8.sh`. The tables the steps read are named in `paths.py`.

| step | group | name | script | switches |
|---|---|---|---|---|
| 010 | extraction | build_manifests_v8 | oximetry_extraction/X1_extraction/scripts/build_manifests.py |  |
| 012 | extraction | probe12_reproduction_local | oximetry_extraction/X1_extraction/scripts/extract_v8.py, oximetry_extraction/X1_extraction/scripts/compare_to_reextract_final.py |  |
| 013 | extraction | pilot300_split_ec2 | oximetry_extraction/X1_extraction/ec2/wait_and_fetch.py, oximetry_extraction/X1_extraction/scripts/check_pilot_split.py |  |
| 014 | extraction | full3622_split_ec2 | oximetry_extraction/X1_extraction/ec2/wait_and_fetch.py |  |
| 015 | extraction | assemble_split_prepap | oximetry_extraction/X1_extraction/scripts/assemble_split.py |  |
| 020 | extraction | pilot300_light_local | oximetry_extraction/X1_extraction/scripts/extract_v8.py, oximetry_extraction/X1_extraction/scripts/check_pilot_light.py |  |
| 021 | extraction | freeze_hb_definition | oximetry_extraction/X1_extraction/scripts/freeze_hb_definition.py |  |
| 022 | extraction | full19173_light_ec2 | oximetry_extraction/X1_extraction/ec2/wait_and_fetch.py |  |
| 023 | extraction | assemble_light | oximetry_extraction/X1_extraction/scripts/assemble_light.py |  |
| 024 | extraction | reconcile_50_second_path | oximetry_extraction/X1_extraction/scripts/second_path_check.py |  |
| 030 | tables and definitions | measure_families_v8 | oximetry_extraction/E1_live/scripts/build_measure_families_v8.py |  |
| 040 | tables and definitions | disease_definitions_v8 | oximetry_extraction/X2_outcomes/scripts/test_definitions_v8.py |  |
| 041 | tables and definitions | rebuild_outcomes_v8 | oximetry_extraction/X2_outcomes/scripts/rebuild_outcomes_v8.py, oximetry_extraction/X2_outcomes/scripts/check_outcomes_v8.py |  |
| 050 | tables and definitions | build_v8_tables | oximetry_extraction/E2_tables/scripts/build_v8_tables.py |  |
| 051 | tables and definitions | invariants_v8 | oximetry_extraction/E2_tables/scripts/invariants_v8.py |  |
| 052 | tables and definitions | reconcile_sample_vs_source | oximetry_extraction/E2_tables/scripts/reconcile_v8.py |  |
| 100 | shard-dependent steps | rebuild_cpap_stage | analysis/numbers/rebuild_cpap_stage.py |  |
| 101 | the inherited chain | freeze_01_cohorts | analysis/numbers/freeze_01_cohorts.py |  |
| 102 | the inherited chain | freeze_02_demographics_sanity | analysis/numbers/freeze_02_demographics.py |  |
| 103 | the inherited chain | freeze_03_primary | analysis/numbers/freeze_03_primary.py | SUPERSEDED penalized specification, no printed number, exits unless T90_RUN_SUPERSEDED=1 |
| 104 | the inherited chain | freeze_04_cpap_stage | analysis/numbers/freeze_04_cpap_stage.py | SUPERSEDED penalized specification, no printed number, exits unless T90_RUN_SUPERSEDED=1 |
| 105 | the inherited chain | build_ranking_v3_current | analysis/numbers/build_ranking_v3.py | T90_FULL_NIGHTS_ONLY=1 |
| 106 | the inherited chain | ranking_v3_noisefloor | analysis/numbers/ranking_v3_noisefloor.py | T90_FULL_NIGHTS_ONLY=1 |
| 107 | the inherited chain | ranking_v3_report_sanity | analysis/numbers/ranking_v3_report.py | T90_FULL_NIGHTS_ONLY=1 |
| 108 | the inherited chain | refresh_master_ranking_v3 | analysis/numbers/refresh_master_ranking_v3.py | T90_FULL_NIGHTS_ONLY=1 |
| 109 | the inherited chain | refresh_master_cohort | analysis/numbers/refresh_master_cohort.py |  |
| 110 | the inherited chain | run_ranking_noapnea | analysis/numbers/run_ranking_noapnea.py | T90_FULL_NIGHTS_ONLY=1 |
| 111 | the inherited chain | run_primary_unpenalized | analysis/numbers/run_primary_unpenalized.py |  |
| 112 | the inherited chain | run_all_v2 | analysis/numbers/run_all_v2.py |  |
| 113 | the inherited chain | run_crosstab_v2 | analysis/numbers/run_crosstab_v2.py |  |
| 114 | the inherited chain | run_extras_v2 | analysis/numbers/run_extras_v2.py |  |
| 115 | the inherited chain | run_extras2_v2 | analysis/numbers/run_extras2_v2.py |  |
| 116 | the inherited chain | run_combo_v2 | analysis/numbers/run_combo_v2.py | T90_FULL_NIGHTS_ONLY=1 |
| 117 | the inherited chain | run_treatment_v2 | analysis/numbers/run_treatment_v2.py |  |
| 118 | the inherited chain | run_gaps_v2 | analysis/numbers/run_gaps_v2.py |  |
| 119 | the inherited chain | run_mortality_v2 | analysis/numbers/run_mortality_v2.py |  |
| 120 | the inherited chain | run_minutes_vs_pct_v2 | analysis/numbers/run_minutes_vs_pct_v2.py | T90_FULL_NIGHTS_ONLY=1 |
| 121 | shard-dependent steps | run_wake_sleep_healthy_refresh | analysis/numbers/run_wake_sleep_healthy.py |  |
| 122 | shard-dependent steps | rebuild_nonresponder_phenotype | analysis/numbers/rebuild_nonresponder_phenotype.py |  |
| 123 | the inherited chain | compute_alenfig3 | figures/figures_alen/scripts/compute_alenfig3.py |  |
| 124 | the inherited chain | fits_one_exposure | figures/New_Figures/_NOT_THESE_workfiles/scripts/fits_one_exposure.py |  |
| 125 | shard-dependent steps | fits_sleep_unadjusted | figures/New_Figures/_NOT_THESE_workfiles/scripts/fits_sleep_unadjusted.py |  |
| 126 | the inherited chain | fits_severe_apnea_negctrl | figures/New_Figures/_NOT_THESE_workfiles/scripts/fits_severe_apnea_negctrl.py |  |
| 127 | the inherited chain | fits_phenotype_combined_v2 | figures/New_Figures/_NOT_THESE_workfiles/scripts/fits_phenotype_combined_v2.py |  |
| 128 | the inherited chain | prep_figure1_profile | figures/New_Figures/_NOT_THESE_workfiles/scripts/prep_figure1_profile.py |  |
| 129 | the inherited chain | make_efig18_shortsleep_oxygen | figures/New_Figures/_NOT_THESE_workfiles/scripts/make_efig18_shortsleep_oxygen.py | T90_FULL_NIGHTS_ONLY=1 |
| 130 | the inherited chain | lagladder | figures/New_Figures/_NOT_THESE_workfiles/scripts/lagladder.py |  |
| 1305 | the inherited chain | lagladder_fullnights | figures/New_Figures/_NOT_THESE_workfiles/scripts/lagladder.py | T90_FULL_NIGHTS_ONLY=1 |
| 131 | the inherited chain | compute_fig3_within_stratum | figures/FINAL_FIGURES_2026-08-14/_scripts/compute_fig3_within_stratum.py |  |
| 132 | shard-dependent steps | freeze_doseresponse | figures/figures_alen/scripts/freeze_doseresponse.py |  |
| 133 | the inherited chain | q2_noapnea_freq_model | analysis/sleep_variability/q2_noapnea_freq/model.py |  |
| 134 | the inherited chain | q8_after_oxygen_joint | analysis/sleep_variability/q8_after_oxygen/joint.py | T90_FULL_NIGHTS_ONLY=1 |
| 135 | the inherited chain | dur_x_oxygen_model | analysis/sleep_variability/dur_x_oxygen/model.py | T90_FULL_NIGHTS_ONLY=1 |
| 136 | the inherited chain | dur_x_oxygen_6cell_tst300_model | analysis/sleep_variability/dur_x_oxygen_6cell_tst300/model.py | T90_FULL_NIGHTS_ONLY=1 |
| 137 | the inherited chain | dur_x_oxygen_6cell_model | analysis/sleep_variability/dur_x_oxygen_6cell/model.py | T90_FULL_NIGHTS_ONLY=1 |
| 138 | the inherited chain | dur_x_oxygen_healthy_model | analysis/sleep_variability/dur_x_oxygen_healthy/model.py | T90_FULL_NIGHTS_ONLY=1 |
| 139 | shard-dependent steps | dur_x_oxygen_healthy_attack | analysis/sleep_variability/dur_x_oxygen_healthy/attack.py |  |
| 140 | the inherited chain | shortsleep_before_after_t90 | analysis/sleep_variability/shortsleep_before_after_t90/run_shortsleep_before_after_t90.py |  |
| 141 | the inherited chain | tst_ceiling_model | analysis/sleep_variability/tst_ceiling/model.py | T90_FULL_NIGHTS_ONLY=1 |
| 142 | the inherited chain | tst_ceiling_attack | analysis/sleep_variability/tst_ceiling/attack.py | T90_FULL_NIGHTS_ONLY=1 |
| 143 | the inherited chain | make_efig15_tst_ceiling | figures/New_Figures/_NOT_THESE_workfiles/scripts/make_efig15_tst_ceiling.py | T90_FULL_NIGHTS_ONLY=1 |
| 144 | shard-dependent steps | stage_specific_build | analysis/sleep_variability/stage_specific/build.py |  |
| 145 | shard-dependent steps | stage_specific_model_general | analysis/sleep_variability/stage_specific/model_general.py |  |
| 146 | shard-dependent steps | stage_specific_attack | analysis/sleep_variability/stage_specific/attack.py |  |
| 1461 | shard-dependent steps | stage_specific_model_healthy_lane | figures/ROUND30_2026-09-08/V7_L5b_PAP_SUPP/Supp_Fig20/scripts/model_healthy_v7_lane.py |  |
| 1463 | shard-dependent steps | stage_specific_negcontrols_v5 | figures/New_Figures/_NOT_THESE_workfiles/scripts/fits_stage_negcontrols_v5.py |  |
| 147 | the inherited chain | pap_residual_run_extension | analysis/sleep_variability/pap_residual_continuous/run_extension.py |  |
| 148 | the inherited chain | pap_residual_delta_run | analysis/sleep_variability/pap_residual_continuous/delta_run.py |  |
| 149 | the inherited chain | habitual_outcomes_v2_build | analysis/sleep_variability/habitual_outcomes_v2/build.py |  |
| 150 | the inherited chain | habitual_outcomes_v2_model_habitual | analysis/sleep_variability/habitual_outcomes_v2/model_habitual.py |  |
| 151 | the inherited chain | habitual_outcomes_v2_attack | analysis/sleep_variability/habitual_outcomes_v2/attack.py |  |
| 152 | shard-dependent steps | oxygen_profile_assemble | analysis/sleep_variability/oxygen_profile/assemble.py |  |
| 153 | shard-dependent steps | oxygen_profile_b01_profile | analysis/sleep_variability/oxygen_profile/b01_profile.py |  |
| 154 | the inherited chain | placebo_restoration_rate_b91 | analysis/sleep_variability/final_verification/compute_b5_b9.py |  |
| 155 | the inherited chain | round12_L4_refit_fig5b | figures/ROUND12_2026-09-03/L4_fig5/L4_refit_fig5b.py |  |
| 156 | the inherited chain | round18_j1_build_fig3b | figures/ROUND18_2026-09-04/LC_panels/work/j1_build.py |  |
| 157 | the inherited chain | round15_L2_km_compute_fig1d | figures/ROUND15_2026-09-04/L2_fig1d/L2_km_compute.py |  |
| 158 | the inherited chain | fits_graded_newcontrols | figures/New_Figures/_NOT_THESE_workfiles/scripts/fits_graded_newcontrols.py |  |
| 159 | the inherited chain | fits_demographics_rank | figures/New_Figures/_NOT_THESE_workfiles/scripts/fits_demographics_rank.py | T90_FULL_NIGHTS_ONLY=1 |
| 160 | the inherited chain | dC_side_s00_positive_control | analysis/sleep_variability/dC_side_analyses/s00_positive_control.py | T90_FULL_NIGHTS_ONLY=1 |
| 161 | the inherited chain | dC_side_s01_organ_decomposition | analysis/sleep_variability/dC_side_analyses/s01_organ_decomposition.py | T90_FULL_NIGHTS_ONLY=1 |
| 162 | shard-dependent steps | dC_side_s02_sleep_vs_whole | analysis/sleep_variability/dC_side_analyses/s02_sleep_vs_whole.py |  |
| 163 | the inherited chain | dC_side_s03_combinations | analysis/sleep_variability/dC_side_analyses/s03_combinations.py | T90_FULL_NIGHTS_ONLY=1 |
| 164 | the inherited chain | dC_side_s04_nonlinear | analysis/sleep_variability/dC_side_analyses/s04_nonlinear.py | T90_FULL_NIGHTS_ONLY=1 |
| 165 | the inherited chain | dC_side_s05_restricted_ranking | analysis/sleep_variability/dC_side_analyses/s05_restricted_ranking.py | T90_FULL_NIGHTS_ONLY=1 |
| 166 | the inherited chain | dC_side_s06_agebasis | analysis/sleep_variability/dC_side_analyses/s06_agebasis.py | T90_FULL_NIGHTS_ONLY=1 |
| 167 | shard-dependent steps | dC_side_s07_hazard_ratios | analysis/sleep_variability/dC_side_analyses/s07_hazard_ratios.py |  |
| 168 | the inherited chain | dC_side_s08_hr_matrix | analysis/sleep_variability/dC_side_analyses/s08_hr_matrix.py |  |
| 169 | the inherited chain | dC_side_s09_ranked_metrics | analysis/sleep_variability/dC_side_analyses/s09_ranked_metrics.py |  |
| 170 | the inherited chain | dC_side_s99_assemble | analysis/sleep_variability/dC_side_analyses/s99_assemble.py |  |
| 171 | the inherited chain | dC_per_disease_p00_positive_control | analysis/sleep_variability/dC_per_disease/p00_positive_control.py | T90_FULL_NIGHTS_ONLY=1 |
| 172 | the inherited chain | dC_per_disease_p06_lag2_gate | analysis/sleep_variability/dC_per_disease/p06_lag2_gate.py | T90_FULL_NIGHTS_ONLY=1 |
| 173 | the inherited chain | dC_per_disease_p07_points | analysis/sleep_variability/dC_per_disease/p07_points.py | T90_FULL_NIGHTS_ONLY=1 |
| 174 | the inherited chain | dC_per_disease_p01_bootstrap_standard | analysis/sleep_variability/dC_per_disease/p01_bootstrap.py | T90_FULL_NIGHTS_ONLY=1 |
| 175 | the inherited chain | dC_per_disease_p01_bootstrap_lag2 | analysis/sleep_variability/dC_per_disease/p01_bootstrap.py | T90_FULL_NIGHTS_ONLY=1 |
| 176 | the inherited chain | dC_per_disease_p02_assemble | analysis/sleep_variability/dC_per_disease/p02_assemble.py | T90_FULL_NIGHTS_ONLY=1 |
| 177 | the inherited chain | dC_per_disease_p03_figure | analysis/sleep_variability/dC_per_disease/p03_figure.py | T90_FULL_NIGHTS_ONLY=1 |
| 178 | the inherited chain | dC_per_disease_p03b_split | analysis/sleep_variability/dC_per_disease/p03b_split.py | T90_FULL_NIGHTS_ONLY=1 |
| 179 | the inherited chain | dC_per_disease_p04_report | analysis/sleep_variability/dC_per_disease/p04_report.py | T90_FULL_NIGHTS_ONLY=1 |
| 180 | the inherited chain | dC_per_disease_p05_verify | analysis/sleep_variability/dC_per_disease/p05_verify.py | T90_FULL_NIGHTS_ONLY=1 |
| 181 | the inherited chain | round17_p00_selftest | figures/ROUND17_T90_GT10_CONCORDANCE/p00_selftest.py | T90_FULL_NIGHTS_ONLY=1 |
| 182 | the inherited chain | round17_p01_points | figures/ROUND17_T90_GT10_CONCORDANCE/p01_points.py | T90_FULL_NIGHTS_ONLY=1 |
| 183 | the inherited chain | round17_p02_bootstrap_lag0 | figures/ROUND17_T90_GT10_CONCORDANCE/p02_bootstrap.py | T90_FULL_NIGHTS_ONLY=1 |
| 184 | the inherited chain | round17_p02_bootstrap_lag2 | figures/ROUND17_T90_GT10_CONCORDANCE/p02_bootstrap.py | T90_FULL_NIGHTS_ONLY=1 |
| 185 | the inherited chain | round17_p03_assemble | figures/ROUND17_T90_GT10_CONCORDANCE/p03_assemble.py | T90_FULL_NIGHTS_ONLY=1 |
| 186 | the inherited chain | round17_p04_figure | figures/ROUND17_T90_GT10_CONCORDANCE/p04_figure.py | T90_FULL_NIGHTS_ONLY=1 |
| 187 | the inherited chain | round17_p05_report | figures/ROUND17_T90_GT10_CONCORDANCE/p05_report.py | T90_FULL_NIGHTS_ONLY=1 |
| 188 | shard-dependent steps | dC_clinical_run_clinical_baselines | analysis/sleep_variability/dC_clinical_baselines/run_clinical_baselines.py | T90_FULL_NIGHTS_ONLY=1 |
| 189 | shard-dependent steps | dC_clinical_run_tst_baselines | analysis/sleep_variability/dC_clinical_baselines/run_tst_baselines.py | T90_FULL_NIGHTS_ONLY=1 |
| 190 | the inherited chain | dC_clinical_run_home_baselines | analysis/sleep_variability/dC_clinical_baselines/run_home_baselines.py | T90_FULL_NIGHTS_ONLY=1 |
| 191 | the inherited chain | dC_clinical_run_external_dc_tst_sanity | analysis/sleep_variability/dC_clinical_baselines/run_external_dc_tst.py |  |
| 192 | the inherited chain | dC_clinical_run_external_dc_agesex_sanity | analysis/sleep_variability/dC_clinical_baselines/run_external_dc_agesex.py |  |
| 193 | shard-dependent steps | dC_clinical_make_comparison | analysis/sleep_variability/dC_clinical_baselines/make_comparison.py |  |
| 194 | shard-dependent steps | dC_clinical_make_three_settings_workbook | analysis/sleep_variability/dC_clinical_baselines/make_three_settings_workbook.py |  |
| 195 | the inherited chain | t85_a01_banded | analysis/sleep_variability/t85_analysis/a01_banded.py |  |
| 196 | the inherited chain | t85_a02_persd | analysis/sleep_variability/t85_analysis/a02_persd.py |  |
| 197 | the inherited chain | t85_a03_dC | analysis/sleep_variability/t85_analysis/a03_dC.py | T90_FULL_NIGHTS_ONLY=1 |
| 198 | the inherited chain | t85_a04_report | analysis/sleep_variability/t85_analysis/a04_report.py |  |
| 199 | the inherited chain | t_ladder_b01_banded | analysis/sleep_variability/t_ladder_up/b01_banded.py |  |
| 200 | the inherited chain | t_ladder_b02_persd | analysis/sleep_variability/t_ladder_up/b02_persd.py |  |
| 201 | the inherited chain | t_ladder_b03_dC | analysis/sleep_variability/t_ladder_up/b03_dC.py | T90_FULL_NIGHTS_ONLY=1 |
| 202 | the inherited chain | t_ladder_b04_summary | analysis/sleep_variability/t_ladder_up/b04_summary.py |  |
| 203 | the inherited chain | t_ladder_b05_figure | analysis/sleep_variability/t_ladder_up/b05_figure.py |  |
| 204 | the inherited chain | t_ladder_b06_report | analysis/sleep_variability/t_ladder_up/b06_report.py |  |
| 211 | new or modified analyses of the v8 recalculation | run_treatment_v3 | analysis/numbers/run_treatment_v3.py |  |
| 212 | new or modified analyses of the v8 recalculation | refit_fig5b_adjusted | figures/ROUND12_2026-09-03/L4_fig5/L4_refit_fig5b.py |  |
| 220 | new or modified analyses of the v8 recalculation | run_ph_all_v8 | analysis/numbers/run_ph_all_v8.py |  |
| 230 | new or modified analyses of the v8 recalculation | ranking_sleepT90 | analysis/numbers/build_ranking_v3.py | T90_FULL_NIGHTS_ONLY=1 T90_COLUMN=spo2_pct_below_90_sleep |
| 231 | new or modified analyses of the v8 recalculation | primary_sleepT90 | analysis/numbers/run_primary_unpenalized.py | T90_COLUMN=spo2_pct_below_90_sleep |
| 232 | new or modified analyses of the v8 recalculation | ranking_without_ohs | analysis/numbers/ranking_v3_no_ohs.py |  |
| 233 | new or modified analyses of the v8 recalculation | fig2c_top10_rule | analysis/numbers/fig2c_selection_v8.py |  |
| 240 | new or modified analyses of the v8 recalculation | outcomes_added_disclosure | oximetry_extraction/outcomes/summarise_outcomes_v8.py |  |
| 235 | new or modified analyses of the v8 recalculation | primary_fullnights | analysis/numbers/run_primary_unpenalized.py | T90_FULL_NIGHTS_ONLY=1 |
| 236 | new or modified analyses of the v8 recalculation | ranking_allnights_splitcov | analysis/numbers/build_ranking_v3_splitcov.py |  |
| 237 | new or modified analyses of the v8 recalculation | habitual_concordance_fullnights | analysis/sleep_variability/habitual_outcomes_v2/model_habitual.py | T90_FULL_NIGHTS_ONLY=1 |
| 250 | new or modified analyses of the v8 recalculation | run_causal_tests_all | analysis/numbers/run_causal_tests_all.py |  |
| 205 | the inherited chain | make_negcontrols_final | analysis/numbers/make_negcontrols_final.py |  |
| 251 | new or modified analyses of the v8 recalculation | negcontrols_panel_v8 | analysis/numbers/negcontrols_panel_v7.py |  |
| 252 | new or modified analyses of the v8 recalculation | fit_joint_t90_tst_landmark | analysis/sleep_variability/reorg_fits_efig4_lag2/fit_joint_t90_tst_landmark.py |  |
| 253 | new or modified analyses of the v8 recalculation | habitual_outcomes_v2_model_crossed | analysis/sleep_variability/habitual_outcomes_v2/model_crossed.py |  |
| 254 | new or modified analyses of the v8 recalculation | habitual_ref79_model_ref79 | analysis/sleep_variability/habitual_ref79/model_ref79.py |  |
| 260 | new or modified analyses of the v8 recalculation | freeze_05_external_v8_3 | analysis/numbers/freeze_05_external.py |  |
| 261 | new or modified analyses of the v8 recalculation | external_sleepduration_v8_3 | analysis/sleep_variability/external_sleepduration/run_external_sleepduration.py |  |
| 262 | new or modified analyses of the v8 recalculation | external_sleepduration_persd_v8_3 | analysis/sleep_variability/external_sleepduration/run_external_sleepduration_persd.py |  |

Steps of the recalculation that are tooling of that machine and are not in the repository: 001 preflight_v8, 002 recover_reex_code, 011 write_extract_v8, 031 sweep_197_to_141, 053 repoint_v8, 054 park_caches_v8, 300 compare_numbers_v8, 301 classify_v8, 302 results_for_alen, 303 provenance_audit_v8.
