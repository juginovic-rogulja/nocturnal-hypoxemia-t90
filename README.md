# Nocturnal hypoxemia and long-term health: analysis code

Code for the paper by Alen Juginović and Dragana Rogulja on nocturnal hypoxemia (T90, the share of the night spent below 90% oxygen saturation) and future disease. The study computed 141 sleep-study parameters from 19,173 overnight laboratory polysomnograms of adults referred to two academic hospitals, followed each patient in the hospital record for 54 medical conditions and death, ranked the 141 parameters by how much each improves the held-out discrimination of future disease beyond age and sex, estimated hazard ratios for T90, and validated the associations in two community cohorts (the Sleep Heart Health Study and the Osteoporotic Fractures in Men Study). The code implements standard published methods: Welch's periodogram for spectral power (Welch 1967), spindle detection by a band-passed envelope threshold (Mölle et al. 2002), slow-oscillation detection by band-passed troughs with an amplitude criterion (Massimini et al. 2004), spindle-to-slow-oscillation coupling by the Hilbert phase (Helfrich et al. 2018), heart rate variability by the Task Force standard (Task Force of the European Society of Cardiology and the North American Society of Pacing and Electrophysiology 1996), the hypoxic burden (Azarbarzin et al. 2019), Cox proportional hazards regression (lifelines) and Benjamini-Hochberg correction (statsmodels).

The repository holds code only. No recording, table, number file or figure is included (see Data access).

## Layout

```
paths.py                 the one module every script reads its locations from (environment variables, documented defaults)
requirements.txt         the pinned Python packages
signal_features/         the April 2026 feature pipeline: the reader and the four feature families (brain spectral power,
                         brain microstructure, heart rate and variability, oximetry), the single-recording driver, the
                         cohort-scale EC2 driver and the merge that wrote aggregated_features_v5.csv
oximetry_extraction/     the September 2026 extraction chain that produced the frozen tables of the paper:
    X1_extraction/reex_code/   the same four feature families with the reader adapted to the regenerated file layout,
                               the per-night extractor (extract_night.py) and its EC2 worker
    X1_extraction/scripts/     the v8 pass: the pre-treatment window of split nights, the native-rate oximetry measures,
                               hypoxic burden, the periodic limb movement index, sleep and REM latency, the assembly of
                               the passes, the frozen hypoxic-burden definition, the synthetic tests
    X1_extraction/ec2/         the EC2 launch, the code tarball and the fetch of the v8 pass
                               (reex_code/SHA256.txt hashes the files as they ran, before the path rewrites of this
                               repository: extract_night.py, h5_reader.py, heart_features.py and lung_features.py still
                               match it, worker.py, compare_lib.py, brain_features.py and microstructure.py do not)
    cpap_stage/, X4_groupP/    the stage-resolved oximetry extraction on and off positive airway pressure (three designs)
    X2_outcomes/scripts/       the outcome definitions from diagnostic codes and the rebuild of first-diagnosis dates
    E2_tables/scripts/         the builder of the four frozen tables and its invariants
    E1_live/scripts/           the builder of measure_families.csv
    outcomes/                  the outcome summary
analysis/
    numbers/                 the statistical scripts of the paper (cohort definition, disease definitions, the measure
                             families, the ranking, the hazard ratios, the graded categories, the negative controls, the
                             cross-tabulation, the sleep-duration models, the treatment comparison, the external validation,
                             the proportional-hazards tests, the combination models, the Figure 2c row rule)
    sleep_variability/       the side analyses the paper reports (discrimination gain per disease and against clinical
                             baselines, sleep duration crossed with oxygenation, habitual sleep, stage-specific T90, the
                             saturation-threshold ladder, the overnight oxygen profile, residual hypoxemia on treatment)
    sleep_only_t90/          the exploratory sleep-only T90 ranking of 2 September 2026 (see the note below)
figures/                     the builders of the 37 figure sheets, kept in the round folders they ran in (ROUND54 built the
                             six main figures, ROUND49 to ROUND52 the Extended Data and Supplementary sheets), with the
                             shared modules the lane scripts import from earlier rounds
```

Some scripts that compute numbers live in the figure tree because they ran there (New_Figures, figures_alen, FINAL_FIGURES_2026-08-14, ROUND12, ROUND15, ROUND17, ROUND18 and ROUND30). They are listed in the report that accompanies the repository and are part of the analysis.

## Data access

The primary cohort comes from the Brain Data Science Platform (bdsp.io), available to qualified investigators under a data use agreement. The Sleep Heart Health Study and the Osteoporotic Fractures in Men Study are distributed by the National Sleep Research Resource (sleepdata.org) to investigators who complete its data access application. This repository therefore contains code and no data: no recording, no derived table, no number file, no figure and no identifier.

## Software

Python 3.14.6 with NumPy 2.4.2, SciPy 1.16.3, pandas 2.3.3, lifelines 0.30.1, statsmodels 0.14.6, patsy 1.0.2, joblib 1.5.3, duckdb 1.5.2, h5py 3.16.0, s3fs 2026.4.0, pyarrow 23.0.1 and matplotlib 3.10.8 (the versions read from the environment the analysis ran in). The figure assembly and verification additionally use pypdf 6.7.1, PyMuPDF 1.27.2.2, fonttools 4.61.1, pillow 12.1.1 and openpyxl 3.1.5, and one outcome-model script uses scikit-learn 1.8.0. `requirements.txt` pins them all.

```
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
```

Tools that are not Python packages and are not part of this repository: Ghostscript (`gs`, renders and text censuses of the figure sheets), tesseract (optical character checks of rendered sheets), Google Chrome (the schematic renderer), Node with the `docx` package (`npm install docx`, the figure document of `make_docx_r44.js`), and the AWS command line (the EC2 extraction).

Two Python packages the figure builders import are separate tools and are not included here. The schematic bands of Figures 1a, 1f, 3d, 4e, 5d and 6 and of Supplementary Fig. 16 were drawn with FigureLab (the `bioglyph` package), a separate drawing tool. The repository holds the scripts that call it (`figures/ROUND54_2026-09-28/lanes/LSKETCH/scripts/build_r54.py` and the round-40 and round-44 builders it extends). The data panels of some Extended Data and Supplementary sheets import the `prism_plotter` style package. Set `T90_FIGURELAB_DIR` and `T90_PRISM_PLOTTER_SRC` to their locations.

## Locations: paths.py

Every script starts with a two-line bootstrap that imports `paths.py` from the repository root, and no script names an absolute path of its own. `paths.py` reads its locations from environment variables, each with a documented default (the module docstring lists them). The main ones:

| variable | meaning | default |
|---|---|---|
| `T90_ROOT` | the data tree: frozen tables, extraction work, caches | `~/T90_work` |
| `T90_TABLES_DIR` | the four frozen v8 tables | `$T90_ROOT/data_frozen_v8_2026-09` |
| `T90_NUMBERS_DIR` | the number files the analysis writes and the figures read | `analysis/numbers` (beside the scripts) |
| `T90_FIGURE_ROOT` | the figure tree (scripts, base sheets, lane outputs) | `figures/` |
| `T90_PSG_DIR`, `BDSP_S3_ACCESS_POINT` | the recordings, local or on the platform's S3 access point | none for the access point |
| `T90_OMOP_CACHE_DIR` | the hospital-record cache | `$T90_ROOT/omop_cache` |
| `T90_SHHS_DIR`, `T90_MROS_DIR` | the National Sleep Research Resource files | under `$T90_ROOT/nsrr` |

Shell scripts read the same variables. Define them once with

```
eval "$(python3 paths.py)"
```

## How the stages run

1. Extraction of the signal parameters. `signal_features/` is the April 2026 implementation that computed the 122 signal-derived parameters (60 spectral, 39 microstructure, 16 heart, 7 oximetry) for the sleep atlas: `extract_one_patient.py` on one recording, `aws_ec2_launcher.sh` over the cohort (its worker is embedded in the instance user data), `aggregate_features.py` and `merge_v5.py` into `aggregated_features_v5.csv`. In September 2026 every night was re-extracted from the source files with `oximetry_extraction/X1_extraction/reex_code/` (the same feature functions, byte for byte apart from one comment line in `brain_features.py`, with the reader adapted to the regenerated file layout), run by `reex_code/ec2/worker.py`. The values in the frozen tables come from that re-extraction.
2. The v8 pass (`oximetry_extraction/X1_extraction/scripts/`). `build_manifests.py` writes the night lists, `extract_v8.py --families full` computes every parameter on the pre-treatment window of the 3,622 split nights and `--families light` computes the native-rate oximetry measures, hypoxic burden, the limb movement index and the latencies on all 19,173 nights (`launch_v8.sh` on EC2, `wait_and_fetch.py` to collect), `check_pilot_split.py` and `check_pilot_light.py` gate the pilots, `freeze_hb_definition.py` freezes the per-site hypoxic-burden windows, `assemble_split.py` and `assemble_light.py` assemble the passes, `second_path_check.py` recomputes a sample by an independent implementation, and `tests/test_synthetic.py` runs on a synthetic recording (`python tests/test_synthetic.py -v`, no data needed).
3. Outcomes (`X2_outcomes/scripts/`): `test_definitions_v8.py`, then `rebuild_outcomes_v8.py` (first-diagnosis dates from the hospital-record cache), then `check_outcomes_v8.py`.
4. Tables (`E2_tables/scripts/build_v8_tables.py --split ... --light ... --outcomes ... --out $T90_TABLES_DIR`), then `reconcile_v8.py` and `invariants_v8.py`.
5. Stage-resolved oximetry (`cpap_stage/`, `X4_groupP/`): `build_worklists.py`, `extract.py` on EC2, `mask_collapsed.py` and `assemble_groupP.py`, then `analysis/numbers/rebuild_cpap_stage.py` (which runs `cpap_stage/analyse.py` with `CPAP_STAGE_OUT` set).
6. Analysis (`analysis/numbers/`, outputs into `T90_NUMBERS_DIR`). `cohort_spec.py` defines the cohort, the tables and the counts and is imported by every script. `RUN_ORDER.md` lists every step of the recalculation in the order it ran, with the switches each step set. In brief: `rebuild_cpap_stage.py`, `freeze_01_cohorts.py` to `freeze_04_cpap_stage.py`, `build_ranking_v3.py` with `ranking_v3_noisefloor.py`, `ranking_v3_report.py`, `refresh_master_ranking_v3.py` and `refresh_master_cohort.py`, `run_ranking_noapnea.py`, `run_primary_unpenalized.py`, `run_all_v2.py`, `run_crosstab_v2.py`, `run_extras_v2.py`, `run_extras2_v2.py`, `run_combo_v2.py`, `run_treatment_v2.py`, `run_gaps_v2.py`, `run_mortality_v2.py`, `run_minutes_vs_pct_v2.py`, `run_wake_sleep_healthy.py`, `rebuild_nonresponder_phenotype.py`, the number-writing scripts in the figure tree, the `analysis/sleep_variability/` lanes (each folder in its numbered order, `dC_side_analyses/s00_positive_control.py` and `dC_per_disease/p00_positive_control.py` first), `make_negcontrols_final.py`, `run_treatment_v3.py`, `run_ph_all_v8.py`, the sleep-period T90 set, `ranking_v3_no_ohs.py`, `fig2c_selection_v8.py`, `summarise_outcomes_v8.py`, the full-nights and split-covariate sets, `run_causal_tests_all.py`, `negcontrols_panel_v7.py`, and last the external validation (`freeze_05_external.py`, `external_sleepduration/`). Two environment switches select the analysis population: `T90_FULL_NIGHTS_ONLY=1` for the ranking family and for every analysis that uses a sleep amount (the 15,551 full diagnostic nights), and `T90_COLUMN=spo2_pct_below_90_sleep` for the sleep-period T90 sensitivity set. `build_ranking_v3_splitcov.py` is the all-nights ranking with split-night status as a covariate.
7. Figures (`figures/`). Each lane folder holds numbered scripts (`01_build...`, `02_compose...`, `03_verify...`) or a runner (`run_chain_r54.sh`, `run_sheet.py`, `run_ed6.py`, `run_ed7.py`, `run_l5a.py`, `build_sheet.py`). ROUND54 built the six main figures (`lanes/L1` to `L5`, `LSKETCH` for the schematic bands and Figure 6, `raster_r54.py` for the raster twins of Figures 3 to 5, `compose_v31.py` and `final_check_r54.py` for the assembled set, `make_docx_r44.js` for the figure document). The Extended Data and Supplementary sheets were last rebuilt in ROUND49 (lanes `LED-A`, `LED-B`, `LSUPP-A` to `LSUPP-D`), with ROUND50 (`LSUPP-A`, `LSUPP-C`, `LSUPP-D`), ROUND51 (`LED-A`, `LED-B`) and ROUND52 (`LED-B`, `LSUPP-D`) superseding some sheets, and ROUND54 `LED7C` adding a panel to Extended Data Fig. 7.

## What can be run without the data

The synthetic tests of the v8 extraction and the byte-compilation and import of every module. The extraction needs the recordings (locally under `T90_PSG_DIR` or on the platform's S3 access point with credentials), the analysis needs the frozen tables and the hospital-record cache, the external validation needs the National Sleep Research Resource files, and the figure lanes need the number files and the base sheets and drawn-value records of the earlier rounds that they verify against.

## Note on the sleep-only T90 folder

`analysis/sleep_only_t90/` holds the exploratory ranking of 2 September 2026 in which T90 measured during sleep alone was entered beside the recording-wide measures. The sleep-period T90 sensitivity set reported in the paper was computed afterwards by the v8 chain (`T90_COLUMN=spo2_pct_below_90_sleep`, the `_sleepT90` outputs of `build_ranking_v3.py` and `run_primary_unpenalized.py`).

## References for the methods

Welch PD. The use of fast Fourier transform for the estimation of power spectra. IEEE Trans Audio Electroacoust, 1967, volume 15, pages 70 to 73. Mölle M, Marshall L, Gais S, Born J. Grouping of spindle activity during slow oscillations in human non-rapid eye movement sleep. J Neurosci, 2002, volume 22, pages 10941 to 10947. Massimini M, Huber R, Ferrarelli F, Hill S, Tononi G. The sleep slow oscillation as a traveling wave. J Neurosci, 2004, volume 24, pages 6862 to 6870. Helfrich RF, Mander BA, Jagust WJ, Knight RT, Walker MP. Old brains come uncoupled in sleep. Neuron, 2018, volume 97, pages 221 to 230. Task Force of the European Society of Cardiology and the North American Society of Pacing and Electrophysiology. Heart rate variability: standards of measurement, physiological interpretation and clinical use. Circulation, 1996, volume 93, pages 1043 to 1065. Azarbarzin A, Sands SA, Stone KL, et al. The hypoxic burden of sleep apnoea predicts cardiovascular disease-related mortality. Eur Heart J, 2019, volume 40, pages 1149 to 1157.

## License

MIT, copyright 2026 Alen Juginović and Dragana Rogulja.
