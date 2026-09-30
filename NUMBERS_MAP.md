# How to check a number the paper prints

`NUMBERS_MAP.tsv`, beside this file, has one row per number printed in the Abstract, the Results and the
Discussion. Each row says which output file the number comes out of, how to find the row inside that file,
which column or json key holds it, which script of this repository wrote that file, the command that
rebuilds it, and whether the value in the file matches the printed number.

A reviewer who wants to check one number reads its row, opens the named output file at the named row and
column, and compares. A reviewer who wants to regenerate the number runs the command in `run_command`.

## The columns

| column | what it holds |
|---|---|
| `section` | Abstract, Results or Discussion, the section the number is printed in |
| `number` | the number exactly as the paper prints it, including the comma and the per cent sign |
| `claim` | the short phrase the number asserts |
| `output_file` | the output the value is read from (see **Where the outputs live**) |
| `row_key` | how to find the row: a column test for a csv, a dotted key for a json, or `-` for a scalar |
| `column` | the column, the dotted json key, or the rule that derives the value from the rows |
| `script` | the script of this repository that writes that output, repository-relative |
| `run_command` | the command that rebuilds it, with the environment switch the step needs and its `RUN_ORDER.md` step id |
| `verified` | what happened when the value was read out of the file (see **The verdicts**) |

## Where the outputs live

`output_file` uses the environment-variable names of `paths.py`, so a path is read the same way a script
reads it:

| prefix | meaning | default |
|---|---|---|
| a bare file name | a file in `$T90_NUMBERS_DIR`, the analysis output folder | `analysis/numbers/` |
| `$T90_SV_ROOT/...` | an output of an `analysis/sleep_variability/` lane | beside the lane scripts |
| `$T90_FIGURE_ROOT/...` | a figure round's own output (a panel's drawn-value record or a lane roll-up) | `figures/` |
| `$T90_TABLES_DIR/...` | one of the four frozen v8 tables, or its provenance sidecar | `$T90_ROOT/data_frozen_v8_2026-09` |
| `$REPO/...` | a file shipped in this repository (`measure_families.csv`, `disease_definitions.py`) | the repository root |

Set them once, as `README.md` describes:

```
export T90_TABLES_DIR=/path/to/data_frozen_v8_2026-09
export T90_NUMBERS_DIR=/path/to/where/outputs/should/go
export T90_SV_ROOT=/path/to/sleep_variability/outputs
export T90_FIGURE_ROOT=/path/to/figure/rounds
eval "$(python3 paths.py)"
```

Then a `run_command` runs from the repository root. `run_command` carries the switch the step sets, because
getting it wrong silently changes the cohort: `T90_FULL_NIGHTS_ONLY=1` selects the 15,551 full diagnostic
nights, no switch keeps all 19,173 recordings, and `T90_COLUMN=spo2_pct_below_90_sleep` swaps the exposure
for the sleep-period sensitivity set. The `# step NNN` comment is the row's step in `RUN_ORDER.md`.

`measure_families.csv` is a repository file, not an output, and `analysis/numbers/cohort_spec.py` addresses
it through `T90_NUMBERS_DIR`. A reviewer who repoints that variable must copy
`analysis/numbers/measure_families.csv` into the folder it points at, or the ranking report and two
side-analysis lanes will not find it.

## The verdicts

`verified` reports what was read out of the file, never an inference. It takes one of these shapes:

| shape | meaning |
|---|---|
| `YES` | the file holds the printed value exactly, to the last digit the paper shows |
| `YES (a setting, proved present in the named output file)` | the number is a cut-off or a level rather than an estimate, and the named file proves the analysis used it |
| `YES (bound holds: the file has X, under/over the printed N)` | the paper states the number as a bound ("at most", "no more than", "within"), and the file's value satisfies it |
| `YES (the file has X; the paper prints it as an approximation, not a rounded value)` | the paper writes "~140" or "about 140" for the 141 in the file |
| `ROUNDS (X)` | the file holds `X`, which equals the printed number after half-up rounding to the printed precision |
| `MISMATCH (...)` | the file does not hold the printed number, with the value it does hold |
| `CITED (...)` | a value quoted from the literature, with its reference number: no output of this code produces it |

Rounding is half up at the precision the paper prints, which is the rule the figure builders use
(`v14lib.halfup`). A value that carries more digits than the paper shows is `ROUNDS`, not `YES`: `1.60`
against `1.6018551142658137` is a correct print, not an exact one, and the map says so.

## Which numbers are in the map, and which are not

The map covers every number in the Abstract, the Results and the Discussion, including the figure legends,
which sit inside the Results.

Some numbers appear many times. `19,173`, `141`, `52`, `53`, `55`, `95%`, `1 SD`, `1,894`, `18`, `11` and
above all the `90` of "T90" recur in almost every paragraph. The map carries one row per distinct claim, not
one row per occurrence, as the `claim` column makes explicit.

Four numbers in the Discussion are `CITED` because they are taken from the literature and no code
here produces them: the published discrimination gains `0.039` for hemoglobin A1c and `0.019` for NT-proBNP
(references 26 and 27), and `845,569` Medicare sleep studies in 2014 with the `1` million a year today
(reference 44).

Nothing else was left out. Fifteen numerals in the target text are Supplementary Table or Figure numbers
that read as values once a sentence is scanned ("Supplementary Tables 19 and 20"), and one is a pair of
reference superscripts that joins onto the measure's name ("T90.47,48"). Those are cross-references, not
results.

## Traps a reviewer will hit

**Two files hold the same quantity and disagree, because one is older.** The Results sentence "In 4,112
patients studied before and during positive airway pressure (PAP) therapy, median T90 decreased from 8.0%
to 0.4%" is `treatment_v2.json` (step 117), `effect_by_stage.sleep`, on 4,112 patients. The v8.2 rebuild of
the stage-resolved table computes the same quantity in `cpap_stage.json` (step 104) as 9.46% falling to
0.41% on 3,747 patients, because that pass masks nights whose staging collapsed. The two are different
samples, not the same number twice. The map names `treatment_v2.json`, which is the file the sentence was
written from.

**`measure_families.csv` carries a stale ranking column.** Its `dC_ranking_v3` for `spo2_pct_below_90` is
`0.0203387549724979`, while `ranking_v3.csv` has `0.0201150933051042`. Both round to the printed `0.020`.
The ranking file is the source; the family lookup is for the family names.

**The threshold ladder's own summary counts over 53 outcomes, the paper over 52.** The Discussion prints
"18 to 38 of the 52 ranked outcomes" for the saturation-threshold ladder. `LADDER_SUMMARY.csv` and
`persd_summary.json` store 18, 28, 34, 35, 40 and 38 for the six rungs, because their family is the 53 real
graded outcomes. Counted over the 52 rows of `ladder_persd.csv` that carry `in_ranking_48 == True`, which is
the denominator the sentence names, the rungs give 18, 27, 33, 34, 38 and 37, so the range over the five
alternative thresholds is 18 to 38 as printed. The map names the 52-row rule.

**`gaps_v2.json` holds two graded counts and the paper prints the stricter one.** "A graded increase across
the four T90 categories was observed for 27 of 53 conditions" is `graded.n_monotone_top_significant`, the
conditions whose hazard ratios never fall across the four categories and whose top category's confidence
interval clears 1. The neighbouring `graded.n_monotone`, which drops the confidence requirement, is 30. A
reviewer who reads the first key they see will find 30, not 27.

**Three scripts write to the current working directory and carry no provenance sidecar.**
`run_all_v2.py`, `run_extras_v2.py` and `run_treatment_v2.py` write `results_v2.json`, `extras_v2.json` and
`treatment_v2.json` with a bare relative filename and do not call `cohort_spec.sidecar`. Run them from the
folder `T90_NUMBERS_DIR` points at, or the output lands elsewhere. Sixty numbers of the map come out of
those three files.

**"at most 0.003 in the C statistic" for total sleep time is a ceiling, and two files round to it.** The map
names `tst_ceiling.json`, `ladder[0].dC_ceiling` = 0.00296: the most total sleep time could gain once its
night-to-night variability is corrected for, which is what "at most" states. The ranking's own average gain
for `TST_min` is -0.0011, below zero, so the printed number is not a gain. A second value also rounds to 0.003
and is not the one: the largest gain over the four short-sleep thresholds in `extras_v2.json` is 0.00321, at
the under-4-hour cut.

**"6 of the 12" external outcomes is an uncorrected count.** The neighbouring sentences say "after correction
for multiple comparisons" and this one does not. Benjamini-Hochberg over the 12 cohort-outcome pairs gives 5.
The 6 is the count the sheet the sentence cites actually draws: `Supp_Fig18_drawn.json`, `rows_b`, the external
rows at `exposure == "t90_adj_ahi"` flagged `significant` on uncorrected P.

**Two counts in the map are rules over rows, not stored cells.** The interquartile bounds `0.67` and `1.07`
of the long-sleep, preserved-oxygen hazard ratios are the 25th and 75th percentiles of the `LN_hr` column
over the 49 primary non-control rows of `dur_x_oxygen_6cell_tst300/results.csv`; the lane stores only the
median and the min-max range. The percentile rule is the linear one (numpy's default), which the printed
`0.67` requires: a plotting-position rule gives `0.66`.

## Regenerating the map

The map is built from a spec that holds, per row, the expression that reads the value out of the file. Both
live outside this repository, with the round that produced them:

```
ROUND55_2026-09-29/working/rows_spec.py          one entry per printed number
ROUND55_2026-09-29/scripts/numbers_map_lib.py    the resolver: file roots, json paths, csv row keys, the
                                                 named count and median rules, the rounding and the verdict
ROUND55_2026-09-29/scripts/build_numbers_map.py  writes NUMBERS_MAP.tsv, taking script and switch from
                                                 RUN_ORDER.md
ROUND55_2026-09-29/scripts/coverage_ledger.py    classifies every printed number of the target text as
                                                 mapped, a repeat, a cross-reference or a citation artefact
ROUND55_2026-09-29/scripts/sidecar_gate.py       checks every output the map reads against the frozen tables'
                                                 shas as they stand now
```

`python3 build_numbers_map.py` rewrites `NUMBERS_MAP.tsv` in place. It reads the frozen outputs and writes
nothing into them.

## Methods parameters

The Methods section prints its own numbers. Those that are settings of the analysis, rather than results,
trace to a line of code rather than to an output table, so they are listed here instead of in the TSV.

Of the 347 numbers the Methods prints, about 271 are cohort and flow counts, Table 1 descriptives, hazard
ratios and cross-references. Those are results, and the ones a reader can check are in the TSV under the
Results section. The 76 below are the settings, each with the line that sets it. Paths beginning `reex/`
expand to `oximetry_extraction/X1_extraction/reex_code/scripts/reex_pipeline/`. The four feature modules
there are byte-identical to the April `signal_features/features/` copies apart from one comment line, and the
September paths are cited because they are what the re-extraction ran.

### Respiratory events and oximetry

| printed | what it sets | file | line |
|---|---|---|---|
| 3% or arousal rule | the hypopnea table per hospital | `oximetry_extraction/X1_extraction/scripts/freeze_hb_definition.py` | 35 |
| events per hour of sleep | the AHI denominator, from decoded staging | `reex/../extract_night.py` | 110 |
| 5, 15, 30 events/h | apnea severity edges | `analysis/numbers/run_all_v2.py` | 70 |
| nadir above the mean | erroneous-oximetry clause 1 | `analysis/numbers/apply_qc_v3.py` | 22 |
| 50 to 80% mean, over 70% below 90, no desaturations | erroneous-oximetry clause 2 | `analysis/numbers/apply_qc_v3.py` | 23 |
| the trace is not resampled | the published oxygen columns are the native-rate reads | `oximetry_extraction/E2_tables/scripts/build_v8_tables.py` | 50 |
| 25, 10, 1, 200 Hz | the oximeter is read at its native rate | `oximetry_extraction/X1_extraction/scripts/extract_v8.py` | 73 |
| below 90% | the T90 threshold | `reex/features/lung_features.py` | 29 |
| 3 percentage points | the ODI3 desaturation depth | `reex/features/lung_features.py` | 34 |
| 4 percentage points | the ODI4 desaturation depth | `oximetry_extraction/X1_extraction/scripts/extract_v8.py` | 87 |
| 30 seconds | the desaturation search window | `reex/features/lung_features.py` | 38 |
| 1 second | the window step | `reex/features/lung_features.py` | 49 |
| 50 to 100% | the saturation artefact-rejection range | `reex/features/lung_features.py` | 22 |
| first percentile | the nadir, on the native trace | `oximetry_extraction/X1_extraction/scripts/extract_v8.py` | 82 |
| first percentile | the published `spo2_nadir_corrected` column is that percentile | `oximetry_extraction/E2_tables/scripts/build_v8_tables.py` | 51 |
| area below the pre-event baseline per hour of sleep | the hypoxic burden | `oximetry_extraction/X1_extraction/scripts/new_measures.py` | 74 |

### Staging, the split-night window and quality

| printed | what it sets | file | line |
|---|---|---|---|
| the pre-treatment window | the split-night cut at the first sustained switch-on | `oximetry_extraction/X1_extraction/scripts/prepap_reader.py` | 73 |
| sleep amounts and composition set missing | the split-night mask column list | `oximetry_extraction/E2_tables/scripts/e2_spec.py` | 64 |
| collapsed to a single stage | the collapsed-staging gate | `oximetry_extraction/X4_groupP/mask_collapsed.py` | 30 |
| once per second | the PAP channel rate | `oximetry_extraction/cpap_stage/extract.py` | 27 |

### Brain spectral power and microstructure

| printed | what it sets | file | line |
|---|---|---|---|
| F3-M2/F4-M1, C3-M2/C4-M1, O1-M2/O2-M1 | the regional derivation pairs | `reex/core/h5_reader.py` | 55 |
| 1 to 4 Hz | delta | `reex/features/brain_features.py` | 39 |
| 4 to 8 Hz | theta | `reex/features/brain_features.py` | 40 |
| 8 to 12 Hz | alpha | `reex/features/brain_features.py` | 41 |
| 12 to 16 Hz | sigma | `reex/features/brain_features.py` | 42 |
| 16 to 30 Hz | beta | `reex/features/brain_features.py` | 43 |
| 30 seconds | the Welch epoch, one segment per epoch, no overlap | `reex/features/brain_features.py` | 68, 54 |
| 10 to 16 Hz | the spindle band-pass | `reex/features/microstructure.py` | 43 |
| two standard deviations | the spindle envelope threshold | `reex/features/microstructure.py` | 61 |
| 0.5 to 3 seconds | the spindle duration limits | `reex/features/microstructure.py` | 70 |
| 0.5 to 1.25 Hz | the slow-oscillation band | `reex/features/microstructure.py` | 110 |
| 0.8 seconds | the minimum trough separation | `reex/features/microstructure.py` | 114 |
| within 1 second | the positive-peak search limit | `reex/features/microstructure.py` | 129 |
| 75 microvolts | the trough-to-peak amplitude floor | `reex/features/microstructure.py` | 131 |
| 0.5 to 1.25 Hz | the Hilbert phase band for coupling | `reex/features/microstructure.py` | 166 |
| within 1 second | the spindle-to-slow-oscillation coupling window | `reex/features/microstructure.py` | 174 |

### Heart rate and variability

| printed | what it sets | file | line |
|---|---|---|---|
| 5 to 15 Hz | the ECG band-pass for R-peak detection | `reex/features/heart_features.py` | 24 |
| 0.4 seconds | the minimum R-peak spacing | `reex/features/heart_features.py` | 27 |
| 300 to 2,000 ms | the R-R interval bounds | `reex/features/heart_features.py` | 37 |
| above 50 ms | the pNN50 threshold | `reex/features/heart_features.py` | 52 |
| 0.04 to 0.15 Hz | the low-frequency band | `reex/features/heart_features.py` | 70 |
| 0.15 to 0.40 Hz | the high-frequency band | `reex/features/heart_features.py` | 71 |
| Lomb-Scargle periodogram | the frequency-domain estimator | `reex/features/heart_features.py` | 69 |

### Outcomes, the cohort and the models

| printed | what it sets | file | line |
|---|---|---|---|
| at least 150 new diagnoses | the outcome inclusion floor | `analysis/numbers/cohort_spec.py` | 175, 194 |
| five negative controls | the control panel | `analysis/numbers/disease_definitions.py` | 185 |
| two controls outside the 52 | contact dermatitis and hemorrhoids excluded from the ranked set | `analysis/numbers/disease_definitions.py` | 202 |
| 141 parameters | the candidate count | `analysis/numbers/cohort_spec.py` | 68 |
| 52 ranked outcomes | the ranked outcome count | `analysis/numbers/cohort_spec.py` | 69 |
| insomnia, restless legs, nocturia, epilepsy left out | the four the sleep study can reveal | `analysis/numbers/build_ranking_v3.py` | 67 |
| fitted in one hospital, scored in the other, both ways | the cross-fit folds | `analysis/numbers/build_ranking_v3.py` | 176 |
| stratified by hospital | the Cox site strata | `analysis/numbers/run_primary_unpenalized.py` | 94 |
| four degrees of freedom | the age natural cubic spline | `analysis/numbers/run_primary_unpenalized.py` | 49 |
| one redundant column dropped | the spline column drop, **in the primary only** | `analysis/numbers/run_primary_unpenalized.py` | 49 |
| rank-based inverse normal | the Blom transform | `analysis/numbers/run_primary_unpenalized.py` | 62 |
| within each hospital | the transform is applied per site | `analysis/numbers/run_primary_unpenalized.py` | 66 |
| 1%, 5%, 10% | the four T90 category edges, `lo < v <= hi` | `analysis/numbers/run_all_v2.py` | 39 |
| the hazard ratio does not fall across the three | the graded rule, part 1 | `analysis/numbers/run_gaps_v2.py` | 50 |
| the top contrast's interval clears 1 | the graded rule, part 2 | `analysis/numbers/run_gaps_v2.py` | 52 |
| at least 10 minutes of scored wake | the wake-and-sleep eligibility floor | `analysis/numbers/run_wake_sleep_healthy.py` | 49 |
| more than 10% of sleep below 90% | the treatment-comparison entry threshold | `analysis/numbers/run_treatment_v3.py` | 71 |
| improved to 10% or less | the responder split | `analysis/numbers/run_treatment_v3.py` | 72 |
| scaled Schoenfeld residual test, rank time transform | the proportional-hazards test | `analysis/numbers/run_ph_all_v8.py` | 107 |
| P below 0.05 | the Schoenfeld violation level | `analysis/numbers/run_ph_all_v8.py` | 45 |
| the first 2 years and after | the time-split windows | `analysis/numbers/run_ph_all_v8.py` | 96 |
| q below 0.05 | the Benjamini-Hochberg criterion | `analysis/numbers/run_wake_sleep_healthy.py` | 50 |
| Benjamini-Hochberg | the step-up implementation | `analysis/numbers/run_treatment_v3.py` | 90 |
| 500 bootstrap replicates | the concordance-gain interval | `analysis/sleep_variability/dC_per_disease/p01_bootstrap.py` | 82 |
| 2.5th and 97.5th percentiles | the interval level and the percentile call | `analysis/sleep_variability/dC_per_disease/p02_assemble.py` | 42, 77 |
| resampled with replacement within each hospital | the bootstrap unit and strata | `analysis/sleep_variability/dC_per_disease/dc_engine.py` | 166 |
| no penalization | the primary per-1-SD models | `analysis/numbers/run_primary_unpenalized.py` | 91 |
| no penalization | the graded-category models | `analysis/numbers/run_all_v2.py` | 114 |
| no penalization | the treatment comparison | `analysis/numbers/run_treatment_v3.py` | 128 |
| a small ridge | 0.01 in the ranking | `analysis/numbers/build_ranking_v3.py` | 182 |
| a small ridge | 0.01 in the wake-and-sleep analysis | `analysis/numbers/run_wake_sleep_healthy.py` | 158 |
| a small ridge | **0.05** in the combination models, not 0.01 | `analysis/numbers/run_combo_v2.py` | 150 |
| Python 3.14 | the interpreter | `requirements.txt`, `README.md` | 1, 52 |

### Where the Methods and the code disagree

Six settings do not read the way the Methods states them. They are listed here rather than buried, because a
reviewer who follows this map will find them.

1. **The external-cohort models are not penalized.** The Methods sentence is "The 141-parameter ranking,
   external-cohort models, wake-and-sleep analysis, combination models, and time-split models used a small
   ridge penalty to stabilize coefficient estimates." `analysis/numbers/freeze_05_external.py` line 58 fits
   `CoxPHFitter(penalizer=0.0)`, its own spec block at line 275 records "penalizer 0", and the frozen output
   `external.json` carries `spec.fit` = "Cox by maximum partial likelihood, penalizer 0".

2. **The time-split models are not penalized either.** `analysis/numbers/run_ph_all_v8.py` line 81 declares
   `split_fit(d, t0, t1, pen=0.0)` and the reported call at line 121 passes no penalty. Lines 154 and 155 run
   0.0 against 0.01 as a side-by-side sensitivity, not as the reported fit. The frozen `ph_all_v8.json`
   records the model as "penalizer 0".

3. **"One redundant column dropped" holds for the primary only.** The 141-parameter ranking keeps all four
   age-spline columns (`build_ranking_v3.py` lines 69 and 119, `dc_engine.py` line 105). The bootstrap that
   supplies the concordance-gain intervals reports the three-column basis while the point estimate it
   brackets comes from the four-column one (`dc_engine.py` line 173).

4. **The ridge is not one value.** 0.01 in the ranking and the wake-and-sleep analysis, 0.05 in the
   combination models (`run_combo_v2.py` line 150), and that file's own note at line 289 says the two are not
   interchangeable.

5. **"No scored desaturation events" is not a scored-event test.** `apply_qc_v3.py` line 23 tests
   `odi3_total == 0`, the pipeline's own 30-second-window desaturation count, not a row count from a scored
   event table.

6. **A superseded primary specification still ships without a guard.** `freeze_03_primary.py` line 86 fits
   the primary associations at `penalizer=0.01`, the specification `run_primary_unpenalized.py` replaced.
   `freeze_04_cpap_stage.py` uses 0.02 and 0.05 for analyses the Methods describes as unpenalized. Unlike
   `build_ranking_v3.py` line 91, which raises on its retired arm, neither carries a marker.

### Settings the code fixes and the Methods does not state

These have a literal in the code and no sentence in the paper. They are listed so a reviewer can find them
rather than infer them: 60 minutes of scored sleep before the AHI and arousal index are defined
(`e2_spec.py` line 75, which is what leaves them undefined on 902 recordings), 5 valid oximetry minutes
(line 76), a 60-minute companion rule for collapsed staging (line 77), the 60 per cent completeness rule and
its companion "more than 10 distinct values" for a candidate measure (`build_ranking_v3.py` line 109), the
event floors of 60 to fit an outcome, 30 per cross-fit fold, 30 per time-split window and 40 for the
treatment family, a floor of 20 per site for the rank transform (`build_ranking_v3.py` line 127), the
100-second hypoxic-burden baseline (`new_measures.py` line 33, with the per-site widths in
`hb_definition.json`, which is data and not in this repository), the sustained-switch-on rule of at least 90
per cent on over the following 30 minutes (`cpap_stage/extract.py` line 143), a Butterworth order of 4 on
every band-pass, a one-minute non-REM floor for spindle detection (`microstructure.py` line 49), and an
R-peak height floor at the 90th percentile of the envelope (`heart_features.py` line 29).

One parameter set could not be traced to a literal. The six implausible-value caps of Methods paragraph 116 (300 and 600
microvolts, 500 ms, two indices at 120 per hour, 13 hours) are applied at `build_v8_tables.py` lines 427 to 437 from
`oximetry_extraction/E2_tables/plausibility_bounds.csv` (column, lo, hi), which ships with the repository and is the
table `e2_spec.BOUNDS_CSV` resolves to when the audit lane's file is absent. The "at least 250 recordings" floor for a
within-group contrast is a reporting rule with no matching literal in the code; a repository-wide search
for 250 returns only the diabetes diagnostic code, two unrelated loop constants and figure geometry.

A comment defect, not a Methods one: the slow-oscillation docstring at `microstructure.py` line 108 says
"period 0.8-2.0 sec" and line 125 says "within 2 sec", while the code at lines 114 and 129 uses 0.8 s and
1.0 s. The code matches the manuscript; the comments do not match the code.
