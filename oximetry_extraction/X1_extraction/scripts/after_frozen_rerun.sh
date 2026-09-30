#!/bin/zsh
# Locations come from the repository's paths.py: run   eval "$(python3 paths.py)"   from the repository root first (it defines the T90_* variables).
# Post-rerun chain: check the frozen pilot, compare it with the default-window run, smoke the light assembly on it, and launch
# the 50-night second path with the frozen definition in the background.
set -uo pipefail
X1=${T90_V8_ROOT}/X1_extraction; X1_SCRIPTS=${T90_X1_SCRIPTS}
PY=${T90_PY:-python3}
cd $X1   # data, work and logs live here; the scripts run from the repository
$PY $X1_SCRIPTS/check_pilot_light.py --jsonl data/light_pilot300_frozen.jsonl --report logs/PILOT300_LIGHT_frozen.md --calib work/hb_window_calibration_frozen_run.json > logs/020b_check_frozen.log 2>&1; echo "check_frozen rc=$?"
grep -E 'identical|Spearman|within 1 min|MISMATCHES|hb_window_source|hypoxic_burden:' logs/020b_check_frozen.log | cut -c1-220
$PY $X1_SCRIPTS/compare_pilot_runs.py data/light_pilot300.jsonl data/light_pilot300_frozen.jsonl > logs/020b_compare_runs.log 2>&1; echo "compare rc=$?"; cat logs/020b_compare_runs.log | cut -c1-300
$PY $X1_SCRIPTS/assemble_light.py --jsonl data/light_pilot300_frozen.jsonl --manifest work/sample_300.csv --out data/light_pilot300_assembled.parquet --expect-n 300 --log-suffix _pilot > logs/023_assemble_light_pilot_smoke.log 2>&1; echo "assemble_light smoke rc=$?"; grep -E 'ASSEMBLE_LIGHT|identical|FAIL|D11|windowed' logs/023_assemble_light_pilot_smoke.log | cut -c1-220
nohup $PY $X1_SCRIPTS/second_path_check.py --pass-jsonl data/light_pilot300_frozen.jsonl --manifest work/sample_300.csv --n 50 --hb-definition work/hb_definition.json --out logs/SECOND_PATH_50.csv > logs/024_second_path_frozen.log 2>&1 &
echo "second path (frozen) launched"
