#!/bin/bash
# Locations come from the repository's paths.py: run   eval "$(python3 paths.py)"   from the repository root first (it defines the T90_* variables).
# Main_Fig4 round 49 (lane L4, 2026-09-26): the +1 pt chain. Four panel builds (page-pinned matplotlib on the numbers files, sidecar
# gated), font-metric alignment of the four panel pages (the round-38 compose idiom), the build manifest, the vector compose with the
# coordinator's band_4e_r49.pdf (letters and title 14 pt), the verifier (census, values, neighbours, render). Sequential, one step at a time.
# usage: run_chain_r49.sh [--band PATH]     (default band: ROUND49/lanes/LSKETCH/work/build/band_4e_r49.pdf)
set -u
PY=${T90_PY:-python3}
L=${T90_FIGURE_ROOT}/ROUND49_2026-09-26/lanes/L4
cd "$L/scripts" || exit 1
mkdir -p "$L/logs"
for st in 01_build_a 01_build_b 01_build_c 01_build_d 01e_align_fonts_r49 01f_build_manifest_r49; do
  $PY -u $st.py > "$L/logs/$st.log" 2>&1; rc=$?; echo "$st rc=$rc"; [ $rc -eq 0 ] || { tail -20 "$L/logs/$st.log"; exit 1; }
done
$PY -u 02_compose_r49.py "$@" > "$L/logs/02_compose_r49.log" 2>&1; rc=$?; echo "02_compose_r49 rc=$rc"; tail -1 "$L/logs/02_compose_r49.log" | cut -c1-300; [ $rc -eq 0 ] || exit 1
$PY -u 03_verify_r49.py > "$L/logs/03_verify_r49.log" 2>&1; rc=$?; echo "03_verify_r49 rc=$rc"; grep -c "^PASS" "$L/logs/03_verify_r49.log"; grep "^FAIL\|^RESULT\|Traceback\|Error" "$L/logs/03_verify_r49.log" | cut -c1-400
echo "FINAL LINE: $(tail -1 "$L/Main_Fig4/verify/checks.txt")"
[ $rc -eq 0 ] && [ "$(tail -1 "$L/Main_Fig4/verify/checks.txt")" = "RESULT ALL PASS" ]
