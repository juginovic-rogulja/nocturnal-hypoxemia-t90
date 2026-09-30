#!/bin/bash
# Locations come from the repository's paths.py: run   eval "$(python3 paths.py)"   from the repository root first (it defines the T90_* variables).
# Main_Fig4 round 54 (lane L4, 2026-09-29): the chain. The layout spec (layout_r54.py, every box from the measured strings), the four
# panel builds (a, b delivered on the sheet, c delivered standalone as the counts panel, d for the record only), font-metric alignment,
# the build manifest, the composes (Main_Fig4 = a | b + band 4e as panel c, the counts page, the gs pdfwrite flat twins), the verifier
# (both checks files). Sequential, one step at a time. usage: run_chain_r54.sh
set -u
PY=${T90_PY:-python3}
L=${T90_FIGURE_ROOT}/ROUND54_2026-09-28/lanes/L4
cd "$L/scripts" || exit 1
mkdir -p "$L/logs"; rm -rf "$L/scripts/__pycache__"
for st in layout_r54 01_build_a_r54 01_build_b_r54 01_build_c_r54 01_build_d_r54 01e_align_fonts_r54 01f_build_manifest_r54 02_compose_r54; do
  $PY -u $st.py > "$L/logs/$st.log" 2>&1; rc=$?; echo "$st rc=$rc"; [ $rc -eq 0 ] || { tail -20 "$L/logs/$st.log"; exit 1; }
done
$PY -u 03_verify_r54.py > "$L/logs/03_verify_r54.log" 2>&1; rc=$?; echo "03_verify_r54 rc=$rc"; grep -c "^PASS" "$L/logs/03_verify_r54.log"; grep "^FAIL\|^RESULT\|Traceback\|Error" "$L/logs/03_verify_r54.log" | cut -c1-400
echo "FINAL LINES: Main_Fig4 $(tail -1 "$L/Main_Fig4/verify/checks.txt") | panel_counts_r54 $(tail -1 "$L/panel_counts_r54/verify/checks.txt")"
[ $rc -eq 0 ] && [ "$(tail -1 "$L/Main_Fig4/verify/checks.txt")" = "RESULT ALL PASS" ] && [ "$(tail -1 "$L/panel_counts_r54/verify/checks.txt")" = "RESULT ALL PASS" ]
