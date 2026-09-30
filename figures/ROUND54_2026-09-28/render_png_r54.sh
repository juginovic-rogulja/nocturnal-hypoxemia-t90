#!/bin/bash
# Locations come from the repository's paths.py: run   eval "$(python3 paths.py)"   from the repository root first (it defines the T90_* variables).
# Round 54 renders: 400-dpi PNGs of the six main sheets (for the manuscript pictures), 300-dpi ED and Supp copies from V30 (unchanged sheets)
R52=${T90_FIGURE_ROOT}/ROUND52_2026-09-26; R54=${T90_FIGURE_ROOT}/ROUND54_2026-09-28; GS=${T90_GS:-gs}; SET="$R54/figures/NEW_FINAL_SET_V31"
mkdir -p "$R54/figures/hires_png_v31" "$R54/figures/png300_v31"
for i in 1 2 3 4 5 6; do $GS -q -dSAFER -dBATCH -dNOPAUSE -dFirstPage=1 -dLastPage=1 -sDEVICE=png16m -r400 -sOutputFile="$R54/figures/hires_png_v31/Main_Fig$i.png" "$SET/Main_Fig$i.pdf"; done
cp -p "$R52/figures/png300_v30/"*.png "$R54/figures/png300_v31/"
$GS -q -dSAFER -dBATCH -dNOPAUSE -dFirstPage=1 -dLastPage=1 -sDEVICE=png16m -r300 -sOutputFile="$R54/figures/png300_v31/ED_Fig07.png" "$SET/ED_Fig07.pdf"   # ED 7 changed in round 54 (panel c)
ls "$R54/figures/hires_png_v31" | wc -l; ls "$R54/figures/png300_v31" | wc -l
