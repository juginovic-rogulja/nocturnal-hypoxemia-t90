#!/bin/bash
# Locations come from the repository's paths.py: run   eval "$(python3 paths.py)"   from the repository root first (it defines the T90_* variables).
# Round 54 finish: collect the lane checks and reports, park the 19:15 package and submission folders of 26 September (round 52 builds)
# into ROUND52/package_superseded, build the 29 September package and submission folder with the V31 figures (documents unchanged), report copy, mirror.
set -e
T90=${T90_ROOT}; R52=$T90/ROUND52_2026-09-26; R54=$T90/ROUND54_2026-09-28; PY=${T90_PY:-python3}; cd "$R54"
head -1 reports/FINAL_CHECK_R54.md | grep -q "RESULT ALL PASS" || { echo "round-52 final check not ALL PASS"; exit 1; }
test -s reports/ROUND54_REPORT.md || { echo "ROUND54_REPORT.md missing"; exit 1; }
: > reports/LANE_checks_r54.txt; : > reports/LANE_reports_r54.md
for c in lanes/*/*/verify/checks.txt; do echo "== $c"; cat "$c"; echo; done >> reports/LANE_checks_r54.txt
for r in lanes/*/REPORT.md; do echo "# $r"; cat "$r"; echo; done >> reports/LANE_reports_r54.md
mkdir -p "$R52/package_superseded"
OLDP="$HOME/Downloads/T90 Paper Package 2026-09-26"; OLDS="$HOME/Downloads/T90 Submission 2026-09-26"
if [ -d "$OLDP" ]; then test ! -e "$R52/package_superseded/T90 Paper Package 2026-09-26 (19-15 build)" && mv "$OLDP" "$R52/package_superseded/T90 Paper Package 2026-09-26 (19-15 build)"; fi
if [ -d "$OLDS" ]; then test ! -e "$R52/package_superseded/T90 Submission 2026-09-26 (19-15 build)" && mv "$OLDS" "$R52/package_superseded/T90 Submission 2026-09-26 (19-15 build)"; fi
$PY assemble_package_r54.py --run | tail -2
P="$HOME/Downloads/T90 Paper Package 2026-09-29"; S="$HOME/Downloads/T90 Submission 2026-09-29"
test ! -e "$S" || { echo "$S exists"; exit 1; }
mkdir -p "$S/Figures - separate files"; for f in "T90 Manuscript.docx" "T90 Supplement.docx" "T90 Extended Data Figure Legends.docx" "T90 Cover Letter.docx"; do cp -p "$P/clean/$f" "$S/$f"; done
cp -p "$P/figures/T90 All Figures V31.pdf" "$S/T90 Figures.pdf"; cp -p "$P/figures/T90 Figures V31 - separate PDFs/"*.pdf "$S/Figures - separate files/"; cp -p "$P/figures/T90 All Figures V31.docx" "$S/T90 Figures (Word).docx"
echo "submission folder built: $(ls "$S/Figures - separate files" | wc -l | tr -d ' ') figure files"
cp -p reports/ROUND54_REPORT.md "$HOME/Downloads/T90 ROUND 54 REPORT 2026-09-29.md"
rsync -a --exclude logs --exclude work --exclude _diag "$R54/" "${T90_MIRROR_ROOT:-/nonexistent}/ROUND54_2026-09-28/" 2>/dev/null || echo "mirror skipped"
echo "FINISH DONE $(date '+%Y-%m-%d %H:%M')"
