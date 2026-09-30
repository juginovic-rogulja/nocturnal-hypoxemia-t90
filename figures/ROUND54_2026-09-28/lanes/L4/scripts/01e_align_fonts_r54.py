#!$T90_PY
"""ROUND 54 (2026-09-29, lane L4, copied from round 49): the four panel pages get the FontDescriptor Ascent/Descent the same Arial face carries on the V13
sheet (the round-28 rule, l4lib.align_font_metrics), work/panel_<x>_raw.pdf -> work/panel_<x>_m.pdf. In round 38 this ran inside
02_compose.py for all four panels (builders a and b also align their own page for the read-back); here it is its own step so the
compose only places pages. Idempotent."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import l4lib as L
WORK = f"{L.LANE}/Main_Fig4/work"
met = L.sheet_font_metrics(json.load(open(f"{WORK}/base_text.json"))["spans"])
for x in "abcd":
    raw, m = f"{WORK}/panel_{x}_raw.pdf", f"{WORK}/panel_{x}_m.pdf"; assert os.path.exists(raw) and L.blocks(raw) > 0, raw
    done = L.align_font_metrics(raw, m, met); print(f"panel {x}: aligned {done} -> {m} ({os.path.getsize(m)} bytes)")
