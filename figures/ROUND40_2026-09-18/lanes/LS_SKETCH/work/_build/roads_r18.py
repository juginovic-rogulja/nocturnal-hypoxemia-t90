"""The exact wording printed on Figure 6's two dashed roads, and the audit that fixed it.

Alen's round-18 instruction: VERIFY BOTH CLAIMS BEFORE PRINTING THEM, and if any single
outcome rises significantly, say so exactly instead of printing a blanket claim.

WHAT WAS CHECKED
----------------
Duration road, from the duration-by-oxygen cross-classification
  $T90_SV_ROOT/
      dur_x_oxygen_6cell_tst300/results.csv
  45 outcomes plus 5 negative controls, reference = 6 to 7 h of sleep with T90 of 1% or less.
  Cells checked: SN (short sleep, T90 of 1% or less) and LN (long sleep, T90 of 1% or less).
    SN  1 of 45 at raw p, and it runs the protective way: death from any cause
        HR 0.724 (0.567 to 0.926), p = 0.0099, q = 0.447. ZERO raised.
    LN  1 of 45 at raw p and raised: gout HR 2.039 (1.222 to 3.402), p = 0.0064, q = 0.289.
    After Benjamini-Hochberg: ZERO outcomes in either cell. Gout does not clear correction,
    and it is one hit where 2.2 are expected by chance across 45 tests.

Apnea road, from the twelve apnea-by-oxygen cells
  $T90_NUMBERS_DIR/alenfig3_cross_v1.json
  9 outcomes plus 5 negative controls, reference = no apnea with T90 of 1% or less.
  Cells checked: mild, moderate and severe apnea, each with T90 of 1% or less.
    mild      0 of 9 raised at raw p (its back-pain negative control is raised, p = 0.023)
    moderate  0 of 9 raised, the cleanest cell in the file
    severe    1 of 9 raised at raw p: type 2 diabetes HR 1.383 (1.065 to 1.796), p = 0.0151
    After Benjamini-Hochberg: ZERO outcomes in any of the three cells. Type 2 diabetes sits
    at q = 0.136, and the drawn Fig 4d cells print 0 of 4 raised for all three.

THE STANDARD, AND WHY THE SECOND LINE EXISTS
--------------------------------------------
The manuscript's Methods fix the criterion: "P values were corrected using the Benjamini-
Hochberg false discovery rate applied within each analysis across its own set of outcomes,
and a q below 0.05 is the criterion printed in the figures." Both figure legends that carry
these analyses key their asterisks to that same q.

Under that standard both claims are TRUE and the round-16 wording survives. Under raw p BOTH
would break, one on gout and one on type 2 diabetes. A bare "no additional risk" would
therefore be a blanket claim a reviewer can break with either number, so each road carries a
second line naming the standard. That keeps the big line instant and the claim exact, and it
applies the SAME standard to both roads, which a bare claim would not.

"additional" became "added" only to shorten the bold line. Nothing else about the message
changed.
"""
ROAD_DURATION = "short or long sleep alone: no added risk"
ROAD_APNEA = "sleep apnea alone: no added risk"
ROAD_QUALIFIER = "no outcome significant after Benjamini-Hochberg correction"
