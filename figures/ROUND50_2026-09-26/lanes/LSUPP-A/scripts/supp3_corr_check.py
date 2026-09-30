#!/usr/bin/env python3
"""In-lane check of the combination_v2.json correlation block (lane V14_L1_RANK, v8.1): which five measures and which order does the matrix carry?
Recomputes Spearman among the five oxygen measures on the v7 cohort exactly as run_combo_v2.py does (within-hospital
rank-normal scores, complete cases on all five) and on the raw scale, then compares the file's matrix to both label orders."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import os, sys, json
os.environ["T90_FULL_NIGHTS_ONLY"] = "1"          # v8.1: run_combo_v2 runs on the 15,551 full diagnostic nights (psv step 116)
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C
OUT = f"{C.LANE}/Supp_Fig03/work/corr_check.json"
T90F = f"{C.DATA_V8}/t90_final.parquet"; MASTER = f"{C.DATA_V8}/master_cohort.csv"
for p in (T90F, MASTER): C.hydrated(p)
sys.path.insert(0, paths.ANALYSIS_DIR)
from cohort_spec import apply_cohort, COHORT_N, DATA_DIR
assert COHORT_N == 15551 and DATA_DIR.endswith("data_frozen_v8_2026-09"), (COHORT_N, DATA_DIR)
NAME = {"spo2_pct_below_90": "Below 90%", "spo2_pct_below_88": "Below 88%", "spo2_nadir_corrected": "Lowest saturation",
        "odi3_total": "3% desaturations", "odi4_total": "4% desaturations", "spo2_mean": "Mean saturation"}
t = pd.read_parquet(T90F)
b = apply_cohort(t); assert len(b) == COHORT_N == 15551
CMB = json.load(open(C.hydrated(f"{C.NUM}/combination_v2.json"))); CMB_SC = C.sidecar_ok(f"{C.NUM}/combination_v2.json")
RK = pd.read_csv(C.hydrated(f"{C.NUM}/ranking_v3.csv"), comment="#")
FAMILY = pd.read_csv(C.hydrated(f"{C.NUM}/measure_families.csv")).set_index("feature")["family"]
# the run_combo_v2 rule: the five oxygen measures with the largest dC in ranking_v3 (TOP5_OXY), in ranking order
FIVE = RK.assign(family=FAMILY.loc[RK.feature].values).query("family == 'Oxygenation'").sort_values("dC", ascending=False).feature.head(5).tolist()
assert len(FIVE) == 5 and all(f in NAME for f in FIVE), FIVE
need = [c for c in FIVE + ["site_id"] if c not in b.columns]
m = pd.read_csv(MASTER, usecols=["BDSPPatientID"] + need, low_memory=False) if need else None
if m is not None:
    b = b.merge(m, on="BDSPPatientID", how="left")
assert b.BDSPPatientID.is_unique
rk_order = RK[RK.feature.isin(FIVE)].sort_values("rank").feature.tolist()      # TOP5_OXY of run_combo_v2.py (ranking order)
print("ranking_v3 order of the five:", rk_order)
Z = {}
for f in FIVE:
    z = pd.Series(np.nan, index=b.index)
    for s, idx in b.groupby("site_id").groups.items():
        v = b.loc[idx, f]; ok = v.notna()
        r = stats.rankdata(v[ok], method="average")
        z.loc[v[ok].index] = stats.norm.ppf((r - 0.375) / (ok.sum() + 0.25))
    Z[f] = z
zf = pd.DataFrame(Z); cc = zf.notna().all(axis=1); n_cc = int(cc.sum())
cm = zf[cc].corr(method="spearman"); raw = b.loc[cc, FIVE].corr(method="spearman")
def reorder(mat, order):
    if not all(o in mat.index for o in order): return np.full((5, 5), np.nan)
    return mat.loc[order, order].values
file_mat = np.array(CMB["correlation_matrix"]); file_raw = np.array(CMB["correlation_matrix_raw_scale"])
file_labels = CMB["correlation_labels"]
inv = {v: k for k, v in NAME.items()}
label_order = [inv[l] for l in file_labels]
file_labels_are_the_five = set(label_order) == set(FIVE)
res = dict(n_complete_cases=n_cc, file_correlation_n=CMB["correlation_n"], file_labels=file_labels, file_labels_are_the_five=file_labels_are_the_five, five=FIVE, ranking_v3_order=rk_order,
           ranking_v3_order_labels=[NAME[f] for f in rk_order], file_label_order_features=label_order,
           max_abs_diff_z_scale=dict(under_ranking_v3_order=float(np.abs(reorder(cm, rk_order) - file_mat).max()),
                                     under_file_label_order=float(np.abs(reorder(cm, label_order) - file_mat).max())),
           max_abs_diff_raw_scale=dict(under_ranking_v3_order=float(np.abs(reorder(raw, rk_order) - file_raw).max()),
                                       under_file_label_order=float(np.abs(reorder(raw, label_order) - file_raw).max())),
           recomputed_z_scale_matrix_ranking_order=[[round(float(x), 3) for x in row] for row in reorder(cm, rk_order)],
           recomputed_raw_scale_matrix_ranking_order=[[round(float(x), 3) for x in row] for row in reorder(raw, rk_order)],
           file_matrix=CMB["correlation_matrix"], mean_abs_corr_recomputed=round(float(np.abs(reorder(cm, rk_order)[np.triu_indices(5, 1)]).mean()), 3),
           mean_abs_corr_file=CMB["mean_abs_correlation"], sidecar=CMB_SC, combination_v2_sha256=C.sha256(f"{C.NUM}/combination_v2.json"),
           recomputed_z_scale_matrix_ranking_order_6dp=[[round(float(x), 6) for x in row] for row in reorder(cm, rk_order)],
           recomputed_raw_scale_matrix_ranking_order_6dp=[[round(float(x), 6) for x in row] for row in reorder(raw, rk_order)],
           printed_2dp_from_file=[[f"{x:.2f}" for x in row] for row in file_mat],
           printed_2dp_from_full_precision=[[f"{float(x):.2f}" for x in row] for row in reorder(cm, rk_order)],
           sign_check="T90 against the lowest saturation must be negative: file cell for the pair under each order",
           t90_vs_nadir_under_ranking_order=float(file_mat[rk_order.index("spo2_pct_below_90"), rk_order.index("spo2_nadir_corrected")]),
           t90_vs_nadir_under_file_labels=(float(file_mat[label_order.index("spo2_pct_below_90"), label_order.index("spo2_nadir_corrected")]) if file_labels_are_the_five else None))
json.dump(res, open(OUT, "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if "matrix" not in k}, indent=1))
print("recomputed z-scale (ranking order):"); print(np.round(reorder(cm, rk_order), 3))
print("file matrix:"); print(file_mat)
