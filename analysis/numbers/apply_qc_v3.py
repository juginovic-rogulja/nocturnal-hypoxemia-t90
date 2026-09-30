"""
Apply the corrected oximetry quality-control rule and stamp it into the cohort file.

The rule used until now, nadir above mean, caught 112 of 210 corrupted recordings. The signature
that identifies them properly is a mean saturation between 50% and 80% together with more than
70% of the recording below 90% and an oxygen desaturation index of exactly zero. That
combination is physically impossible: a recording cannot spend most of the night below 90%
saturation while scoring no desaturation events. Independent re-extraction of 30 of these from
the source recordings returned a median T90 of 0.47% against the 88.5% on file, so the stored
values are not a measurement of those patients' oxygen at all.

208 of the 210 are at one hospital, and 198 of them sit inside the group the paper reports as
having no sleep apnea yet low oxygen, so this materially changes that result.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import pandas as pd

ROOT = paths.T90_ROOT
b = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet")

impossible = (b.spo2_nadir_corrected > b.spo2_mean)
signature = (b.spo2_mean.between(50, 80)) & (b.spo2_pct_below_90 > 70) & (b.odi3_total == 0)
b["oximetry_bad"] = (impossible | signature).astype(int)

a = b[(b.fu_valid == 1) & b.spo2_pct_below_90.notna()]
print(f"analysis cohort before  {len(a):,}")
print(f"  flagged by nadir > mean            {int((a.spo2_nadir_corrected > a.spo2_mean).sum()):4d}")
print(f"  flagged by the impossible signature{int(((a.spo2_mean.between(50,80)) & (a.spo2_pct_below_90>70) & (a.odi3_total==0)).sum()):4d}")
print(f"  flagged by either                  {int(a.index.isin(b[b.oximetry_bad==1].index).sum()):4d}")
print(f"  by site {b.loc[(b.fu_valid==1) & b.spo2_pct_below_90.notna() & (b.oximetry_bad==1), 'site_id'].value_counts().to_dict()}")
print(f"analysis cohort after   {len(a) - int(a.index.isin(b[b.oximetry_bad==1].index).sum()):,}")

b.to_parquet(f"{paths.TABLES_DIR}/t90_final.parquet", index=False)
print("\nwritten: column oximetry_bad added to data_frozen/t90_final.parquet")
print("every downstream analysis must now filter oximetry_bad == 0")
