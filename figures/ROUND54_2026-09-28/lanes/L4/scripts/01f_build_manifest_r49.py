#!$T90_PY
"""ROUND 49 (2026-09-26, lane L4): the build manifest work/build_manifest_r49.json (sha256 of every panel page, raw page, drawn record,
geometry record and builder script); the compose gates the panel pages against it."""
import json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import l4lib as L
W = f"{L.LANE}/Main_Fig4/work"; V = f"{L.LANE}/Main_Fig4/verify"
man = {"lane": "L4 round 49", "written": time.strftime("%Y-%m-%d %H:%M:%S"), "panels": {}}
for x in "abcd":
    m = f"{W}/panel_{x}_m.pdf"; raw = f"{W}/panel_{x}_raw.pdf"; dj = f"{V}/Main_Fig4_{x}_drawn.json"; gj = f"{W}/panel_{x}_geometry.json"
    for p in (m, raw, dj, gj): assert os.path.exists(p) and L.blocks(p) > 0, p
    assert os.path.getmtime(m) >= os.path.getmtime(raw) - 1, (x, "the aligned page is older than the raw page")
    man["panels"][x] = {"file": m, "sha256": L.sha256(m), "bytes": os.path.getsize(m), "mtime": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(os.path.getmtime(m))),
                        "raw_sha256": L.sha256(raw), "drawn_json": dj, "drawn_sha256": L.sha256(dj), "geometry_json": gj, "geometry_sha256": L.sha256(gj),
                        "builder": f"{L.LANE}/scripts/01_build_{x}.py", "builder_sha256": L.sha256(f"{L.LANE}/scripts/01_build_{x}.py")}
    print(x, man["panels"][x]["sha256"][:16], man["panels"][x]["bytes"], "bytes")
json.dump(man, open(f"{W}/build_manifest_r49.json", "w"), indent=1); print("wrote", f"{W}/build_manifest_r49.json")
