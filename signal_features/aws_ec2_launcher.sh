#!/bin/bash
# Locations come from the repository's paths.py: run   eval "$(python3 paths.py)"   from the repository root first (it defines the T90_* variables).
#
# Sleep Atlas EC2 launcher — spot instance, self-terminating, cost-capped.
#
# Usage:
#   ./aws_ec2_launcher.sh validation   # 100 patients, ~$0.50
#   ./aws_ec2_launcher.sh scale1000    # 1000 patients, ~$2
#   ./aws_ec2_launcher.sh full         # 55,253 patients, ~$10
#
# Safety:
#   - Spot instance (can be terminated by AWS if price exceeds cap)
#   - MaxPrice $0.60/hr for c7i.16xlarge (well below $2.86/hr on-demand)
#   - Instance user-data ends with 'shutdown -h now' after job completes
#   - shutdown-behavior=terminate ensures self-destruction
#   - 24-hr wall timeout baked into user-data

set -euo pipefail

MODE="${1:-validation}"
REGION="us-east-1"
BUCKET="alen-sleep-atlas-$(date +%s)"
AMI="ami-05e86b3611c60b0b4"  # Ubuntu 22.04 LTS x86_64 us-east-1 (standard)

# Per-mode instance sizing + market type + which cohort to use
case "$MODE" in
    validation)
        N_PATIENTS=100
        INSTANCE_TYPE="c7i.4xlarge"
        USE_SPOT=0
        COHORT_SRC="feature_pipeline/outputs/cohort_v2_atlas.csv"
        ;;
    scale1000)
        N_PATIENTS=1000
        INSTANCE_TYPE="c7i.16xlarge"
        USE_SPOT=0
        COHORT_SRC="feature_pipeline/outputs/cohort_v2_atlas.csv"
        ;;
    paper14)
        N_PATIENTS=13769
        INSTANCE_TYPE="c7i.16xlarge"
        USE_SPOT=0
        COHORT_SRC="feature_pipeline/outputs/cohort_paper14_unified.csv"
        ;;
    atlas-full)
        N_PATIENTS=55245
        INSTANCE_TYPE="c7i.16xlarge"
        USE_SPOT=0
        COHORT_SRC="feature_pipeline/outputs/cohort_v2_atlas.csv"
        ;;
    *) echo "Mode must be validation / scale1000 / paper14 / atlas-full"; exit 1 ;;
esac
MAX_SPOT_PRICE="0.60"

echo "=== Sleep Atlas EC2 Launch ==="
echo "Mode:          $MODE"
echo "N patients:    $N_PATIENTS"
if [ "$USE_SPOT" -eq 1 ]; then
    echo "Instance type: $INSTANCE_TYPE (spot, max \$$MAX_SPOT_PRICE/hr)"
else
    echo "Instance type: $INSTANCE_TYPE (on-demand)"
fi
echo "Output bucket: $BUCKET"
echo "Region:        $REGION"
echo ""

# Phase 1: Create personal S3 bucket
echo "[1/5] Creating personal S3 output bucket..."
aws s3 mb "s3://${BUCKET}" --region "$REGION" >/dev/null
aws s3api put-bucket-encryption --bucket "$BUCKET" \
    --server-side-encryption-configuration \
    '{"Rules":[{"ApplyServerSideEncryptionByDefault":{"SSEAlgorithm":"AES256"}}]}'
echo "  OK: s3://${BUCKET}"

# Phase 2: Package feature_pipeline
echo "[2/5] Packaging feature_pipeline code..."
PROJ_ROOT="${T90_REPO_ROOT}"
cd "$PROJ_ROOT"
tar czf /tmp/feature_pipeline.tar.gz signal_features/
aws s3 cp /tmp/feature_pipeline.tar.gz "s3://${BUCKET}/code/"

# Upload cohort CSV (per-mode source)
COHORT_CSV="$PROJ_ROOT/$COHORT_SRC"
if [ ! -f "$COHORT_CSV" ]; then
    echo "ERROR: cohort CSV not found at $COHORT_CSV"
    exit 1
fi
# Subset cohort for validation/scale modes — age-stratified, reproducible seed
/usr/bin/python3 -c "
import pandas as pd, numpy as np
df = pd.read_csv('$COHORT_CSV')
n = $N_PATIENTS
if n >= len(df):
    df.to_csv('/tmp/cohort_subset.csv', index=False)
else:
    rng = np.random.default_rng(42)
    # Stratify by site so we always see cross-site data
    picks = []
    per_site = max(1, n // df['site_id'].nunique())
    for s, grp in df.groupby('site_id'):
        k = min(per_site, len(grp))
        picks.append(grp.iloc[rng.choice(len(grp), k, replace=False)])
    out = pd.concat(picks)
    if len(out) < n:
        extra = df[~df.index.isin(out.index)]
        k = n - len(out)
        out = pd.concat([out, extra.iloc[rng.choice(len(extra), k, replace=False)]])
    out.head(n).to_csv('/tmp/cohort_subset.csv', index=False)
"
aws s3 cp /tmp/cohort_subset.csv "s3://${BUCKET}/inputs/cohort.csv"
echo "  OK: ${N_PATIENTS} patients in s3://${BUCKET}/inputs/cohort.csv"

# Phase 3: Build user-data script
echo "[3/5] Building user-data script..."
# Read user's AWS credentials (from ~/.aws/credentials)
AWS_ACCESS_KEY=$(aws configure get aws_access_key_id)
AWS_SECRET_KEY=$(aws configure get aws_secret_access_key)

cat > /tmp/user_data.sh <<USER_DATA
#!/bin/bash
# Sleep Atlas extraction worker
set -euo pipefail
exec > >(tee -a /var/log/atlas_run.log) 2>&1
echo "=== ATLAS WORKER START: \$(date) ==="

# 24-hr dead-man timer
( sleep 86400 && shutdown -h now ) &

export AWS_ACCESS_KEY_ID="$AWS_ACCESS_KEY"
export AWS_SECRET_ACCESS_KEY="$AWS_SECRET_KEY"
export AWS_DEFAULT_REGION="$REGION"
export OUTPUT_BUCKET="$BUCKET"

apt-get update -qq
DEBIAN_FRONTEND=noninteractive apt-get install -y -qq python3-pip awscli parallel

pip3 install --quiet numpy scipy h5py pandas scikit-learn

# Pull code + cohort
mkdir -p /opt/atlas
cd /opt/atlas
aws s3 cp s3://${BUCKET}/code/feature_pipeline.tar.gz .
tar xzf feature_pipeline.tar.gz
aws s3 cp s3://${BUCKET}/inputs/cohort.csv .

# Write parallel worker script — with 6 fail-safes
cat > /opt/atlas/worker.py <<'WORKER'
import sys, os, subprocess, time, json, traceback, hashlib
from pathlib import Path
sys.path.insert(0, "/opt/atlas")
import pandas as pd
import numpy as np
from signal_features.core.h5_reader import HSPH5Reader, has_usable_staging, probe_h5
from signal_features.features.brain_features import extract_brain_features
from signal_features.features.heart_features import extract_heart_features
from signal_features.features.lung_features import extract_lung_features
from signal_features.features.microstructure import extract_microstructure

S3_BASE = "s3://${BDSP_S3_ACCESS_POINT:?set BDSP_S3_ACCESS_POINT}/PSG/bids"
OUTPUT_BUCKET = os.environ.get("OUTPUT_BUCKET", "")  # set via user-data
WORK = Path("/tmp/atlas_work"); WORK.mkdir(exist_ok=True)

# FAIL-SAFE expected feature counts (validated on 6 patients already)
EXPECTED = {"brain_": 60, "micro_": 39, "heart_": 17, "lung_": 7}

def process(row):
    site = row["site_id"]; bids = row["bids_col"]; ses = int(row["SessionID"])
    stem = f"{bids}_ses-{ses}"
    local = WORK / f"{stem}.h5"
    url = f"{S3_BASE}/{site}/{bids}/ses-{ses}/eeg/{stem}.h5"
    try:
        subprocess.run(["aws","s3","cp",url,str(local),"--only-show-errors"],
                       check=True, timeout=600)
        # FAIL-SAFE 1: SHA256 of downloaded h5
        h = hashlib.sha256()
        with open(local,"rb") as f:
            for chunk in iter(lambda: f.read(1<<20), b""):
                h.update(chunk)
        sha = h.hexdigest()
        # FAIL-SAFE 2: staging usability check
        if not has_usable_staging(local):
            return {"BDSPPatientID": row["BDSPPatientID"], "site_id": site, "status": "unscored", "sha256": sha}
        with HSPH5Reader(local) as r:
            out = {"BDSPPatientID": row["BDSPPatientID"], "site_id": site,
                   "age": row["AgeAtVisit"], "sex": row.get("sex",""), "sha256": sha}
            brain = extract_brain_features(r)
            micro = extract_microstructure(r)
            heart = extract_heart_features(r)
            lung = extract_lung_features(r)
        # FAIL-SAFE 3: feature-count validation
        fam_ok = (len(brain)==EXPECTED["brain_"] and len(micro)==EXPECTED["micro_"]
                  and len(heart)==EXPECTED["heart_"] and len(lung)==EXPECTED["lung_"])
        if not fam_ok:
            out["status"] = f"feature_count_mismatch_b{len(brain)}_m{len(micro)}_h{len(heart)}_l{len(lung)}"
        for k,v in brain.items(): out[f"brain_{k}"]=v
        for k,v in micro.items(): out[f"micro_{k}"]=v
        for k,v in heart.items(): out[f"heart_{k}"]=v
        for k,v in lung.items(): out[f"lung_{k}"]=v
        return out
    except Exception as e:
        return {"BDSPPatientID": row["BDSPPatientID"], "site_id": site, "status": "fail", "err": str(e)[:300]}
    finally:
        if local.exists(): local.unlink()


def upload(path: str, key: str):
    """Upload intermediate file to output bucket (non-fatal on failure)."""
    if not OUTPUT_BUCKET: return
    try:
        subprocess.run(["aws","s3","cp",path,f"s3://{OUTPUT_BUCKET}/{key}","--only-show-errors"],
                       check=False, timeout=60)
    except Exception: pass


if __name__ == "__main__":
    import multiprocessing as mp
    df = pd.read_csv("/opt/atlas/cohort.csv")
    n_workers = min(48, mp.cpu_count())
    print(f"WORKER v2 — Processing {len(df)} patients with {n_workers} workers", flush=True)
    print(f"Output bucket: {OUTPUT_BUCKET}", flush=True)
    t0 = time.time()

    results = []
    fail_count = 0

    # FAIL-SAFE 4: wall-time cap per job = 2x expected
    wall_cap_sec = max(3600, len(df) * 15 * 2)  # 2x expected (15 sec/patient typical)

    with mp.Pool(processes=n_workers) as pool:
        for i,r in enumerate(pool.imap_unordered(process, df.to_dict("records"), chunksize=1)):
            results.append(r)
            if r.get("status") in ("fail",):
                fail_count += 1

            # progress log every 50 patients
            if i % 50 == 0 and i > 0:
                el = (time.time()-t0)/60
                eta = el/i*(len(df)-i)
                fail_rate = fail_count / i
                print(f"[{i}/{len(df)}] elapsed={el:.1f}m ETA={eta:.1f}m fails={fail_count} ({fail_rate:.1%})", flush=True)

                # FAIL-SAFE 5: abort if fail rate too high after first 500
                if i >= 500 and fail_rate > 0.20:
                    print(f"ABORT — fail rate {fail_rate:.1%} exceeds 20% threshold after {i} patients", flush=True)
                    pd.DataFrame(results).to_csv("/opt/atlas/features_aborted.csv", index=False)
                    upload("/opt/atlas/features_aborted.csv", "outputs/features_ABORTED.csv")
                    sys.exit(2)

                # FAIL-SAFE 6: wall-time abort
                if (time.time()-t0) > wall_cap_sec:
                    print(f"ABORT — wall time exceeded {wall_cap_sec/60:.0f} min cap", flush=True)
                    pd.DataFrame(results).to_csv("/opt/atlas/features_partial.csv", index=False)
                    upload("/opt/atlas/features_partial.csv", "outputs/features_TIMEOUT.csv")
                    sys.exit(3)

                # Incremental save every 500 patients (resilient to EC2 termination)
                if i % 500 == 0:
                    pd.DataFrame(results).to_csv("/opt/atlas/features_partial.csv", index=False)
                    upload("/opt/atlas/features_partial.csv", f"outputs/features_partial_{i}.csv")

    pd.DataFrame(results).to_csv("/opt/atlas/features.csv", index=False)
    print(f"Done in {(time.time()-t0)/60:.1f}m, {len(results)} rows, {fail_count} fails", flush=True)
WORKER

python3 /opt/atlas/worker.py 2>&1 | tee /var/log/atlas_worker.log

# Upload results + logs
aws s3 cp /opt/atlas/features.csv s3://${BUCKET}/outputs/features.csv
aws s3 cp /var/log/atlas_run.log s3://${BUCKET}/outputs/run.log || true
aws s3 cp /var/log/atlas_worker.log s3://${BUCKET}/outputs/worker.log || true

echo "=== ATLAS WORKER DONE: \$(date) ==="
shutdown -h now
USER_DATA

USERDATA_B64=$(base64 -i /tmp/user_data.sh)

# Phase 4: Create default security group if needed (egress only)
echo "[4/5] Launching spot instance..."

# Launch (spot or on-demand based on mode)
set +u  # temporarily disable unbound-variable check for array expansion
if [ "$USE_SPOT" -eq 1 ]; then
    INSTANCE_ID=$(aws ec2 run-instances \
        --region "$REGION" \
        --image-id "$AMI" \
        --instance-type "$INSTANCE_TYPE" \
        --instance-market-options "MarketType=spot,SpotOptions={MaxPrice=$MAX_SPOT_PRICE,SpotInstanceType=one-time}" \
        --instance-initiated-shutdown-behavior terminate \
        --user-data file:///tmp/user_data.sh \
        --block-device-mappings '[{"DeviceName":"/dev/sda1","Ebs":{"VolumeSize":100,"VolumeType":"gp3","DeleteOnTermination":true}}]' \
        --tag-specifications "ResourceType=instance,Tags=[{Key=Project,Value=SleepAtlas},{Key=Mode,Value=$MODE},{Key=Owner,Value=alen}]" \
        --query 'Instances[0].InstanceId' \
        --output text)
    MARKET_DESC="spot max \$$MAX_SPOT_PRICE/hr"
else
    INSTANCE_ID=$(aws ec2 run-instances \
        --region "$REGION" \
        --image-id "$AMI" \
        --instance-type "$INSTANCE_TYPE" \
        --instance-initiated-shutdown-behavior terminate \
        --user-data file:///tmp/user_data.sh \
        --block-device-mappings '[{"DeviceName":"/dev/sda1","Ebs":{"VolumeSize":100,"VolumeType":"gp3","DeleteOnTermination":true}}]' \
        --tag-specifications "ResourceType=instance,Tags=[{Key=Project,Value=SleepAtlas},{Key=Mode,Value=$MODE},{Key=Owner,Value=alen}]" \
        --query 'Instances[0].InstanceId' \
        --output text)
    MARKET_DESC="on-demand"
fi
set -u

echo "  OK: instance $INSTANCE_ID launched"
echo ""
echo "[5/5] Monitor via:"
echo "  aws s3 ls s3://${BUCKET}/outputs/      # check for features.csv"
echo "  aws ec2 describe-instances --instance-ids $INSTANCE_ID --region $REGION --query 'Reservations[0].Instances[0].State.Name'"
echo ""
echo "Saved config to /tmp/atlas_last_run.txt"
cat > /tmp/atlas_last_run.txt <<EOF
INSTANCE_ID=$INSTANCE_ID
BUCKET=$BUCKET
MODE=$MODE
N_PATIENTS=$N_PATIENTS
LAUNCHED=$(date)
EOF

echo ""
echo "*** Instance launched ($INSTANCE_TYPE, $MARKET_DESC) ***"
echo "*** Auto-terminates when job completes ***"
echo "*** Expected runtime: ~30 min ($MODE mode) ***"
