#!/bin/bash
# Locations come from the repository's paths.py: run   eval "$(python3 paths.py)"   from the repository root first (it defines the T90_* variables).
# T90 v8 extraction worker (self-terminating). No secrets in this file: credentials come from the EC2BDSPAccess instance
# profile. Placeholders __RUN_TAG__ __SHARD__ __N_WORKERS__ __DEADMAN_SEC__ __MANIFEST_KEY__ __TMPFS_GB__ __MODE__ (full|light)
# __NIGHT_TIMEOUT__ are filled by launch_v8.sh.
set -uo pipefail
exec > >(tee -a /var/log/reex_run.log) 2>&1
echo "=== V8 START $(date -u) tag=__RUN_TAG__ mode=__MODE__ shard=__SHARD__ ==="
export AWS_DEFAULT_REGION=us-east-1
export DEBIAN_FRONTEND=noninteractive
BUCKET=${T90_S3_OUTPUT_BUCKET:?set T90_S3_OUTPUT_BUCKET}
PREFIX=outputs/v8_2026-09/__RUN_TAG__
CODE_KEY=code/v8/reex_code_v8.tar.gz
HB_KEY=inputs/v8/hb_definition.json
RUN_TAG=__RUN_TAG__
MODE=__MODE__
TOKEN=$(curl -s --max-time 5 -X PUT "http://169.254.169.254/latest/api/token" -H "X-aws-ec2-metadata-token-ttl-seconds: 21600" || true)
IID=$(curl -s --max-time 5 -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id || echo unknown)
echo "instance $IID"

self_terminate() {
  echo "=== self-terminate $(date -u) reason=$1 ==="
  command -v aws >/dev/null && aws s3 cp /opt/reex/out/${RUN_TAG}.jsonl s3://$BUCKET/$PREFIX/${RUN_TAG}_LAST.jsonl --only-show-errors || true
  command -v aws >/dev/null && aws s3 cp /var/log/reex_run.log s3://$BUCKET/$PREFIX/${RUN_TAG}_run.log --only-show-errors || true
  echo "$1" > /tmp/term.txt; command -v aws >/dev/null && aws s3 cp /tmp/term.txt s3://$BUCKET/$PREFIX/${RUN_TAG}_TERMINATED.txt --only-show-errors || true
  command -v aws >/dev/null && aws ec2 terminate-instances --instance-ids "$IID" --region us-east-1 || true
  sleep 5; shutdown -h now
}
# dead-man timer: terminate no matter what after __DEADMAN_SEC__ seconds
( sleep __DEADMAN_SEC__ && self_terminate deadman ) &

for i in 1 2 3; do apt-get update -qq && break; sleep 10; done
apt-get install -y -qq python3-pip python3-venv curl unzip > /dev/null 2>&1
apt-get install -y -qq awscli > /dev/null 2>&1 || true
if ! command -v aws >/dev/null; then
  cd /tmp && curl -s -o awscliv2.zip "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" && unzip -q awscliv2.zip && ./aws/install >/dev/null 2>&1
fi
if ! command -v aws >/dev/null; then echo "aws cli missing after both install paths"; self_terminate no_aws_cli; fi
echo "aws: $(aws --version 2>&1)"
aws sts get-caller-identity || echo "sts failed"

# heartbeat every 30 s
( while true; do
    N=$(cat /opt/reex/out/${RUN_TAG}.jsonl 2>/dev/null | wc -l)
    printf "%s iid=%s tag=%s mode=%s nights_done=%s load=%s mem_used_mb=%s work_mb=%s\n" \
      "$(date -u +%FT%TZ)" "$IID" "$RUN_TAG" "$MODE" "$N" "$(cut -d' ' -f1-3 /proc/loadavg)" \
      "$(free -m | awk '/Mem:/{print $3}')" "$(du -sm /opt/reex/work 2>/dev/null | cut -f1)" > /tmp/hb.txt
    aws s3 cp /tmp/hb.txt s3://$BUCKET/$PREFIX/heartbeat_${RUN_TAG}.txt --only-show-errors || true
    sleep 30
  done ) &

mkdir -p /opt/reex/work /opt/reex/out
mount -t tmpfs -o size=__TMPFS_GB__G tmpfs /opt/reex/work || echo "tmpfs mount failed, using disk"
python3 -m venv /opt/reex/venv || { echo "venv failed"; self_terminate venv_failed; }
/opt/reex/venv/bin/pip install --quiet --upgrade pip
/opt/reex/venv/bin/pip install --quiet numpy scipy h5py pandas pyarrow s3fs || { echo "pip failed"; self_terminate pip_failed; }
cd /opt/reex
aws s3 cp s3://$BUCKET/$CODE_KEY reex_code_v8.tar.gz --only-show-errors || { echo "code download failed"; self_terminate code_missing; }
tar xzf reex_code_v8.tar.gz
sha256sum reex_code_v8.tar.gz
aws s3 cp s3://$BUCKET/__MANIFEST_KEY__ manifest.csv --only-show-errors || { echo "manifest download failed"; self_terminate manifest_missing; }
aws s3 cp s3://$BUCKET/$HB_KEY hb_definition.json --only-show-errors || { echo "hb_definition download failed"; self_terminate hbdef_missing; }
/opt/reex/venv/bin/python -c "import numpy, scipy, h5py, pandas, s3fs; print('numpy', numpy.__version__, 'scipy', scipy.__version__, 'h5py', h5py.__version__, 'pandas', pandas.__version__, 's3fs', s3fs.__version__)" | tee /opt/reex/out/versions.txt
aws s3 cp /opt/reex/out/versions.txt s3://$BUCKET/$PREFIX/${RUN_TAG}_versions.txt --only-show-errors || true
aws s3 ls "s3://${BDSP_S3_ACCESS_POINT:?set BDSP_S3_ACCESS_POINT}/PSG/bids/I0002/" --page-size 5 2>&1 | head -3

export OUTPUT_BUCKET=$BUCKET OUTPUT_PREFIX=$PREFIX RUN_TAG=$RUN_TAG N_WORKERS=__N_WORKERS__ MANIFEST=/opt/reex/manifest.csv SHARD="__SHARD__" \
       WORK_DIR=/opt/reex/work OUT_DIR=/opt/reex/out NIGHT_TIMEOUT=__NIGHT_TIMEOUT__ HB_DEFINITION=/opt/reex/hb_definition.json \
       OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [ "$MODE" = "light" ]; then
  # periodic partial uploads for the pool-based light pass
  ( while true; do sleep 600; aws s3 cp /opt/reex/out/${RUN_TAG}.jsonl s3://$BUCKET/$PREFIX/${RUN_TAG}_partial.jsonl --only-show-errors || true; done ) &
  /opt/reex/venv/bin/python /opt/reex/scripts/extract_v8.py --manifest /opt/reex/manifest.csv --families light --t-end t_end_sec \
      --source s3fs --workers __N_WORKERS__ --night-timeout __NIGHT_TIMEOUT__ --report-every 200 \
      --hb-definition /opt/reex/hb_definition.json --out /opt/reex/out/${RUN_TAG}.jsonl 2>&1 | tee /var/log/reex_worker.log
  aws s3 cp /opt/reex/out/${RUN_TAG}.jsonl s3://$BUCKET/$PREFIX/${RUN_TAG}_FINAL.jsonl --only-show-errors || true
else
  /opt/reex/venv/bin/python /opt/reex/ec2/worker_v8.py 2>&1 | tee /var/log/reex_worker.log
fi
aws s3 cp /var/log/reex_worker.log s3://$BUCKET/$PREFIX/${RUN_TAG}_worker.log --only-show-errors || true
self_terminate done
