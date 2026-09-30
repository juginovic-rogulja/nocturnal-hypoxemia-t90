#!/usr/bin/env bash
# Locations come from the repository's paths.py: run   eval "$(python3 paths.py)"   from the repository root first (it defines the T90_* variables).
# Poll the group P fleets and fetch their outputs. usage: fetch_groupP.sh <mode> <nshard> [max_wait_min]
# Done when FINISHED_<shard> exists for every shard (or PARTIAL/PRECHECK_FAILED markers appear: reported, not fetched as final).
# Every poll syncs the partial shard outputs and logs into X4/<folder>/raw/ so a dead box still leaves its work.
set -uo pipefail
MODE="${1:?mode}"; NSHARD="${2:?nshard}"; MAXMIN="${3:-120}"
BUCKET="${T90_S3_OUTPUT_BUCKET:?set T90_S3_OUTPUT_BUCKET}"; PREFIX="${GROUPP_PREFIX:-outputs/groupP_v8_2}"; REGION=us-east-1
X4="${T90_V8_ROOT}/X4_groupP"
case "$MODE" in cpapA|cpapP|cpapB) DST="$X4/cpap_stage/raw/$MODE";; t90) DST="$X4/t90_by_stage/raw";; oxy) DST="$X4/oxygen_profile/raw";; *) echo "bad mode"; exit 2;; esac
mkdir -p "$DST" "$X4/logs/ec2"
t0=$(date +%s); last=""; still=0
while true; do
  aws s3 sync "s3://$BUCKET/$PREFIX/$MODE/" "$DST/" --region $REGION --only-show-errors --exclude "*" --include "out_*" --include "FINISHED_*" --include "PARTIAL_*" --include "PRECHECK_FAILED_*" --include "ALIVE_*" >/dev/null 2>&1
  aws s3 sync "s3://$BUCKET/$PREFIX/logs/" "$X4/logs/ec2/" --region $REGION --only-show-errors --exclude "*" --include "${MODE}_*.log" >/dev/null 2>&1
  fin=$(ls "$DST" 2>/dev/null | grep -c "^FINISHED_"); bad=$(ls "$DST" 2>/dev/null | grep -cE "^(PARTIAL_|PRECHECK_FAILED_)"); alive=$(ls "$DST" 2>/dev/null | grep -c "^ALIVE_")
  sizes=$(ls -l "$DST"/out_* 2>/dev/null | awk '{s+=$5} END {print s+0}')
  el=$(( ($(date +%s)-t0)/60 ))
  echo "[$(date +%H:%M:%S) +${el}m] $MODE alive=$alive finished=$fin/$NSHARD bad=$bad out_bytes=$sizes"
  if [ "$fin" -ge "$NSHARD" ]; then echo "DONE $MODE"; exit 0; fi
  if [ "$bad" -gt 0 ]; then echo "MARKERS: $(ls "$DST" | grep -E '^(PARTIAL_|PRECHECK_FAILED_)' | tr '\n' ' ')"; fi
  if [ "$sizes" = "$last" ]; then still=$((still+1)); else still=0; fi; last=$sizes
  if [ "$still" -ge 8 ] && [ "$alive" -gt 0 ] && [ "$el" -gt 10 ]; then echo "STALL: no output growth for 8 polls"; fi
  if [ "$el" -ge "$MAXMIN" ]; then echo "TIMEOUT $MODE after $el min ($fin/$NSHARD finished)"; exit 3; fi
  sleep 90
done
