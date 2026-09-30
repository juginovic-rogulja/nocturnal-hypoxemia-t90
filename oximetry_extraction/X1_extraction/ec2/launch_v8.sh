#!/bin/zsh
# Locations come from the repository's paths.py: run   eval "$(python3 paths.py)"   from the repository root first (it defines the T90_* variables).
# Launch one self-terminating v8 extraction instance (costs money: run only on the owner's go).
# usage: launch_v8.sh <RUN_TAG> <INSTANCE_TYPE> <MANIFEST_KEY> <SHARD i/n or ''> <N_WORKERS> <DEADMAN_SEC> <TMPFS_GB> <NAME> [MODE full|light] [NIGHT_TIMEOUT]
set -euo pipefail
RUN_TAG=$1; ITYPE=$2; MANIFEST_KEY=$3; SHARD=${4:-}; NW=$5; DEADMAN=$6; TMPFS=$7; NAME=$8; MODE=${9:-full}; NT=${10:-1500}
REGION=us-east-1
AMI=ami-05a3e9423ae4d7a19   # Canonical Ubuntu 22.04 LTS amd64, the v7 AMI (2026-09-06/07 runs)
HERE=$(cd "$(dirname "$0")" && pwd)
UD=$HERE/user_data_${RUN_TAG}.sh
sed -e "s#__RUN_TAG__#${RUN_TAG}#g" -e "s#__SHARD__#${SHARD}#g" -e "s#__N_WORKERS__#${NW}#g" -e "s#__DEADMAN_SEC__#${DEADMAN}#g" \
    -e "s#__MANIFEST_KEY__#${MANIFEST_KEY}#g" -e "s#__TMPFS_GB__#${TMPFS}#g" -e "s#__MODE__#${MODE}#g" -e "s#__NIGHT_TIMEOUT__#${NT}#g" \
    "$HERE/user_data_v8_template.sh" > "$UD"
grep -q "AWS_SECRET\|aws_secret_access_key" "$UD" && { echo "secret in user-data, abort"; exit 1; }
aws s3 ls "s3://${T90_S3_OUTPUT_BUCKET:?set T90_S3_OUTPUT_BUCKET}/code/v8/reex_code_v8.tar.gz" >/dev/null || { echo "code tarball not on S3 (run build_code_tarball.sh --upload)"; exit 1; }
aws s3 ls "s3://${T90_S3_OUTPUT_BUCKET:?set T90_S3_OUTPUT_BUCKET}/${MANIFEST_KEY}" >/dev/null || { echo "manifest ${MANIFEST_KEY} not on S3"; exit 1; }
aws s3 ls "s3://${T90_S3_OUTPUT_BUCKET:?set T90_S3_OUTPUT_BUCKET}/inputs/v8/hb_definition.json" >/dev/null || { echo "hb_definition.json not on S3"; exit 1; }
IID=$(aws ec2 run-instances --region $REGION --image-id $AMI --instance-type $ITYPE \
  --iam-instance-profile Name=EC2BDSPAccess \
  --instance-initiated-shutdown-behavior terminate \
  --user-data "file://$UD" \
  --block-device-mappings '[{"DeviceName":"/dev/sda1","Ebs":{"VolumeSize":40,"VolumeType":"gp3","DeleteOnTermination":true}}]' \
  --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=${NAME}},{Key=Project,Value=T90_v8},{Key=RunTag,Value=${RUN_TAG}},{Key=Owner,Value=alen}]" \
  --query 'Instances[0].InstanceId' --output text)
echo "$(date -u +%FT%TZ) $RUN_TAG $ITYPE $IID mode=$MODE workers=$NW deadman=${DEADMAN}s manifest=$MANIFEST_KEY" >> "$HERE/instances_launched.txt"
echo "$IID"
