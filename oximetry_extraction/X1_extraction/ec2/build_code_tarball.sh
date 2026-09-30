#!/bin/zsh
# Locations come from the repository's paths.py: run   eval "$(python3 paths.py)"   from the repository root first (it defines the T90_* variables).
# Build X1/ec2/reex_code_v8.tar.gz (scripts/ + reex_code/ + ec2/worker_v8.py) with SHA256; --upload also puts the tarball,
# the manifests and hb_definition.json on S3 (no cost). Never launches anything.
set -euo pipefail
X1=${T90_V8_ROOT}/X1_extraction; X1_CODE=$(cd "$(dirname "$0")/.." && pwd)   # scripts/, reex_code/ and ec2/ of the repository
BUCKET=${T90_S3_OUTPUT_BUCKET:?set T90_S3_OUTPUT_BUCKET}
STAGE=$(mktemp -d)
mkdir -p $STAGE/scripts $STAGE/ec2 $STAGE/reex_code
cp $X1_CODE/scripts/{x1_paths,provenance,prepap_reader,new_measures,extract_v8}.py $STAGE/scripts/
cp -R $X1_CODE/reex_code/scripts $STAGE/reex_code/
cp $X1_CODE/ec2/worker_v8.py $STAGE/ec2/
find $STAGE -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true
( cd $STAGE && tar czf $X1/ec2/reex_code_v8.tar.gz scripts reex_code ec2 )
shasum -a 256 $X1/ec2/reex_code_v8.tar.gz | tee $X1/ec2/reex_code_v8.tar.gz.sha256
tar tzf $X1/ec2/reex_code_v8.tar.gz | wc -l
rm -rf $STAGE
if [ "${1:-}" = "--upload" ]; then
  aws s3 cp $X1/ec2/reex_code_v8.tar.gz s3://$BUCKET/code/v8/reex_code_v8.tar.gz --only-show-errors
  aws s3 cp $X1/ec2/reex_code_v8.tar.gz.sha256 s3://$BUCKET/code/v8/reex_code_v8.tar.gz.sha256 --only-show-errors
  for f in manifest_split_3622.csv manifest_full_19173_v8.csv manifest_pilot300_split.csv sample_300.csv hb_definition.json; do
    [ -f $X1/work/$f ] && aws s3 cp $X1/work/$f s3://$BUCKET/inputs/v8/$f --only-show-errors && echo "uploaded inputs/v8/$f"
  done
  aws s3 ls s3://$BUCKET/code/v8/; aws s3 ls s3://$BUCKET/inputs/v8/
fi
