#!/usr/bin/env bash

# Usage: ./alignments-move2LRA.sh 14444_p1 path/to/LRA hg38 hifi

sample=$1
lra_path=$2
ref=$3
tech=$4

if [[ $ref =~ "hg38" ]]; then
  ref="GRCh38"
else
  ref="CHM13_v2.0"
fi

if [[ $tech =~ "hifi" ]]; then
  tech="PacBio_HiFi"
else
  echo "unknown tech, try again pls" 2>&1
  exit
fi

# -------- Begin -------- #
target_dir="${lra_path}/${sample}/alignments/${tech}/${ref}"
mkdir -p $target_dir
echo "directory made: $target_dir"
find . -name "*${sample}*.bam*" -type f | while read line; do ext=$(echo $line | grep -Eo "bam.*"); cp $line "${target_dir}/alignment.${ext}"; echo "rm $line"; done
echo "files copied, safe to remove"