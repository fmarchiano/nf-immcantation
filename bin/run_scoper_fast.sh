#!/usr/bin/env bash
# Wrapper for scoper-fast matching the run_scoper.R CLI contract.
# Usage: run_scoper_fast.sh <outname> <method> <linkage> <threshold> <nproc> <tsv1> [tsv2 ...]
set -euo pipefail

outname=$1; method=$2; linkage=$3; threshold=$4; nproc=$5
shift 5
tsv_files=("$@")

if [[ "$method" != "novj" ]]; then
    echo "ERROR: scoper-fast only supports method=novj (spectral), got '$method'" >&2
    exit 1
fi

input_args=()
for f in "${tsv_files[@]}"; do
    input_args+=(--input "$f")
done

scoper-fast \
    "${input_args[@]}" \
    --outname "$outname" \
    --method novj \
    --threads "$nproc"

# Nextflow expects versions.yml without the outname prefix
mv "${outname}_versions.yml" versions.yml
