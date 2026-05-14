#!/usr/bin/env bash
# Extract human IGH reference files from immcantation/suite:4.5.0 into refs/
# Run once from the project root: bash scripts/extract_refs.sh

set -euo pipefail

SRC_IMAGE="immcantation/suite:4.5.0"
DEST="$(cd "$(dirname "$0")/.." && pwd)/refs"

echo "Extracting refs from ${SRC_IMAGE} → ${DEST}"

CID=$(docker create "$SRC_IMAGE")
trap "docker rm -f $CID >/dev/null" EXIT

# --- pRESTO: V-gene FASTA for AssemblePairs ---
mkdir -p "${DEST}/igblast/fasta"
docker cp "${CID}:/usr/local/share/igblast/fasta/imgt_human_ig_v.fasta" \
    "${DEST}/igblast/fasta/"

# --- Change-O: IgBLAST databases for AssignGenes (human IG only, incl. C for AssignGenes -c_region_db) ---
mkdir -p "${DEST}/igblast/database"
for segment in v d j c; do
    docker cp "${CID}:/usr/local/share/igblast/database/imgt_human_ig_${segment}.ndb"  "${DEST}/igblast/database/" 2>/dev/null || true
    docker cp "${CID}:/usr/local/share/igblast/database/imgt_human_ig_${segment}.nhr"  "${DEST}/igblast/database/"
    docker cp "${CID}:/usr/local/share/igblast/database/imgt_human_ig_${segment}.nin"  "${DEST}/igblast/database/"
    docker cp "${CID}:/usr/local/share/igblast/database/imgt_human_ig_${segment}.njs"  "${DEST}/igblast/database/" 2>/dev/null || true
    docker cp "${CID}:/usr/local/share/igblast/database/imgt_human_ig_${segment}.nog"  "${DEST}/igblast/database/"
    docker cp "${CID}:/usr/local/share/igblast/database/imgt_human_ig_${segment}.nos"  "${DEST}/igblast/database/" 2>/dev/null || true
    docker cp "${CID}:/usr/local/share/igblast/database/imgt_human_ig_${segment}.not"  "${DEST}/igblast/database/" 2>/dev/null || true
    docker cp "${CID}:/usr/local/share/igblast/database/imgt_human_ig_${segment}.nsq"  "${DEST}/igblast/database/"
    docker cp "${CID}:/usr/local/share/igblast/database/imgt_human_ig_${segment}.ntf"  "${DEST}/igblast/database/" 2>/dev/null || true
    docker cp "${CID}:/usr/local/share/igblast/database/imgt_human_ig_${segment}.nto"  "${DEST}/igblast/database/" 2>/dev/null || true
done

# --- Change-O: IgBLAST internal framework data + aux file ---
mkdir -p "${DEST}/igblast/internal_data/human"
docker cp "${CID}:/usr/local/share/igblast/internal_data/human/." \
    "${DEST}/igblast/internal_data/human/"

mkdir -p "${DEST}/igblast/optional_file"
docker cp "${CID}:/usr/local/share/igblast/optional_file/human_gl.aux" \
    "${DEST}/igblast/optional_file/"

# --- Change-O: IMGT germlines for MakeDb (IGH only) ---
mkdir -p "${DEST}/germlines/imgt/human/vdj"
for gene in IGHV IGHD IGHJ; do
    docker cp "${CID}:/usr/local/share/germlines/imgt/human/vdj/imgt_human_${gene}.fasta" \
        "${DEST}/germlines/imgt/human/vdj/"
done

echo "Done. Refs extracted to ${DEST}"
echo ""
echo "Directory layout:"
find "${DEST}" -type f | sort
