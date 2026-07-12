#!/usr/bin/env bash
# Build a combined host+phage FASTA and GTF for STAR alignment, tagging every
# phage contig and gene feature with a prefix so origin is recoverable after
# alignment without a second/competitive mapping step.
#
# Usage:
#   ./scripts/build_reference.sh \
#     resources/genomes/PAO1.fna resources/genomes/PAO1.gtf \
#     resources/genomes/PhiKZ.fna resources/genomes/PhiKZ.gtf \
#     phage_ \
#     resources/genomes/combined_host_phage.fna resources/genomes/combined_host_phage.gtf

set -euo pipefail

HOST_FASTA="${1:?host fasta}"
HOST_GTF="${2:?host gtf}"
PHAGE_FASTA="${3:?phage fasta}"
PHAGE_GTF="${4:?phage gtf}"
PREFIX="${5:?phage prefix, e.g. phage_}"
OUT_FASTA="${6:?output combined fasta}"
OUT_GTF="${7:?output combined gtf}"

mkdir -p "$(dirname "$OUT_FASTA")"

# --- FASTA: prefix phage contig headers, then concatenate ---
TMP_PHAGE_FASTA="$(mktemp)"
awk -v prefix="$PREFIX" '/^>/{sub(/^>/, ">" prefix)} {print}' "$PHAGE_FASTA" > "$TMP_PHAGE_FASTA"
cat "$HOST_FASTA" "$TMP_PHAGE_FASTA" > "$OUT_FASTA"
rm -f "$TMP_PHAGE_FASTA"

# --- GTF: prefix phage seqname (col 1) and gene_id/transcript_id values ---
TMP_PHAGE_GTF="$(mktemp)"
awk -v prefix="$PREFIX" 'BEGIN{FS=OFS="\t"} \
  /^#/ {print; next} \
  { $1 = prefix $1; \
    gsub(/gene_id "/, "gene_id \"" prefix, $9); \
    gsub(/transcript_id "/, "transcript_id \"" prefix, $9); \
    print }' "$PHAGE_GTF" > "$TMP_PHAGE_GTF"
cat "$HOST_GTF" "$TMP_PHAGE_GTF" > "$OUT_GTF"
rm -f "$TMP_PHAGE_GTF"

echo "Combined reference written:" >&2
echo "  FASTA: $OUT_FASTA" >&2
echo "  GTF:   $OUT_GTF" >&2
echo "Phage contigs/genes tagged with prefix: '$PREFIX'" >&2
