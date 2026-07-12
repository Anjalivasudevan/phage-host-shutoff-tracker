#!/usr/bin/env bash
# Populate config/samples.tsv srr_accession column with real SRR IDs.
#
# The repo ships with PLACEHOLDER values in samples.tsv because accession
# numbers weren't verifiable from the sandbox this repo was scaffolded in
# (no NCBI network access there). Run this once, on a machine with internet
# access, before running the Snakemake pipeline.
#
# Requires NCBI edirect (https://www.ncbi.nlm.nih.gov/books/NBK179288/):
#   sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
#
# Usage: ./scripts/fetch_run_accessions.sh PRJNA929253 > resources/runinfo.csv
# Then manually cross-reference resources/runinfo.csv against the
# timepoint/replicate metadata (sample title / library name columns) and
# fill in config/samples.tsv by hand - automated matching is deliberately
# NOT done here because GEO sample titles/timepoint labels vary in format
# by submission and a wrong automatic match would silently corrupt the
# time course.

set -euo pipefail

BIOPROJECT="${1:?Usage: $0 <BioProject accession, e.g. PRJNA929253>}"

esearch -db sra -query "${BIOPROJECT}" \
  | efetch -format runinfo

echo "Wrote/streamed RunInfo CSV for ${BIOPROJECT}." >&2
echo "Match Run, SampleName, and LibraryName columns against GEO GSM titles" >&2
echo "to fill in config/samples.tsv timepoint/replicate assignments." >&2
