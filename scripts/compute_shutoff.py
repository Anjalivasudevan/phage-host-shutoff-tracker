#!/usr/bin/env python3
"""
Parse featureCounts gene-level count tables (one per sample, or one combined
matrix) and compute host vs phage read fractions per sample, using the
phage_prefix from config.yaml to classify each gene by origin.

Usage:
  python scripts/compute_shutoff.py \
    --counts results/counts/*.featureCounts.txt \
    --samples config/samples.tsv \
    --phage-prefix phage_ \
    --out results/shutoff_timecourse.tsv \
    --log results/logs/compute_shutoff.log
"""
import argparse
import glob
import logging
import sys
from pathlib import Path

import pandas as pd


def load_featurecounts(path: str) -> pd.Series:
    """featureCounts output: comment line, then header, then Geneid ... count col (last col)."""
    df = pd.read_csv(path, sep="\t", comment="#")
    if df.shape[1] < 7:
        raise ValueError(f"{path}: unexpected featureCounts format ({df.shape[1]} cols)")
    gene_ids = df["Geneid"]
    counts = df.iloc[:, -1]  # last column is the sample's count column
    return pd.Series(counts.values, index=gene_ids.values)


def sample_id_from_path(path: str) -> str:
    return Path(path).name.split(".")[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", nargs="+", required=True,
                     help="One or more featureCounts output files (glob-expanded by shell or here).")
    ap.add_argument("--samples", required=True, help="config/samples.tsv")
    ap.add_argument("--phage-prefix", default="phage_")
    ap.add_argument("--out", required=True)
    ap.add_argument("--log", default=None)
    args = ap.parse_args()

    handlers = [logging.StreamHandler(sys.stderr)]
    if args.log:
        Path(args.log).parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(args.log))
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", handlers=handlers)
    log = logging.getLogger("compute_shutoff")

    # expand any glob patterns that weren't already expanded by the shell
    count_files = []
    for pattern in args.counts:
        matches = glob.glob(pattern)
        count_files.extend(matches if matches else [pattern])
    if not count_files:
        log.error("No count files found matching: %s", args.counts)
        sys.exit(1)

    samples = pd.read_csv(args.samples, sep="\t")
    samples = samples[samples["srr_accession"] != "PLACEHOLDER_FILL_FROM_RUN_SELECTOR"]
    if samples.empty:
        log.warning("samples.tsv still has only PLACEHOLDER accessions - "
                     "run scripts/fetch_run_accessions.sh and fill it in first.")

    rows = []
    for path in sorted(count_files):
        sid = sample_id_from_path(path)
        counts = load_featurecounts(path)

        is_phage = counts.index.to_series().str.startswith(args.phage_prefix)
        phage_reads = int(counts[is_phage].sum())
        host_reads = int(counts[~is_phage].sum())
        total = phage_reads + host_reads

        n_unmatched = 0  # kept for symmetry; every gene is classified as host or phage by construction
        if total == 0:
            log.warning("%s: zero total counted reads, skipping", sid)
            continue

        row = {
            "sample_id": sid,
            "host_reads": host_reads,
            "phage_reads": phage_reads,
            "total_reads": total,
            "host_fraction": host_reads / total,
            "phage_fraction": phage_reads / total,
        }

        meta = samples[samples["sample_id"] == sid]
        if not meta.empty:
            row["timepoint_min"] = int(meta.iloc[0]["timepoint_min"])
            row["replicate"] = int(meta.iloc[0]["replicate"])
        else:
            log.warning("%s: not found in samples.tsv, timepoint/replicate left blank", sid)
            row["timepoint_min"] = None
            row["replicate"] = None

        rows.append(row)
        log.info("%s: host=%.3f phage=%.3f (n=%d reads)", sid, row["host_fraction"], row["phage_fraction"], total)

    if not rows:
        log.error("No samples produced usable counts - nothing to write.")
        sys.exit(1)

    out_df = pd.DataFrame(rows).sort_values(["timepoint_min", "replicate"])
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, sep="\t", index=False)
    log.info("Wrote %s (%d samples)", args.out, len(out_df))


if __name__ == "__main__":
    main()
