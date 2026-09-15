# phage-host-shutoff

Track bacterial host transcriptional shutoff during lytic bacteriophage infection,
using RNA-seq time-course data.

**Question:** as a lytic phage takes over the cell, what fraction of total mRNA
reads come from the host genome vs. the phage genome, and how fast does that
ratio flip?

## Dataset

Default config targets **GSE223979** (BioProject PRJNA929253):
*Pseudomonas aeruginosa* PAO1 infected with giant phage **PhiKZ**, sampled at
−2, 0, 2, 4, 6, 8, and 10 minutes post-infection, 2 biological replicates per
time point. Small, high-depth, and captures the very earliest window of host
shutdown.

Swappable — `config/config.yaml` also has a commented-out block for
**GSE58494** (PhiKZ, 0–35 min, coarser time course) if you want a longer arc
instead of high early-time resolution.

## Pipeline

```
SRA FASTQ (per timepoint/replicate)
   -> cutadapt (adapter/quality trim)
   -> FastQC (QC report)
   -> STAR align to combined host+phage reference
   -> featureCounts (combined GTF, gene-level)
   -> compute_shutoff.py (host vs phage read fraction per sample)
   -> plot_shutoff.R (stacked-area timecourse + per-replicate line plot)
```

Reference is a single FASTA/GTF built by concatenating the host and phage
genomes, with phage contigs/genes prefixed (`phage_`) so origin can be
recovered post-alignment without a second alignment pass.

## Repo layout

```
config/
  config.yaml          # accessions, paths, params - edit this first
  samples.tsv           # sample sheet (SRR id, timepoint, replicate)
envs/
  environment.yml       # conda env: star, cutadapt, fastqc, subread, snakemake, R/DESeq2 not required here (no DE test), tidyverse
scripts/
  build_reference.sh    # concatenate host+phage fasta/gtf, tag phage features
  compute_shutoff.py    # parse featureCounts output -> host/phage fraction table
  plot_shutoff.R         # stacked-area + line plot of shutoff over time
Snakefile
```

## Setup

```bash
conda env create -f envs/environment.yml
conda activate phage-shutoff

# pull genomes (host PAO1 + PhiKZ) and accession list - see config/config.yaml
# then:
snakemake -j 8 --use-conda
```



## Output

`results/shutoff_timecourse.tsv` - one row per sample: timepoint, replicate,
host_reads, phage_reads, host_fraction, phage_fraction.

`results/plots/shutoff_stacked_area.png` - the headline figure: host fraction
collapsing / phage fraction rising across the time course.

`results/plots/shutoff_by_replicate.png` - same data as separate replicate
lines, to sanity-check reproducibility before trusting the averaged plot.

## Notes / gotchas

- featureCounts is run in gene-level mode against the combined GTF; genes
  with no `phage_`/host prefix match are dropped with a warning rather than
  silently miscounted - check `results/logs/compute_shutoff.log` if your
  fractions look off.
- PhiKZ has an unusually large genome (~280 kb) and its own RNA polymerase
  genes fire very early, so don't be surprised if the phage fraction is
  already non-trivial at the 0 min timepoint (adsorption/injection isn't
  instantaneous across the culture).
- This is fraction-of-reads shutoff tracking, not per-gene differential
  expression. If you want per-gene host shutdown dynamics on top of this,
  that's a natural follow-up repo (DESeq2 time-course, same idea as your
  HNSCC tumor/normal pseudobulk work, just with infection time instead of
  tumor status as the covariate).

## License

MIT
# phage-host-shutoff-tracker
