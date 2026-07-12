import pandas as pd

configfile: "config/config.yaml"

samples = pd.read_csv(config["sample_sheet"], sep="\t")
samples = samples[samples["srr_accession"] != "PLACEHOLDER_FILL_FROM_RUN_SELECTOR"]
SAMPLE_IDS = list(samples["sample_id"])

def srr_for(wildcards):
    return samples.loc[samples.sample_id == wildcards.sample, "srr_accession"].iloc[0]

rule all:
    input:
        "results/shutoff_timecourse.tsv",
        "results/plots/shutoff_stacked_area.png",
        "results/plots/shutoff_by_replicate.png"

# ---------------------------------------------------------------------------
# Reference
# ---------------------------------------------------------------------------

rule build_reference:
    input:
        host_fa=config["host_genome_fasta"],
        host_gtf=config["host_genome_gtf"],
        phage_fa=config["phage_genome_fasta"],
        phage_gtf=config["phage_genome_gtf"]
    output:
        fa=config["combined_fasta"],
        gtf=config["combined_gtf"]
    params:
        prefix=config["phage_prefix"]
    shell:
        "bash scripts/build_reference.sh {input.host_fa} {input.host_gtf} "
        "{input.phage_fa} {input.phage_gtf} {params.prefix} {output.fa} {output.gtf}"

rule star_index:
    input:
        fa=config["combined_fasta"],
        gtf=config["combined_gtf"]
    output:
        directory(config["star_index_dir"])
    params:
        overhang=config["star_sjdb_overhang"]
    threads: config["threads"]
    shell:
        "mkdir -p {output} && "
        "STAR --runMode genomeGenerate --genomeDir {output} "
        "--genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} "
        "--sjdbOverhang {params.overhang} --runThreadN {threads} "
        "--genomeSAindexNbases 9"  # small bacterial+phage genome, default 14 is oversized

# ---------------------------------------------------------------------------
# Per-sample: download -> trim -> QC -> align -> count
# ---------------------------------------------------------------------------

rule download_fastq:
    output:
        r1="resources/fastq/{sample}_1.fastq.gz",
        r2="resources/fastq/{sample}_2.fastq.gz"
    params:
        srr=srr_for
    shell:
        """
        mkdir -p resources/fastq
        prefetch {params.srr} -O resources/sra_tmp
        fasterq-dump resources/sra_tmp/{params.srr} -O resources/fastq --split-files -e {threads}
        gzip -f resources/fastq/{params.srr}_1.fastq resources/fastq/{params.srr}_2.fastq
        mv resources/fastq/{params.srr}_1.fastq.gz {output.r1}
        mv resources/fastq/{params.srr}_2.fastq.gz {output.r2}
        """

rule trim:
    input:
        r1="resources/fastq/{sample}_1.fastq.gz",
        r2="resources/fastq/{sample}_2.fastq.gz"
    output:
        r1="results/trimmed/{sample}_1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_2.trimmed.fastq.gz"
    log:
        "results/logs/cutadapt_{sample}.log"
    shell:
        "mkdir -p results/trimmed results/logs && "
        "cutadapt -q 20 --minimum-length 25 -a AGATCGGAAGAGC -A AGATCGGAAGAGC "
        "-o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log} 2>&1"

rule fastqc:
    input:
        r1="results/trimmed/{sample}_1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_2.trimmed.fastq.gz"
    output:
        directory("results/fastqc/{sample}")
    shell:
        "mkdir -p {output} && fastqc -o {output} {input.r1} {input.r2}"

rule align:
    input:
        r1="results/trimmed/{sample}_1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_2.trimmed.fastq.gz",
        index=config["star_index_dir"]
    output:
        bam="results/aligned/{sample}/Aligned.sortedByCoord.out.bam"
    params:
        outdir="results/aligned/{sample}/"
    threads: config["threads"]
    shell:
        "mkdir -p {params.outdir} && "
        "STAR --runThreadN {threads} --genomeDir {input.index} "
        "--readFilesIn {input.r1} {input.r2} --readFilesCommand zcat "
        "--outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.outdir} && "
        "samtools index {output.bam}"

rule count:
    input:
        bam="results/aligned/{sample}/Aligned.sortedByCoord.out.bam",
        gtf=config["combined_gtf"]
    output:
        "results/counts/{sample}.featureCounts.txt"
    threads: config["threads"]
    shell:
        "mkdir -p results/counts && "
        "featureCounts -T {threads} -p --countReadPairs -a {input.gtf} "
        "-o {output} {input.bam}"

# ---------------------------------------------------------------------------
# Aggregate + plot
# ---------------------------------------------------------------------------

rule compute_shutoff:
    input:
        counts=expand("results/counts/{sample}.featureCounts.txt", sample=SAMPLE_IDS),
        samples=config["sample_sheet"]
    output:
        "results/shutoff_timecourse.tsv"
    log:
        "results/logs/compute_shutoff.log"
    params:
        prefix=config["phage_prefix"]
    shell:
        "python scripts/compute_shutoff.py --counts {input.counts} "
        "--samples {input.samples} --phage-prefix {params.prefix} "
        "--out {output} --log {log}"

rule plot_shutoff:
    input:
        "results/shutoff_timecourse.tsv"
    output:
        "results/plots/shutoff_stacked_area.png",
        "results/plots/shutoff_by_replicate.png"
    shell:
        "Rscript scripts/plot_shutoff.R {input} results/plots/"
