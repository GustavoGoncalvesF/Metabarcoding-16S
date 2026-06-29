configfile: "config.yaml"

import os
from glob import glob

RAW_DIR = config["raw_dir"]

# Reading samples ID
SAMPLES = sorted([
    os.path.basename(f).replace("_R1.fastq.gz", "")
    for f in glob(f"{RAW_DIR}/*_R1.fastq.gz")
])

rule all:
    input:
        # stats
        "results/stats/seqkit_stats_raw.txt",
        "results/stats/seqkit_stats_post_cutadapt.txt",
        "results/stats/seqkit_stats_post_filter.txt",

        # trimmed files
        expand("results/trimmed/{sample}_R1.final.fastq.gz", sample=SAMPLES),
        expand("results/trimmed/{sample}_R2.final.fastq.gz", sample=SAMPLES),

        # FASTQC files
        expand("results/fastqc/raw/{sample}_{read}_fastqc.html", sample=SAMPLES, read=["R1", "R2"]),
        expand("results/fastqc/raw/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=["R1", "R2"]),
        expand("results/fastqc/filtered/{sample}_{read}_fastqc.html", sample=SAMPLES, read=["R1", "R2"]),
        expand("results/fastqc/filtered/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=["R1", "R2"]),

        # MULTIQC files
        "results/multiqc/raw/multiqc_report.html",
        directory("results/multiqc/raw/multiqc_data"),
        "results/multiqc/filtered/multiqc_report.html",
        directory("results/multiqc/filtered/multiqc_data")

###################################################
# 1. STATS RAW
###################################################
rule seqkit_stats_raw:
    input:
        R1 = expand("{raw_dir}/{sample}_R1.fastq.gz", raw_dir=RAW_DIR, sample=SAMPLES),
        R2 = expand("{raw_dir}/{sample}_R2.fastq.gz", raw_dir=RAW_DIR, sample=SAMPLES)
    output:
           "results/stats/seqkit_stats_raw.txt"
    shell:
        """
        seqkit stats {input.R1} {input.R2} > {output}
        """

###################################################
# 2. RUN RAW FASTQC
###################################################
rule fastqc_raw:
    input:
        f"{RAW_DIR}/{{sample}}_{{read}}.fastq.gz"
    output:
        html="results/fastqc/raw/{sample}_{read}_fastqc.html",
        zip="results/fastqc/raw/{sample}_{read}_fastqc.zip"
    shell:
        """
        fastqc {input} -o results/fastqc/raw
        """

###################################################
# 3. RUN RAW MULTIQC
###################################################
rule multiqc_raw:
    input:
        expand("results/fastqc/raw/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=["R1", "R2"])
    output:
        html="results/multiqc/raw/multiqc_report.html"
        data=directory("results/multiqc/raw/multiqc_data")
    shell:
        """
        multiqc results/fastqc/raw -o results/multiqc/raw --force
        """


###################################################
# 4. CUTADAPT (remove primers)
###################################################
rule cutadapt:
    input:
        R1=f"{RAW_DIR}/{{sample}}_R1.fastq.gz",
        R2=f"{RAW_DIR}/{{sample}}_R2.fastq.gz"
    output:
        R1="results/trimmed/{sample}_R1.cut.fastq.gz",
        R2="results/trimmed/{sample}_R2.cut.fastq.gz"
    params:
        fwd=config["cutadapt"]["fwd_primer"],
        rev=config["cutadapt"]["rev_primer"],
        threads=config["cutadapt"]["threads"]
    shell:
        """
        cutadapt \
        -j {threads} \
        -g {params.fwd} \
        -G {params.rev} \
        --discard-untrimmed \
        -o {output.R1} -p {output.R2} \
        {input.R1} {input.R2}
        """


###################################################
# 5. STATS POST CUTADAPT
###################################################
rule seqkit_stats_cutadapt:
    input:
        R1 = expand("results/trimmed/{sample}_R1.cut.fastq.gz", sample=SAMPLES),
        R2 = expand("results/trimmed/{sample}_R2.cut.fastq.gz", sample=SAMPLES)
    output:
        "results/stats/seqkit_stats_post_cutadapt.txt"
    shell:
        """
        seqkit stats {input.R1} {input.R2} > {output}
        """


###################################################
# 6. LENGTH FILTER (via cutadapt)
###################################################
rule length_filter:
    input:
        R1="results/trimmed/{sample}_R1.cut.fastq.gz",
        R2="results/trimmed/{sample}_R2.cut.fastq.gz"
    output:
        R1="results/trimmed/{sample}_R1.final.fastq.gz",
        R2="results/trimmed/{sample}_R2.final.fastq.gz"
    params:
        minlen=config["cutadapt"]["min_length"],
        threads=config["cutadapt"]["threads"]
    shell:
        """
        cutadapt \
        -j {threads} \
        --minimum-length {params.minlen} \
        -o {output.R1} -p {output.R2} \
        {input.R1} {input.R2}
        """

###################################################
# 7. STATS FINAL
###################################################
rule seqkit_stats_final:
    input:
        R1 = expand("results/trimmed/{sample}_R1.final.fastq.gz", sample=SAMPLES),
        R2 = expand("results/trimmed/{sample}_R2.final.fastq.gz", sample=SAMPLES)
    output:
        "results/stats/seqkit_stats_post_filter.txt"
    shell:
        """
        seqkit stats {input.R1} {input.R2} > {output}
        """
###################################################
# 8. RUN FILTERED FASTQC
###################################################
rule fastqc_filtered:
    input:
        "results/trimmed/{sample}_{read}.final.fastq.gz"
    output:
        html="results/fastqc/filtered/{sample}_{read}_fastqc.html",
        zip="results/fastqc/filtered/{sample}_{read}_fastqc.zip"
    shell:
        """
        fastqc {input} -o results/fastqc/filtered
        """

###################################################
# 9. RUN FILTERED MULTIQC
###################################################
rule multiqc_filtered:
    input:
        expand("results/fastqc/filtered/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=["R1", "R2"])
    output:
        html="results/multiqc/filtered/multiqc_report.html"
        data=directory("results/multiqc/filtered/multiqc_data")
    shell:
        """
        multiqc results/fastqc/filtered -o results/multiqc/filtered --force
        """

