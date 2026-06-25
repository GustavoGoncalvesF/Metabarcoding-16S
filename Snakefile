configfile: "config.yaml"

import os
from glob import glob

RAW_DIR = config["raw_dir"]

# detecta R1 automaticamente
SAMPLES = sorted([
    os.path.basename(f).replace("_R1.fastq.gz", "")
    for f in glob(f"{RAW_DIR}/*_R1.fastq.gz")
])

rule all:
    input:
        # stats
        "results/seqkit_stats_raw.txt",
        "results/seqkit_stats_post_cutadapt.txt",
        "results/seqkit_stats_post_filter.txt",

        # outputs finais
        expand("results/trimmed/{sample}_R1.final.fastq.gz", sample=SAMPLES),
        expand("results/trimmed/{sample}_R2.final.fastq.gz", sample=SAMPLES)


###################################################
# 1. STATS RAW
###################################################
rule seqkit_stats_raw:
    input:
        expand("{raw_dir}/{sample}_R1.fastq.gz", raw_dir=RAW_DIR, sample=SAMPLES)
    output:
        "results/seqkit_stats_raw.txt"
    shell:
        """
        seqkit stats {RAW_DIR}/*.fastq.gz > {output}
        """


###################################################
# 2. CUTADAPT (remove primers)
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
        rev=config["cutadapt"]["rev_primer"]
    shell:
        """
        cutadapt \
        -g {params.fwd} \
        -G {params.rev} \
        --discard-untrimmed \
        -o {output.R1} -p {output.R2} \
        {input.R1} {input.R2}
        """


###################################################
# 3. STATS POST CUTADAPT
###################################################
rule seqkit_stats_cutadapt:
    input:
        expand("results/trimmed/{sample}_R1.cut.fastq.gz", sample=SAMPLES)
    output:
        "results/seqkit_stats_post_cutadapt.txt"
    shell:
        """
        seqkit stats results/trimmed/*_cut.fastq.gz > {output}
        """


###################################################
# 4. LENGTH FILTER (via cutadapt)
###################################################
rule length_filter:
    input:
        R1="results/trimmed/{sample}_R1.cut.fastq.gz",
        R2="results/trimmed/{sample}_R2.cut.fastq.gz"
    output:
        R1="results/trimmed/{sample}_R1.final.fastq.gz",
        R2="results/trimmed/{sample}_R2.final.fastq.gz"
    params:
        minlen=config["cutadapt"]["min_length"]
    shell:
        """
        cutadapt \
        --minimum-length {params.minlen} \
        -o {output.R1} -p {output.R2} \
        {input.R1} {input.R2}
        """

###################################################
# 5. STATS FINAL
###################################################
rule seqkit_stats_final:
    input:
        expand("results/trimmed/{sample}_R1.final.fastq.gz", sample=SAMPLES)
    output:
        "results/seqkit_stats_post_filter.txt"
    shell:
        """
        seqkit stats results/trimmed/*_final.fastq.gz > {output}
        """