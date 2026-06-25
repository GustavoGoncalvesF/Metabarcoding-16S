SAMPLES = ["AL10_S342"]

rule all:
    input:
        "results/seqkit_stats.txt"

rule seqkit_stats:
    input:
        expand("raw/{sample}_R1.fastq.gz", sample=SAMPLES)
    output:
        "results/seqkit_stats.txt"
    shell:
        """
        seqkit stats raw/*fastq.gz > {output}
        """