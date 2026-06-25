rule all:
    input:
        "results/hello.txt"

rule hello:
    output:
        "results/hello.txt"
    shell:
        "echo 'Snakemake funcionando no servidor!' > {output}"