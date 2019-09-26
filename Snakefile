include: "rules/common.smk"
include: "rules/preprocess.smk"
include: "rules/assembly.smk"
include: "rules/scaffold_corr.smk"
include: "rules/qc.smk"

rule all:
    input:
        # "results/assembly/SRR292770-scaffolds.fasta",
        # "results/assembly/SRR292770-contigs.fasta",
        # "results/SRR292770-params.tsv",
        # "results/assembly/SRR292770-spades.bam",
        # "results/assembly/SRR292770-spades.bam.bai",
        # "results/assembly/SRR292770-contigs.fasta.bwt",
        # "results/assembly/corrected/SRR292770.fasta",
        # "results/assembly/corrected/SRR292770.changes",
        "qc/SRR292770-report.html"