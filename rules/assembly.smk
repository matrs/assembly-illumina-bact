
rule assembly:
    input:
        # estimating_kmers("results/{sample}-seqtk_stats.txt"),
        r1="results/merge/{sample}.notCombined_1.fastq.gz",
        r2="results/merge/{sample}.notCombined_2.fastq.gz",
        merged="results/merge/{sample}.extendedFrags.fastq.gz"
    output:
        "results/assembly/{sample}-scaffolds.fasta",
        "results/assembly/{sample}-contigs.fasta"
    log:
        "log/{sample}-spades.log"
    threads:
        config["ncores"]
    conda:
        "../envs/env.yaml"
    run:
        kmers = estimating_kmers("results/"+sample_name+"-params.tsv")
        shell("spades.py --pe1-1 {input.r1} --pe1-2 {input.r2} --only-assembler --threads {threads} "
        "--memory 8 -o results/assembly/ -k {kmers} --s2 {input.merged} &> {log}; "
        "mv results/assembly/scaffolds.fasta results/assembly/{sample_name}-scaffolds.fasta; "
        "mv results/assembly/contigs.fasta results/assembly/{sample_name}-contigs.fasta")