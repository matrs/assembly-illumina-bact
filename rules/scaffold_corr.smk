rule bwa_index:
    input:
        # estimating_kmers("results/{sample}-seqtk_stats.txt"),
        "results/assembly/{sample}-contigs.fasta"
    output:
        "results/assembly/{sample}-contigs.fasta.sa",
        "results/assembly/{sample}-contigs.fasta.bwt",
        "results/assembly/{sample}-contigs.fasta.pac",
        "results/assembly/{sample}-contigs.fasta.ann",
        "results/assembly/{sample}-contigs.fasta.amb"
    log:
        "log/{sample}-bwa_index.log"
    priority:1
    conda:
        "../envs/env.yaml"
    shell:
        "bwa index {input} &> {log}"


rule faidx:
    input:
        "results/assembly/{sample}-contigs.fasta"
    output:
        "results/assembly/{sample}-contigs.fasta.fai"
    log:
        "log/{sample}-faidx.log"
    priority:1
    conda:
        "../envs/env.yaml"
    shell:
        "samtools faidx {input} &> {log}"


rule produce_bam:
    input:
        r1=get_fastqs(sample_name)["r1"],
        r2=get_fastqs(sample_name)["r2"],
        contigs="results/assembly/{sample}-contigs.fasta",
        fai="results/assembly/{sample}-contigs.fasta.fai",
        bai="results/assembly/{sample}-contigs.fasta.bwt"#not explicitly passed but needed
    output:
        "results/assembly/{sample}-spades.bam"
    log:
        "log/{sample}-produce_bam.log"
    params:
        mem = int((config["mem"]/2) * 1024 / max( 1, int(0.25 * config["ncores"]))),
        cpus=max( 1, int(0.25 * config["ncores"]))
    conda:
        "../envs/env.yaml"
    shell:
        ("bwa mem -v 3 -x intractg -t {threads} {input.contigs} {input.r1} {input.r2} | " 
        "samclip --ref {input.fai} | samtools sort --threads {params.cpus} " 
        "-m {params.mem}m --reference {input.contigs} -o {output} &> {log}")

rule sam_index:
    input:
        "results/assembly/{sample}-spades.bam"
    output:
        "results/assembly/{sample}-spades.bam.bai"
    log:
        "log/{sample}-sam_index.log"
    priority:1
    conda:
        "../envs/env.yaml"
    shell:
        "samtools index {input} &> {log}"

# Under shell directive, {wildcards} aren't accessible, so {sample_name}, global variable from preprocess.smk, is used
rule assembly_correction:
    input:
        bam="results/assembly/{sample}-spades.bam",
        fasta="results/assembly/{sample}-contigs.fasta",
        bai="results/assembly/{sample}-spades.bam.bai" #not explicitly passed but needed
    output:
        "results/assembly/corrected/{sample}.fasta",
        "results/assembly/corrected/{sample}.changes"
    log:
        "log/{sample}-assembly_corr.log"
    params:
        min_bq=config["min_bq"],
        min_mq=config["min_mq"]
    threads:
        config['ncores']
    conda:
        "../envs/env.yaml"
    shell:
        ("pilon --genome {input.fasta} --frags {input.bam} --minmq {params.min_mq} "
        "--minqual {params.min_bq} --fix bases --output {sample_name} --threads {threads} "
        "--changes --mindepth 0.25 --outdir results/assembly/corrected &> {log}")






