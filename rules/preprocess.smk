from pathlib import Path
configfile: "config.yaml"


rule reads_statistics:
    input:
        r1=get_fastqs(sample_name)["r1"]
    output:
        "results/{sample}-seqtk_stats.txt"
    log:
        "log/{sample}-seqtk.log"
    priority:1
    conda:
        "../envs/env.yaml"
    shell:
        "seqtk  fqchk -q3 {input.r1} > {output}"

#get_reads_stats(file)
rule genome_size:
    input:
        r1=get_fastqs(sample_name)["r1"]
    output:
        "results/{sample}-mash.log" #The log is the target
    priority:1
    conda:
        "../envs/env.yaml"
    shell:
        "mash sketch -k 32 -m 3 -r {input.r1} -o results/sketch_r1  &> {output}" #-m ==> Minimum copies 
        #of each k-mer required to pass noise filter for reads.

#This comes from glob_wildcards

rule params_table:
    input:
        stats="results/{sample}-seqtk_stats.txt",
        gsize="results/{sample}-mash.log"
    output:
        "results/{sample}-params.tsv"
    log:
        "log/{sample}-params.log"
    priority:1
    script:
        "../scripts/get_params.py"


rule reads_correction:
    input:
        r1=get_fastqs(sample_name)["r1"], 
        r2=get_fastqs(sample_name)["r2"],
        pars="results/{name}-params.tsv"
    output:
        "results/lighter/{name}_1.cor.fq.gz",
        "results/lighter/{name}_2.cor.fq.gz"
    log:
        "log/{name}-lighter.log"
    threads:
        config["ncores"]
    conda:
        "../envs/env.yaml"
    run:
        records = read_params("results/{}-params.tsv".format(sample_name))
        shell("lighter -od results/lighter -r {input.r1} -r {input.r2} -K 32 "
        "{records[gsize]} -t {threads} -maxcor 1 &> {log}")


rule merge_pairs:
    input:
        r1="results/lighter/{sample}_1.cor.fq.gz",
        r2="results/lighter/{sample}_2.cor.fq.gz",
        pars="results/{sample}-params.tsv"
    output:
        r1="results/merge/{sample}.notCombined_1.fastq.gz",
        r2="results/merge/{sample}.notCombined_2.fastq.gz",
        merged="results/merge/{sample}.extendedFrags.fastq.gz"
    log:
        "log/{sample}-flash.log"
    params:
        min_overlap=config["min_overlap"]
    threads:
        config["ncores"]
    conda:
        "../envs/env.yaml"
    run:
        stats=read_params("results/{}-params.tsv".format(sample_name))#("results/"+sample_name+"-params.tsv")
        avg_len = stats["avg_len"]
        shell("flash -m {params.min_overlap} -M {avg_len} -d results/merge/ -o {sample_name} -z " 
        "-t {threads} {input.r1} {input.r2} &> {log}")






