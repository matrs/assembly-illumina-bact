rule sam_stats:
    input:
        "results/assembly/{sample}-spades.bam"
    output:
        "qc/{sample}-sam_stats.txt"
    # log:
    #     "log/{sample}-sam_stats.log"
    conda:
        "../envs/env.yaml"
    shell:
            "samtools stats {input}  > {output}"

rule sam_flagstat:
    input:
        "results/assembly/{sample}-spades.bam"
    output:
        "qc/{sample}-sam_flagstat.txt"
    # log:
    #     "log/{sample}-sam_flagstat.log"
    conda:
        "../envs/env.yaml"
    shell:
        "samtools flagstat {input} > {output}"

rule busco:
    input:
        "results/assembly/corrected/{sample}.fasta"
    output:
        "qc/busco/full_table_{sample}.tsv",
        "qc/busco/missing_busco_list_{sample}.tsv",
        "qc/busco/short_summary_{sample}.txt"
    log:
        "log/{sample}-busco.log"
    threads:
        config["ncores"]
    params:
        db=config["busco_db"],
        busco= config["params"]["busco"]
    conda:
        "../envs/env.yaml"
    shell:
        ("run_busco -i {input} -o {sample_name} --lineage {params.db} --mode genome "
        "--cpu {threads} {params.busco} --force &> {log} ; cp -r run_{sample_name}/* "
        "qc/busco; rm -rf run_{sample_name}")


rule quast:
    input:
        nocorr="results/assembly/{sample}-contigs.fasta",
        corr="results/assembly/corrected/{sample}.fasta"
    output:
        "qc/quast-{sample}/report.html",
        "qc/quast-{sample}/report.tsv",
        "qc/quast-{sample}/report.pdf"
    log:
        "log/{sample}-quast.log"
    threads:
        config["ncores"]
    params:
        ref=config["ref_genome"]
    conda:
        "../envs/env.yaml"
    shell:
        ("quast -t  -{threads} -R  {params.ref} {input.nocorr} {input.corr} -o qc/quast-{sample_name} &> {log}")


rule multiqc:
    input:
        "qc/busco/short_summary_{sample}.txt",
        "qc/quast-{sample}/report.tsv",
        "qc/{sample}-sam_stats.txt",
        "qc/{sample}-sam_flagstat.txt"
    output:
        "qc/{sample}-report.html"
    log:
        "log/{sample}-multiqc.log"
    conda:
        "../envs/env.yaml"
    wrapper:
        "0.35.1/bio/multiqc"
