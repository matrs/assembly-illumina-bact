# Snakemake workflow: assembly-illumina-bacteria

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.6-brightgreen.svg)](https://snakemake.bitbucket.io)


This pipeline is designed to assemble short genomic reads from bacteria produced by illumina sequencers. It's based on the [shovill](https://github.com/tseemann/shovill) pipeline and uses `spades` to generate an assembly and `pilon` to correct it. It automatically computes k-mer sizes to pass to `spades` and then align the original reads to the assembly using `bwa mem`, passing the resulting `bam` file to `pilon`. Furthermore, performs quality control steps using `quast`, `busco` and `multiqc` to summarize general statistics. See the picture of the DAG at the end of this document for more details.
## Authors

* Jose Maturana (@matrs)

## Usage

### Simple

#### Step 1: Install workflow

If you simply want to use this workflow, download and extract this repository.
If you intend to modify and further extend this workflow or want to work under version control, fork this repository as outlined in [Advanced](#advanced). The latter way is recommended.

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository and, if available, its DOI (see above).

#### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`.

#### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

If you not only want to fix the software stack but also the underlying OS, use

    snakemake --use-conda --use-singularity

in combination with any of the modes above.
See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

Use 

### Advanced

The following recipe provides established best practices for running and extending this workflow in a reproducible way.

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to the desired working directory for the concrete project/run on your machine.
3. [Create a new branch](https://git-scm.com/docs/gittutorial#_managing_branches) (the project-branch) within the clone and switch to it. The branch will contain any project-specific modifications (e.g. to configuration, but also to code).
4. Modify the config, and any necessary sheets (and probably the workflow) as needed.
5. Commit any changes and push the project-branch to your fork on github.
6. Run the analysis.
7. Optional: Merge back any valuable and generalizable changes to the [upstream repo](https://github.com/snakemake-workflows/assembly-illumina-bact) via a [**pull request**](https://help.github.com/en/articles/creating-a-pull-request). This would be **greatly appreciated**.
8. Optional: Push results (plots/tables) to the remote branch on your fork.
9. Optional: Create a self-contained workflow archive for publication along with the paper (snakemake --archive).
10. Optional: Delete the local clone/workdir to free space.

### Pipeline's directed acyclic graph

![DAG](/dag_rules.svg)
