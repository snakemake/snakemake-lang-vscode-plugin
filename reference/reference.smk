# Modularization
# https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#modularization
include: "path/to/other.snakefile"

# General reference
rule reference_rule:
    input:
        "input/main_input.csv",
        expand("input/sample_input_{sample}.csv"], sample=range(10))
    output:
        "output/main_output.csv",
        protected("output/protected_output.csv"),
        temp("output/temp_output.csv"),
        expand("output/sample_output_{sample}.csv"], sample=range(10))
    conda:
        "envs/environment.yml"
    params:
        threads = "4"
    shell:
        "fastqc -o data/fastqc -t {params.threads} {input}"

# Benchmarking
# https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#benchmarking
rule bwa_map:
    input:
        "data/genome.fa",
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        temp("mapped_reads/{sample}.bam")
    params:
        rg="@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_mem/{sample}.log"
    benchmark:
        "benchmarks/{sample}.bwa.benchmark.txt"
    threads: 8
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"


# Wildcard constraints
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards
rule complex_conversion:
    input:
        "{dataset}/inputfile"
    output:
        "{dataset}/file.{group}.txt"
    wildcard_constraints:
        dataset="\d+"
    shell:
        "somecommand --group {wildcards.group}  < {input}  > {output}"


# Tool wrappers
# https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#tool-wrappers
rule bwa_mem:
    input:
        ref="data/genome.fa",
        sample=lambda wildcards: config["samples"][wildcards.sample]
    output:
        temp("mapped_reads/{sample}.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        "-R '@RG\tID:{sample}\tSM:{sample}'"
    threads: 8
    wrapper:
        "0.15.3/bio/bwa/mem"

# Input Functions and unpack()
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#input-functions-and-unpack
def myfunc(wildcards):
    return { 'foo': '{wildcards.token}.txt'.format(wildcards=wildcards)

rule:
    input: unpack(myfunc)
    output: "someoutput.{token}.txt"
    shell: "..."
