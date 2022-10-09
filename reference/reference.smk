# Modularization
# https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#modularization
include: "path/to/other.snakefile"

# Onstart, onsuccess and onerror handlers
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#onstart-onsuccess-and-onerror-handlers
onstart:
    print("Workflow started")

onsuccess:
    print("Workflow finished, no errors")

onerror:
    print("An error occurred")

# Containerization
container: "docker://something/here"
containerized: "docker://something/here"

# Pep schema
pepfile: "path/to/pepconfig.yaml"
pepschema: "https://pepschema.org"

# Handling Ambiguous Rules
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#handling-ambiguous-rules
ruleorder: rule1 > rule2 > rule3

# Local Rules
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#local-rules
localrules: all, foo

# Standard Configuration
# https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#standard-configuration
configfile: "path/to/config.json

# Sub-Workflows
# https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#sub-workflows
subworkflow otherworkflow:
    workdir:
        "../path/to/otherworkflow"
    snakefile:
        "../path/to/otherworkflow/Snakefile"
    configfile:
        "path/to/custom_configfile.yaml"

rule a:
    input:
        otherworkflow("test.txt")
    output: ...
    shell:  ...

# Data-dependent conditional execution
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution
checkpoint somestep:
    input:
        "samples/{sample}.txt"
    output:
        "somestep/{sample}.txt"
    shell:
        "somecommand {input} > {output}"

# General reference
rule reference_rule:
    input:
        "input/main_input.csv",
        ancient("path/to/inputfile"),
        expand("input/sample_input_{sample}.csv"], sample=range(10))
    output:
        "output/main_output.csv",
        directory("path/to/outputdir"),
        protected("output/protected_output.csv"),
        temp("output/temp_output.csv"),
        pipe("test.{i}.txt"),
        expand("output/sample_output_{sample}.csv"], sample=range(10)),
        report("fig1.png", category="main category")
    conda:
        "envs/environment.yml"
    params:
        threads = "4"
    priority: 50
    shadow: "shallow"
    group: "mygroup"
    singularity:
        "docker://something/here"
    container:
        "docker://something/here"
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
    return { 'foo': '{wildcards.token}.txt'.format(wildcards=wildcards) }

rule:
    input: unpack(myfunc)
    output: "someoutput.{token}.txt"
    shell: "..."

# Flag Files
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#flag-files
rule mytask:
    output: touch("mytask.done")
    shell: "mycommand ..."

# Common-Workflow-Language (CWL) support
# https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#common-workflow-language-cwl-support
rule samtools_sort:
    input:
        input="mapped/{sample}.bam"
    output:
        output_name="mapped/{sample}.sorted.bam"
    params:
        threads=lambda wildcards, threads: threads,
        memory="4G"
    threads: 8
    cwl:
        "https://github.com/common-workflow-language/workflows/blob/"
        "fb406c95/tools/samtools-sort.cwl"
