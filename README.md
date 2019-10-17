# Snakemake Language Support

Provides basic language support and snippets for [Snakemake](https://snakemake.readthedocs.io) files (Snakefile, *.smk).

This is an unofficial extension and I am not affiliated with the [Snakemake](https://snakemake.readthedocs.io) project.
But since I use Snakemake daily (it is truly an excellent tool), I wanted better support for it in VSCode and I in publishing this plugin I hope others might get some value from it as well.

This is still very much in alpha, so it's likely that not everything will work.
Feedback, suggestions, and contributions are very welcome!

## Features

- Syntax definitions based on Python with added Snakemake grammar
- Language rules based on Python with added Snakemake grammar
- Rule snippet

## TODO

- [ ] Identify Snakemake string substitutions (e.g. `"command {input} > {output}"`)
- [ ] Identify wildcard constraints inside Snakemake string substitutions (e.g. `"sorted_reads/{sample,[A-Za-z0-9]+}.bam"`)
