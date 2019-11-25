# Snakemake Language Support

Provides basic language support and snippets for [Snakemake](https://snakemake.readthedocs.io) files (Snakefile, *.smk).

This is an unofficial extension and I am not affiliated with the [Snakemake](https://snakemake.readthedocs.io) project, but since I use Snakemake daily (it is truly an excellent tool) I wanted better support for it in VSCode and in publishing this plugin I hope that others might get some value from it as well.

Feedback, suggestions, and contributions are very welcome!

## Features

- Syntax definitions based on Python, with added Snakemake grammar
- Language rules
- Snippets

## Example

![Snakemake syntax highlighting example](misc/example.png)

<!--
Needs to be published with:
vsce publish --baseContentUrl https://gitlab.com/alping/vscode-snakemake/raw/master
-->

Example taken from [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html#summary).

## Supported Syntax

<details>

 <summary>Rule types</summary>

  - rule
  - subworkflow
  - checkpoint

</details>

<details>

 <summary>Configs</summary>

  - include
  - configfile
  - snakefile
  - workdir
  - ruleorder
  - onsuccess
  - onerror
  - onstart

</details>

<details>

 <summary>Rule Parameters</summary>

  - input
  - output
  - params
  - log
  - priority
  - shadow
  - group
  - benchmark
  - message
  - threads
  - resources
  - version
  - run
  - shell
  - script
  - cwl
  - conda
  - wildcard_constraints
  - wrapper

</details>

<details>

 <summary>Functions</summary>

  - expand
  - unpack
  - ancient
  - directory
  - temp
  - protected
  - touch
  - pipe

</details>

## TODO

- [ ] Indentation rules
- [ ] Recognize string substitutions `"command {input}"`
- [ ] Recognize wildcard constraints inside string substitutions `"{sample,[A-Za-z0-9]+}"`

## Snakemake Support for other Editors

- [Vim](https://github.com/snakemake/snakemake/tree/master/misc/vim)
