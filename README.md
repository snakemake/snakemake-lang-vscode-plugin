# Snakemake Language Support

Provides basic language support for [Snakemake](https://snakemake.readthedocs.io) files (Snakefile, *.smk).
Feedback, suggestions, and contributions are very welcome!

This project has been started by Peter Alping, and can be considered a fork of [this repository](https://gitlab.com/alping/vscode-snakemake).

## Features

- Syntax definitions based on Python, with added Snakemake keywords
- Language rules based on Python
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

<summary>[Keywords and Functions](src/keywords.yaml)</summary>

- Configurations
  - these keywords specified in the format `<keyword>: <values>` within a snakefile.
  - for example:
    ```snakemake
    configfile: "config.yaml"

    onsuccess:
        print("All jobs finished successfully.")

    include: "workflow/report.smk"
    ```

- Rules and Modules
  - `rule <rulename>:` defines a rule with its parameters, followed by an indented body
  - `checkpoint <rulename>:` is a special type of `rule` that allows dynamic determination of workflow steps based on its output.
  - `module <modulename>:` declares a module repository from which rules can be imported
  - `subworkflow` is now deprecated

- Rule Parameters
  - Used in `rule` and `checkpoint` blocks to define inputs, outputs, and other settings.
  - You can re- `used rule <rulename> as <newrulename> ... with` different parameters, but you cannot change the way it runs, including:
    - run
    - shell
    - script
    - notebook
    - wrapper
    - template_engine
    - cwl

- Module Parameters:
  - Define where and how to import rules when using the `module` directive.
  - Example:
    ```snakemake
    module remote:
        snakefile: "local.smk"
        config: {**config, "remote": True}
        replace_prefix: "local/"
        prefix: "remote/"

    use rule * from remote as remote_*
    ```

- Global Variables:
  - The following classes are available without explicit import from the snakemake module:
    - `Path`
    - `WorkflowError`
  - Important variables used across workflows include:
    - snakemake
    - `rules`: reference rules via `rules.<rulename>`
    - workflow
    - `checkpoints`: access outputs of checkpoint rules
    - storage
    - access
    - scatter
    - gather

- Job Parameters:
  - These parameters are only available during job execution and should be used within `run` or `shell` blocks
    - input
    - output
    - params
    - wildcards
    - threads
    - resources
    - log
    - config

- Functions:
  - Built-in helper functions are available for use without needing to import them

</details>

## TODO

- [ ] recognize files provided as string after `include:`, `conda`, `snakefile`, and allow quickly jumping to them
- [ ] add bash-color to docstring in `shell:` block
- [ ] Indentation rules (really tricky for some reason)
- [ ] Recognize wildcard constraints inside string substitutions: `"{sample,[A-Za-z0-9]+}"`

## Snakemake Support for other Editors

- [Vim](https://github.com/snakemake/snakemake/tree/master/misc/vim)
