# Change Log

## [0.5.1](https://github.com/snakemake/snakemake-lang-vscode-plugin/compare/v0.5.0...v0.5.1) (2025-03-11)


### Bug Fixes

* added missing keywords for highlighting ([#30](https://github.com/snakemake/snakemake-lang-vscode-plugin/issues/30)) ([e623a48](https://github.com/snakemake/snakemake-lang-vscode-plugin/commit/e623a483fc5af15aea9c9d6cf199262dfa84ff7a))

## [0.5.0](https://github.com/snakemake/snakemake-lang-vscode-plugin/compare/v0.4.0...v0.5.0) (2024-10-20)


### Features

* Add `retries` keyword ([#31](https://github.com/snakemake/snakemake-lang-vscode-plugin/issues/31)) ([ea88750](https://github.com/snakemake/snakemake-lang-vscode-plugin/commit/ea887508abb81f5ed151da95e19358f8a2957e45))

## [0.4.0](https://github.com/snakemake/snakemake-lang-vscode-plugin/compare/v0.3.0...v0.4.0) (2023-08-07)


### Features

* Add new `localrule` ruleparam to keywords.yaml ([#24](https://github.com/snakemake/snakemake-lang-vscode-plugin/issues/24)) ([b046590](https://github.com/snakemake/snakemake-lang-vscode-plugin/commit/b046590fb61d50a176649a336e6037451e3ba055))
* auto closing rule for triple double quotes ([0043178](https://github.com/snakemake/snakemake-lang-vscode-plugin/commit/004317818211ce74501de26bb9eb4a1d400795de))


### Bug Fixes

* restrictions to function keywords (word boundary) ([#16](https://github.com/snakemake/snakemake-lang-vscode-plugin/issues/16)) ([0043178](https://github.com/snakemake/snakemake-lang-vscode-plugin/commit/004317818211ce74501de26bb9eb4a1d400795de))

## [0.3.0](https://github.com/snakemake/snakemake-lang-vscode-plugin/compare/v0.2.0...v0.3.0) (2023-08-07)


### Features

* add template_engine keyword ([16cd7d4](https://github.com/snakemake/snakemake-lang-vscode-plugin/commit/16cd7d4e02c4ac8685fec46c4d61f84ec90827d3))

## [0.2.0](https://github.com/snakemake/snakemake-lang-vscode-plugin/compare/v0.1.8...v0.2.0) (2022-10-09)


### Features

* update keywords and implement autoindentation ([#17](https://github.com/snakemake/snakemake-lang-vscode-plugin/issues/17)) ([7762354](https://github.com/snakemake/snakemake-lang-vscode-plugin/commit/776235401264754bd14a16a944f70f0e45da1cda))

## [1.1.0](https://github.com/ftabaro/snakemake-lang-vscode-plugin/compare/v1.0.0...v1.1.0) (2022-08-26)


### Features

* Implement onEnterRules ([e77976f](https://github.com/ftabaro/snakemake-lang-vscode-plugin/commit/e77976ff72656fb1e72f026d93bc6368311d5662))

## 1.0.0 (2022-08-25)


### Bug Fixes

* add more rule keywords ([e99bc9c](https://github.com/ftabaro/snakemake-lang-vscode-plugin/commit/e99bc9cff84d3c4f9045487a57e62aad26767a96))
* adjusted release type ([8f15738](https://github.com/ftabaro/snakemake-lang-vscode-plugin/commit/8f15738c3e54c3655a0a5a2c0099380b097ee76f))

### [0.1.8](https://github.com/snakemake/snakemake-lang-vscode-plugin/compare/v0.1.7...v0.1.8) (2022-03-30)


### Bug Fixes

* adjusted release type ([8f15738](https://github.com/snakemake/snakemake-lang-vscode-plugin/commit/8f15738c3e54c3655a0a5a2c0099380b097ee76f))

### [0.1.7](https://github.com/snakemake/snakemake-lang-vscode-plugin/compare/v0.1.6...v0.1.7) (2022-03-30)


### Bug Fixes

* add more rule keywords ([e99bc9c](https://github.com/snakemake/snakemake-lang-vscode-plugin/commit/e99bc9cff84d3c4f9045487a57e62aad26767a96))

## 0.1.1 (2019-11-26)

- Updated Changelog

## 0.1.0 (2019-11-26)

- Updated build process to use templating for the keyword regular expressions
- Added script snippet
- Bumped version number to 0.1.0 as the extension seems reasonably stable now

## 0.0.6 (2019-10-21)

- Updated --baseContentUrl to get example image to work

## 0.0.5 (2019-10-18)

- New keyword: singularity
- Added example image to README

## 0.0.4 (2019-10-18)

- Fixed faulty Rule snippet
- Converted definitions to YAML source
- Removed indentation rules until they work properly

## 0.0.3 (2019-10-17)

- New keywords: subworkflow, checkpoint, configfile, snakefile, ruleorder, localrules, onsuccess, onerror, onstart, priority, shadow, group, cwl
- New functions: ancient, directory, touch, pipe

## 0.0.2 (2019-10-17)

- New keywords: conda, wildcard_constraints, wrapper
- New functions: expand, unpack, temp, protected
- Support for unnamed rules

## 0.0.1 (2019-10-15)

- Initial release
