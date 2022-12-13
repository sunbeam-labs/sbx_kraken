<img src="https://github.com/sunbeam-labs/sunbeam/blob/stable/docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# sbx_kraken

<!-- badges: start -->
[![Tests](https://github.com/sunbeam-labs/sbx_kraken/actions/workflows/tests.yml/badge.svg)](https://github.com/sunbeam-labs/sbx_kraken/actions/workflows/tests.yml)
[![Super-Linter](https://github.com/sunbeam-labs/sbx_kraken/actions/workflows/linter.yml/badge.svg)](https://github.com/sunbeam-labs/sbx_kraken/actions/workflows/linter.yml)
<!-- badges: end -->

A [Sunbeam](https://github.com/sunbeam-labs/sunbeam) extension for taxonomic assignment of reads to databases using Kraken. You can get pre-built kraken databases at the [Kraken homepage](http://ccb.jhu.edu/software/kraken/) or build your own then specify the path to your database in the config.

## Installation

To install, activate your conda environment (using the name of your environment) and use `sunbeam extend`:

    conda activate <i>sunbeamX.X.X</i>
    sunbeam extend https://github.com/sunbeam-labs/sbx_kraken.git

## Usage

To generate a tsv-format report, create a project, specify your database, and use the `all_classify` target:

    sunbeam init --data_fp /path/to/reads/ /path/to/project/
    sunbeam config modify -i -f /path/to/project/sunbeam_config.yml -s 'sbx_kraken: {{kraken_db_fp: {/path/to/db}}}'
    sunbeam run --profile /path/to/project/ all_classify

N.B. For sunbeam versions <4 the last command will be something like `sunbeam run --configfile /path/to/project/sunbeam_config.yml all_classify`.

## Configuration

  - kraken_db_fp: Is the filepath to a kraken database

## Legacy Installation

For sunbeam versions <3 or if `sunbeam extend` isn't working, you can use `git` directly to install an extension:

    git clone https://github.com/sunbeam-labs/sbx_kraken.git extensions/sbx_kraken

and then include it in the config for any given project with:

    cat extensions/sbx_kraken/config.yml >> /path/to/project/sunbeam_config.yml
