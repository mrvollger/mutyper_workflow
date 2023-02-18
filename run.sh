#!/usr/bin/env bash
set -euo pipefail
snakemake -k -p --configfile config/config.yaml --use-conda --cores 120 --notemp --rerun-incomplete $@

