#!env bash
set -o pipefail

Rscript -e "library(devtools); pak::pak('shahcompbio/signals', upgrade=FALSE, dependencies=FALSE)"
