#!env bash
set -o pipefail

Rscript -e "library(devtools); pak::pak('shahcompbio/signals@v0.14.1', upgrade=FALSE, dependencies=FALSE)"
