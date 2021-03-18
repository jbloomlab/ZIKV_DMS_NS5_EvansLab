#!/bin/bash

set -e

. "/fh/fast/bloom_j/software/miniconda3/etc/profile.d/conda.sh"
conda activate ZIKV_DMS

RESULTSDIR="results/summary"

mkdir -p $RESULTSDIR

declare -a nbs=(
                "dms_analysis.ipynb"
                )

for nb in "${nbs[@]}"
do
    echo "Running $nb"

    jupyter nbconvert \
        --to notebook \
        --execute \
        --inplace \
        --ExecutePreprocessor.timeout=-1 \
        $nb

    echo "Converting $nb to Markdown"
    jupyter nbconvert \
        --output-dir $RESULTSDIR \
        --to markdown \
        $nb
done
