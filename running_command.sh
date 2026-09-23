#!/bin/sh
# Use the snakemake of the active conda env, not a stale one earlier on PATH.

if [ -z "$CONDA_PREFIX" ]; then
    echo "No conda environment active. Run: conda activate Ribo_DT" >&2
    exit 1
fi

SNAKEMAKE="$CONDA_PREFIX/bin/snakemake"

if [ ! -x "$SNAKEMAKE" ]; then
    echo "snakemake not found in $CONDA_PREFIX/bin - is Ribo_DT the active env?" >&2
    exit 1
fi

"$SNAKEMAKE" -s Snakefile -j 999 --cluster-config cluster.json \
    --cluster "sbatch --cpus-per-task {cluster.n} --time {cluster.time} --mem {cluster.mem}"
