#!/bin/sh
# Run snakemake from the ACTIVE conda environment.
# A stale "pip install --user" snakemake (5.4.2, spack python 3.6) sits in
# ~/.local/bin and can win on PATH, failing with a pkg_resources/urllib3 error.
# Calling $CONDA_PREFIX/bin/snakemake explicitly sidesteps whatever PATH order
# the shell ends up with.

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
