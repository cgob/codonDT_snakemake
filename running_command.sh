#!/bin/sh
snakemake  -s Snakefile  -j 999 --cluster-config cluster.json --cluster "sbatch --cpus-per-task {cluster.n}  --time {cluster.time} --mem {cluster.mem}"
