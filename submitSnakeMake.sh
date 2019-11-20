#!/bin/bash

conda activate RecBlat

snakemake \
    -kp \
    --use-conda \
    --latency-wait 60 \
    --ri \
    -j 999 \
    --cluster-config /project2/gilad/juanvazquez/projects/smRecSearch/cluster.json \
    -c "sbatch \
        --mem={cluster.mem} \
        --nodes={cluster.n} \
        --tasks-per-node={cluster.tasks} \
        --partition=broadwl \
        --job-name={cluster.name} \
	--output={cluster.logfile}" \
    $*
