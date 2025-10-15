#!/bin/bash

# Calculate number of cores to use (1 less than available; 1 if only 1 available)
CORES=$(( $(nproc) - 1 ))
if [ "$CORES" -lt 1 ]; then
  CORES=1
fi

# Run docker container
# Mounts local 'data', 'data-raw' and 'R' folders from repo
# Sets working directory as 'project'
# Runs snakemake pipeline w. appropriate number of cores
docker run --rm \
  -v "$PWD/data-raw":/project/data-raw \
  -v "$PWD/R":/project/R \
  -v "$PWD/data":/project/data \
  -v "$PWD/.snakemake":/project/.snakemake \
  -v "$PWD/snakefile":/project/snakefile \
  -w /project \
  sorenjessen/b2a_res_single_fiber \
  /opt/conda/bin/conda run -n snakemake_env --no-capture-output \
  snakemake --cores "$CORES"