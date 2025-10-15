# PowerShell version of the docker run command
docker run --rm `
  -v "${PWD}\data-raw:/project/data-raw" `
  -v "${PWD}\R:/project/R" `
  -v "${PWD}\data:/project/data" `
  -v "${PWD}\.snakemake:/project/.snakemake" `
  -v "${PWD}\snakefile:/project/snakefile" `
  -w /project `
  sorenjessen/b2a_res_single_fiber `
  /opt/conda/bin/conda run -n snakemake_env --no-capture-output snakemake --cores ([Environment]::ProcessorCount - 1)
