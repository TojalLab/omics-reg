#!/bin/bash
mkdir -p joblogs
exec /conda/bin/snakemake \
	-s Snakefile.py \
	--cluster "sbatch -J {params.jobname} --parsable --nodes=1 --cpus-per-task={threads} --time=500:00:00 --mem={resources.mem_gb}G -o joblogs/%x.%j.out -e joblogs/%x.%j.err" \
	-k -p -j20 \
	--rerun-incomplete \
	--use-conda \
	--rerun-triggers mtime \
	--cluster-cancel scancel \
	$@
