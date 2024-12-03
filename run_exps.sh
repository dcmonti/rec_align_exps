#!/bin/bash

rm -rf output
bash get_HLA_full.sh
snakemake -s generate_HLA_reads.smk -c 4
bash split_reads.sh
snakemake rg-vs-ra.smk -c 4
snakemake align_all_HLA.smk -c 4
snakemake -s sars_cov2_reads_graphs.smk -c 1
snakemake -s sars_cov2_align.smk -c 4