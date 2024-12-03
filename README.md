
### Data
In our experiments we selected 5 species from [panX](https://pangenome.org) but these scripts can be used with any specie.

```bash
# Get a .msa from pangenome.org
wget https://data.master.pangenome.org/dataset/Burkholderia_pseudomallei/core_gene_alignments.zip
unzip -q core_gene_alignments.zip
cd core_genes
# Remove proteins
rm *_aa_aln.fa
# Removes .msa with fewer strains (not needed)
# n=... ; grep -c "^>" *.fa | grep -v ":${n}" | cut -f1 -d':' | xargs rm
# Removes / and - from headers since they break everything
for fa in $(ls *.fa) ; do sed -i "s/\//-/g;s/|/-/g" $fa ; samtools faidx $fa ; done
cd ..
```

### Software
1. Install [RecGraph](https://github.com/AlgoLab/RecGraph)
2. All other dependencies are available on conda
```
mamba create -c bioconda -n rg-exps snakemake-minimal make_prg pandas seaborn biopython graphaligner vg odgi pggb samtools
```

### Experiment 1
In this experiment, we build graphs using `make_prg`, simulate recombinants using minimum path cover, and then align the recombinants to a reduced graph using `RecGraph` and `GraphAligner`.

```bash
# Select 100 random genes from the core_genes directory previously created (see Data section)
python3 scripts/select_random_genes.py core_genes 100 > core_genes_random100.csv

# Prepare folder with selected genes. Cleans .msa from IUPAC
bash copy_selected_genes.sh core_genes core_genes_random100.csv core_genes_random100

# Build graphs and compute minimal path cover (-> recombinants + reduced graph)
snakemake -s makegraphs.smk -c 32 -p --config seqsdir=core_genes_random100 recgraph=[/PATH/TO/RECGRAPH/BIN]

# Extract mosaics
bash get_mosaics.sh core_genes_random100

# Add noise to mosaics
snakemake -s addnoise.smk -c 16 -p --config seqsdir=core_genes_random100

# Align mosaics back to reduced graph
snakemake -s align.smk -c 16 -p --config seqsdir=core_genes_random100 recgraph=[/PATH/TO/RECGRAPH/BIN]

# results are in core_genes_random100.results.txt
```

### Experiment 2
In this experiment, we build graphs using `make_prg`, align FASTA files in `/data/cdifficile/`
(Adjust cores and memory depending on your setup.)

```bash
snakemake -s clost_diff.smk --use-conda -p --cores 16 --resources mem_mb=100000 -- output/cdifficile/simulated_recgraph_alone.csv
snakemake -s clost_diff.smk --use-conda -p --cores 16 --resources mem_mb=100000 -- output/cdifficile/full.csv

# results are in output/cdifficile/{simulated_recgraph_alone.csv,full.csv}
```
