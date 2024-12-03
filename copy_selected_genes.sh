#!/bin/bash

indir=$1
csv=$2
odir=$3

mkdir $odir
cut -f1 -d',' $csv | while read gene
do
    python3 scripts/clean_fa_from_iupac.py $indir/$gene > $odir/$gene
    samtools faidx $odir/$gene
done