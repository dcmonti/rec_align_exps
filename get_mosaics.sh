#!bin/bash

indir=$1
WD=$PWD/$(dirname $0) # FIXME

for log in $(ls ${indir}/*/MPCSIM/recombinations.log)
do
    DDIR=$(dirname $log)
    mosaics=$(grep "is a mosaic" ${log} | cut -f3 -d" ")
    for mosaic in $mosaics
    do
        if [[ ! -f $DDIR/$mosaic/sequence.0.fa ]]
        then
            mkdir -p $DDIR/$mosaic
            samtools faidx $DDIR/../seqs.noindels.fa $mosaic > $DDIR/$mosaic/sequence.fa
            cp $DDIR/$mosaic/sequence.fa $DDIR/$mosaic/sequence.0.fa
        fi
    done
done