import sys
import random
from Bio import SeqIO
from Bio.Data import IUPACData
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main():
    fa_path = sys.argv[1]
    NUC = ["A", "C", "G", "T", "N"]
    IUPAC = IUPACData.ambiguous_dna_values
    for record in SeqIO.parse(fa_path, "fasta"):
        newseq = ""
        for c in record.seq:
            c = c.upper()
            if c == "-" or c in NUC:
                newseq += c
            else:
                newseq += random.choice(IUPAC[c])
        record = SeqRecord(
            Seq(newseq), id=record.id, name=record.name, description=record.description
        )
        SeqIO.write(record, sys.stdout, "fasta")


if __name__ == "__main__":
    main()
