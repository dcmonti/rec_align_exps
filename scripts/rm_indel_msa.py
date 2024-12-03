import sys

from Bio import SeqIO


def main():
    for record in SeqIO.parse(sys.argv[1], "fasta"):
        print(f">{record.id}")
        print(record.seq.replace("-", ""))


if __name__ == "__main__":
    main()
