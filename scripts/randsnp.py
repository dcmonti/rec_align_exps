from Bio import SeqIO


def randsnp(n):
    nucl = ["A", "C", "G", "T"]
    nucl.remove(n)
    return random.choice(nucl)


def main(fpath, perc):
    if fpath == "-":
        lines = sys.stdin.readlines()
    else:
        # with open(fpath) as fin:
        #     lines = fin.readlines()
        record = next(SeqIO.parse(fpath, "fasta"))
        lines = [record.id, str(record.seq)]
    name = lines[0].strip()
    seq = list(lines[1].strip())
    print(f">{name}")  # .RANDSNP-{perc}%")
    perc = float(perc) / 100
    randixs = np.random.choice(len(seq), int(len(seq) * perc), replace=False)
    for ix in randixs:
        seq[ix] = randsnp(seq[ix])
    print("".join(seq))


if __name__ == "__main__":
    import sys
    import numpy as np
    import random

    main(*sys.argv[1:])
