import sys
import os
import glob
import pandas as pd
from Bio import SeqIO


def main():
    indir = sys.argv[1]
    n = int(sys.argv[2])

    data = []
    for fa in glob.glob(os.path.join(indir, "*.fa")):
        rec = next(SeqIO.parse(fa, "fasta"))
        fname = fa.split("/")[-1]
        l = len(rec)
        data.append([fname, l])
    df = pd.DataFrame(data, columns=["Gene", "l"])

    q1 = df["l"].quantile(0.25)
    q2 = df["l"].quantile(0.75)

    df = df[(df["l"] >= q1) & (df["l"] <= q2)]
    fdf = df.sample(n=n)
    fdf.to_csv(sys.stdout, index=False, header=False)


if __name__ == "__main__":
    main()
