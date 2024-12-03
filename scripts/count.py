import sys
import os
import glob
from Bio import SeqIO
import statistics


def parse_gfa(gfa_path):
    seqs = {}
    size = 0
    longest_path = 0
    nodes, edges, paths = 0, 0, 0
    if os.path.isfile(gfa_path):
        for line in open(gfa_path):
            if line.startswith("S"):
                _, idx, seq = line.strip("\n").split("\t")
                seqs[idx] = seq
                nodes += 1
                size += len(seq)
            elif line.startswith("L"):
                edges += 1
            elif line.startswith("P"):
                _, idx, P, _ = line.strip("\n").split("\t")
                path_len = sum([len(seqs[p[:-1]]) for p in P.split(",")])
                if path_len > longest_path:
                    longest_path = path_len
                paths += 1
    return nodes, edges, paths, size, longest_path


def print_row(L, label):
    print(label, min(L), max(L), round(statistics.mean(L), 1), round(statistics.stdev(L), 1), sep="\t") #, statistics.median(L))


def main():
    in_dir = sys.argv[1]
    L = []
    NREC = []
    INSIZES = []
    INPATHS = []
    REDSIZES = []
    REDPATHS = []

    for d in glob.glob(os.path.join(in_dir, "*/")):
        fa_path = d[:-1] + ".fa"
        log_path = os.path.join(d, "MPCSIM", "recombinations.log")
        in_gfa_path = os.path.join(d, "graph.tsorted.wpaths.gfa")
        red_gfa_path = os.path.join(d, "MPCSIM", "graph.gfa")

        l = len(next(SeqIO.parse(fa_path, "fasta")))
        L.append(l)

        nrec = 0
        if os.path.isfile(log_path):
            for line in open(log_path):
                nrec += "is a mosaic with" in line
        NREC.append(nrec)

        innodes, inedges, inpaths, insize, inlongestpath = parse_gfa(in_gfa_path)
        INSIZES.append(insize)
        INPATHS.append(inpaths)
        
        rednodes, rededges, redpaths, redsize, redlongestpath = parse_gfa(
            red_gfa_path
        )
        REDSIZES.append(redsize)
        REDPATHS.append(redpaths)

    print("")
    print("Metric", "Min", "Max", "Mean", "std", sep="\t")
    print_row(L, "L")
    print_row(NREC, "NREC")
    print_row(INSIZES, "INSIZES")
    print_row(INPATHS, "INPATHS")
    print_row(REDSIZES, "RESIZES")
    print_row(REDPATHS, "REPATHS")
    print("")


if __name__ == "__main__":
    main()
