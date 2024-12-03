import sys
import os
import glob
import Levenshtein
import re
import numpy as np
from Bio import SeqIO

from path_fact import get_min_recs


def parse_truth(fpath):
    truth = {}
    for line in open(fpath):
        subpath1, subpath2 = [], []
        if "is a mosaic with" in line:
            pathname = line.split(" ")[2]
            L = eval(" ".join(line.strip("\n").split(" ")[10:]))
            pid1, pid2 = L[0][1], L[-1][1]
            first = True
            for (n1, n2), p in L[1:-1]:
                if p == pid2:
                    first = False
                if first:
                    subpath1.append(str(n1))
                else:
                    subpath2.append(str(n1))
            subpath2.append(str(n2))
            # print("T: ", pid1, pid2, subpath1, subpath2)
            truth[pathname] = (pid1 - 1, pid2 - 1, subpath1, subpath2)
    return truth


def get_seq(record):
    seq, p = "", -1
    # FIXME assuming one record
    seq = record.strip("\n").split("\t")[13]
    # p = int(record.split("\t")[-1])
    return seq


def get_seq_ga(record, nodes):
    p1, p2 = record.split("\t")[7:9]
    p1, p2 = int(p1), int(p2)
    path = get_path(record)
    seq = "".join(nodes[nidx] for nidx in path)
    return seq[p1 : p2 + 1]


def jaccard(S1, S2):
    return round(len(S1 & S2) / len(S1 | S2), 2)


def parse_gfa(gfa_path):
    nodes, paths = {}, {}
    for line in open(gfa_path):
        if line.startswith("S"):
            _, idx, seq = line.strip("\n").split("\t")
            nodes[idx] = seq
        elif line.startswith("P"):
            _, idx, P, _ = line.strip("\n").split("\t")
            path = [p[:-1] for p in P.split(",")]
            paths[idx] = path
    return nodes, paths


def get_path(record):
    return record.split("\t")[5][1:].split(">")


def get_score(record):
    return float(re.findall("score: (-?[0-9]+[\.0-9]*)", record)[0])


def get_recombination(record):
    pid1, pid2 = -1, -1  # paths ID
    subpath1, subpath2 = [], []  # subpaths
    p1, p2 = -1, -1  # positions on paths (prefix/suffix)
    # FIXME assuming one record
    path = record.split("\t")[5]
    path = [int(p) for p in path[1:].split(">")]
    info = record.split("\t")[12]
    if "recombination path" in info:
        # paths
        pid1, pid2 = info.split(" ")[3:5]
        pid1 = int(pid1)
        pid2 = int(pid2[:-1])
        # nodes
        n1, n2 = info.split(" ")[6:8]
        n1, p1 = n1.split("[")
        n1 = int(n1)
        p1 = int(p1[:-1])
        n2, p2 = n2.split("[")
        n2 = int(n2)
        p2 = int(p2[:-2])
        first = True
        for p in path:
            if p == n2:
                first = False
            if first:
                subpath1.append(p)
            else:
                subpath2.append(p)
    return pid1, subpath1, p1, pid2, subpath2, p2


def main():
    in_dir = sys.argv[1]
    print(
        "Aligned",
        "Tool",
        "Gene",
        "Mosaic",
        "PercErr",
        "Score",
        "IsInputPath",
        "AlignedSeqLen",
        "TrueSeqLen",
        "DeltaLen",
        "AlignedSeqLen==TrueSeqLen",
        "ED",
        "RecClass",
        "Rec",
        "TrueRec",
        "Jaccard",
        "MinExpRec",
        "Path",
        "TruePath",
        "Seq",
        "TrueSeq",
        sep="\t",
    )
    for d in glob.glob(os.path.join(in_dir, "*/")):
        gene = d.split("/")[-2]

        gfa_path = os.path.join(d, "MPCSIM", "graph.gfa")
        if not os.path.isfile(gfa_path):
            continue
        gfa_nodes, gfa_paths = parse_gfa(gfa_path)

        truth_path = os.path.join(d, "MPCSIM", "recombinations.log")
        truth = parse_truth(truth_path)

        for gaf_path in glob.glob(os.path.join(d, "MPCSIM", "*", "p*", "*.gaf")):
            mosaic, p, fname = gaf_path.split("/")[-3:]
            p = int(p[1:])
            fname = fname.split(".")[0]
            if fname == "recgraph-9" or fname == "recgraph-5" or fname == "recgraph-8":
                continue

            seq_path = os.path.join(
                os.path.join(d, "MPCSIM", mosaic, f"p{p}", "sequence.fa")
            )
            true_seq = str(next(SeqIO.parse(seq_path, "fasta")).seq)

            true_pid1, true_pid2, true_subpath1, true_subpath2 = truth[mosaic]
            true_path = true_subpath1 + true_subpath2
            true_path_str = "-".join(true_path)

            flag = "X"
            score = -1
            is_input_path = False
            aligned_seq = ""
            ed = len(true_seq)
            hclass = "."
            jacc = 0
            min_rec = -1

            gaf_record = open(gaf_path).readline()
            if gaf_record != "":
                flag = "O"
                aligned_path = get_path(gaf_record)
                aligned_path_str = "-".join(aligned_path)

                # Check if aligned path is in input .gfa (subpath)
                is_input_path = False
                for path in gfa_paths.values():
                    if aligned_path_str in "-".join(path):
                        is_input_path = True
                        break

                # alignment score
                if fname != "graphaligner":
                    score = get_score(gaf_record)

                # recombination
                pid1, subpath1, p1, pid2, subpath2, p2 = -1, [], -1, -1, [], -1
                hclass = "."  # hit class for recombination
                if fname == "recgraph-8":
                    pid1, subpath1, p1, pid2, subpath2, p2 = get_recombination(
                        gaf_record
                    )
                    if pid1 == true_pid1 and pid2 == true_pid2:
                        hclass = "BOTH"
                    elif pid1 == true_pid1:
                        hclass = "FIRST"
                    elif pid2 == true_pid2:
                        hclass = "SECOND"
                    else:
                        if pid1 == -1 and pid2 == -1:
                            hclass = "NONE"
                        else:
                            hclass = "WRONG"

                # ED between path and input seq
                if fname == "graphaligner":
                    aligned_seq = get_seq_ga(gaf_record, gfa_nodes)
                else:
                    aligned_seq = get_seq(gaf_record)
                ed = Levenshtein.distance(aligned_seq, true_seq)

                # Jaccard between path and input path
                jacc = jaccard(set(aligned_path), set(true_path))

                # how many recombinations can explain the aligned path?
                min_rec = get_min_recs(gfa_path, gaf_path)

            print(
                flag,
                fname,
                gene,
                mosaic,
                p,
                score,
                is_input_path,
                len(true_seq),
                len(aligned_seq),
                len(aligned_seq) - len(true_seq),
                aligned_seq == true_seq,
                ed,
                hclass,
                f"{pid1}>{pid2}" if pid1 != -1 else ".",
                f"{true_pid1}>{true_pid2}",
                jacc,
                min_rec,
                aligned_path_str,
                true_path_str,
                aligned_seq,
                true_seq,
                sep="\t",
            )


if __name__ == "__main__":
    main()
