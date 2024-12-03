import sys
import os
import glob


def parse_gfa(gfa_path):
    seqs = {}
    size = 0
    longest_path = 0
    nodes, edges, paths = 0, 0, 0
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


def parse_time(time_path):
    time, ram = 0, 0
    for line in open(time_path):
        if "Elapsed (wall clock)" in line:
            t = line.strip().split(" ")[-1]
            m, sec = t.split(":")
            time = int(m) * 60 + float(sec)
        elif "Maximum resident set" in line:
            r = line.strip().split(" ")[-1]
            ram = float(r)*1e-6 # TODO: CHECKME
    return round(time/60, 2), round(ram, 2)


def main():
    in_dir = sys.argv[1]

    print("Tool", "Time (min)", "RAM (GB)", sep="\t")
    # print("Nodes", "Edges", "Paths", "Size", "LongestPath", "Time", sep="\t")
    for d in glob.glob(os.path.join(in_dir, "*/")):
        gene = d.split("/")[-2]

        # gfa_path = os.path.join(d, "MPCSIM", "graph.gfa")
        # if not os.path.isfile(gfa_path):
        #     continue
        # nodes, edges, paths, size, longest_path = parse_gfa(gfa_path)

        for time_path in glob.glob(
            os.path.join(d, "MPCSIM", "*", "p*", "*.time")
        ):
            tname = time_path.split("/")[-1].split(".")[0]
            if tname in ["graphaligner", "recgraph-8", "recgraph-0", "recgraph-4"]:
                time, ram = parse_time(time_path)
                # print(nodes, edges, paths, size, longest_path, time, ram, sep="\t")
                print(tname, time, ram, sep="\t")


if __name__ == "__main__":
    main()
