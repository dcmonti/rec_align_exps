import sys
from collections import defaultdict


def print_real_successors(idx, edges, to_remove, new_idx):
    for s in edges[idx]:
        if s not in to_remove:
            print("L", new_idx, "+", s, "+", "0M", sep="\t")
        else:
            print_real_successors(s, edges, to_remove, new_idx)


def print_new_sources(
    idx, vertices, edges, not_sources, to_remove, next_idx, already_printed=set()
):
    for s in edges[idx]:
        if s not in to_remove:
            if s in not_sources and s not in already_printed:
                v = vertices[s].split("\t")
                v[1] = str(next_idx)
                print("\t".join(v), end="")
                print_real_successors(s, edges, to_remove, next_idx)
                next_idx += 1
                already_printed.add(s)
        else:
            next_idx = print_new_sources(
                s, vertices, edges, not_sources, to_remove, next_idx, already_printed
            )
    return next_idx


def main():
    gfa_path = sys.argv[1]

    vertices = {}
    to_remove = set()
    predecessors = {}
    successors = {}
    for line in open(gfa_path):
        if not line.startswith("S"):
            continue
        _, idx, seq, _ = line.split("\t")
        vertices[idx] = line

        if seq == "*":
            to_remove.add(idx)
            predecessors[idx] = []
            successors[idx] = []
        else:
            print(line, end="")

    all_edges = defaultdict(list)
    not_sinks = set()
    not_sources = set()

    for line in open(gfa_path):
        if not line.startswith("L"):
            continue
        _, idx1, _, idx2, _, _ = line.split("\t")
        all_edges[idx1].append(idx2)

        if idx1 in to_remove:
            successors[idx1].append(idx2)
        elif idx2 in to_remove:
            predecessors[idx2].append(idx1)
        else:
            not_sinks.add(idx1)
            not_sources.add(idx2)
            print(line, end="")

    for idx in successors:
        while any([x in to_remove for x in successors[idx]]):
            new_successors = []
            for sidx in successors[idx]:
                if sidx in to_remove:
                    new_successors.extend(successors[sidx])
                else:
                    new_successors.append(sidx)
            successors[idx] = new_successors
    for idx in predecessors:
        while any([x in to_remove for x in predecessors[idx]]):
            new_predecessors = []
            for sidx in predecessors[idx]:
                if sidx in to_remove:
                    new_predecessors.extend(predecessors[sidx])
                else:
                    new_predecessors.append(sidx)
            predecessors[idx] = new_predecessors

    # adds new edges
    next_idx = 10000
    for idx in to_remove:
        if len(predecessors[idx]) == 0:
            # Vertex idx was a source, add all successors as new sources
            next_idx = print_new_sources(
                idx, vertices, all_edges, not_sources, to_remove, next_idx
            )
            continue

        if len(successors[idx]) == 0:
            continue
        for p in predecessors[idx]:
            for s in successors[idx]:
                print("L", p, "+", s, "+", "0M", sep="\t")

    # for line in gfa_path:
    #     if line.startswith("P"):
    #         _, idx1, _, idx2, _, _ = line.split("\t")
    #         if idx1 in to_remove or idx2 in to_remove:
    #             pass


if __name__ == "__main__":
    main()
