#!/usr/bin/env python3

import collections
import sys

def edges_from_path(path: str):
    vertices = tuple([-2] + [ int(v.rstrip("+-")) for v in path.split(",") ] + [-1])
    return (list(zip(vertices, vertices[1:])), vertices)

path_ids = []   # To keep the input order
paths = {}
path_set = {}
edge_counts = collections.Counter()

vlens = { -1: 0, -2: 0 }

dup_paths = {}

for line in sys.stdin:
    line = line.rstrip()
    if line.startswith('S'):
        prefix, vid, seq = line.split(sep=None)
        vlens[int(vid)] = len(seq)

    if not line.startswith('P'):
        print(line)
        continue
    
    prefix, pid, path, rest = line.split(sep=None, maxsplit=3)

    edges, vertices = edges_from_path(path)
    if vertices in path_set:
        print("# Path", pid, "is a duplicate of path", path_set[vertices], file=sys.stderr)
        dup_paths[pid] = ((prefix, pid, path, rest), edges)
        continue
    path_ids.append(pid)
    path_set[vertices] = pid

    paths[pid] = ((prefix, pid, path, rest), edges)
    edge_counts.update(edges)

to_delete = []
for pid in path_ids:
    edges = paths[pid][1]
    # Check if counters are all gt 1
    if all((edge_counts[edge]>1 for edge in edges)):
        to_delete.append(pid)
        edge_counts.subtract(edges)

to_delete = set(to_delete)
to_keep = paths.keys() - to_delete

print("# Input paths:", len(paths) + len(dup_paths), file=sys.stderr)
print("# Duplicated paths:", len(dup_paths), file=sys.stderr)
print("# Redundant paths:", len(to_delete), file=sys.stderr)
print("# Output paths: ", len(paths) - len(to_delete), file=sys.stderr)

path_id = {
    pid:i+1 for i,pid in enumerate(sorted(to_keep))
}
edge_to_paths = collections.defaultdict(set)
for pid in sorted(to_keep):
    line, edges = paths[pid]
    for edge in edges:
        edge_to_paths[edge].add(pid)
    print(f"# PATH {path_id[pid]} -> {pid}", file=sys.stderr)
    print("\t".join(line))

# Compute the set of paths
candidates_cache = {}
def compute_candidates(edges, ce, min_dist, base_dist = None):
    base_dist = base_dist if base_dist is not None else vlens[edges[ce][0]]
    key = (edges[ce], base_dist)
    if key in candidates_cache:
        return candidates_cache[key]
    cum_len = base_dist

    candidates = set(to_keep)
    for i in range(ce, len(edges)):
        edge = edges[i]
        candidates.intersection_update(edge_to_paths[edge])
        cum_len += vlens[edge[1]]
        if cum_len >= min_dist:
            break
    if cum_len < min_dist:
        candidates_cache[key] = set()
        return set()
    candidates_cache[key] = candidates
    return candidates

def ok_path(edges, ce, curr, path, curr_dist, max_recomb, min_dist):
    cpid = path_id[curr]
    path = path + [cpid]
    if ce >= len(edges):
        if curr_dist < min_dist:
            return None
        return path[1:]
    edge = edges[ce]
    vertex = edge[0]
    curr_dist += vlens[vertex]

    if curr in edge_to_paths[edge]:
        result = ok_path(edges, ce+1, curr, path, curr_dist, max_recomb, min_dist)
        if result is not None:
            return result

    # We have to introduce a recombination
    if curr_dist < min_dist or max_recomb == 0:
        return None
    base_dist = min(curr_dist - min_dist, vlens[vertex])
    candidates = compute_candidates(edges, ce, min_dist) - {curr}
    for candidate in sorted(candidates):
        result = ok_path(
            edges, ce + 1, candidate, path, base_dist, max_recomb - 1, min_dist
        )
        if result is not None:
            return result
    return None


MIN_DIST = 150
MAX_RECOMB = 1
for pid in sorted(to_delete):
    edges = paths[pid][1]

    candidates_cache = {}
    candidates = compute_candidates(edges, 0, MIN_DIST)

    result = None
    for max_recomb in range(1, MAX_RECOMB + 1):
        print(
            f"# Trying to explain path {pid} with {max_recomb} recombinations...",
            file=sys.stderr,
        )
        for candidate in sorted(candidates):
            result = ok_path(edges, 0, candidate, [], 0, max_recomb, MIN_DIST)
            if result is not None:
                break
        if result is not None:
            break
    if result is None:
        print(
            f"# Path {pid} is a NOT a mosaic of the path cover with minimum recombination distance of {MIN_DIST} and max no of recombination of {MAX_RECOMB}",
            file=sys.stderr,
        )
    else:
        print(
            f"# Path {pid} is a mosaic with {max_recomb} recombinations of {list(zip(edges, result))}",
            file=sys.stderr,
        )
        # print(f"# Path {pid} is a mosaic with {max_recomb} recombinations of {result}", file=sys.stderr)
        assert len(result) == len(edges)
