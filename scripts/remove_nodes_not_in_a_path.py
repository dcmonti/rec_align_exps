import sys

paths = []
for line in open(sys.argv[1]):
    if not line.startswith("P"):
        continue
    paths.append([x[:-1] for x in line.strip("\n").split("\t")[2].split(",")])

for line in open(sys.argv[1]):
    if line.startswith("H") or line.startswith("P"):
        print(line, end="")
    elif line.startswith("S"):
        _, idx, *_ = line.strip("\n").split("\t")
        if any([idx in path for path in paths]):
            print(line, end="")
    elif line.startswith("L"):
        _, idx1, _, idx2, _, _ = line.strip("\n").split("\t")
        if any([idx1 in path for path in paths]) and any(
            [idx2 in path for path in paths]
        ):
            print(line, end="")
