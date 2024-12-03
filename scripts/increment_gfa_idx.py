import sys


def main():
    for line in open(sys.argv[1]):
        if line.startswith("S"):
            t, idx, *rest = line.strip("\n").split("\t")
            idx = str(int(idx) + 1)
            print("\t".join([t, idx] + rest))
        elif line.startswith("L"):
            t, idx1, s1, idx2, s2, c = line.strip("\n").split("\t")
            idx1 = str(int(idx1) + 1)
            idx2 = str(int(idx2) + 1)
            print(t, idx1, s1, idx2, s2, c, sep="\t")
        else:
            print(line, end="")


if __name__ == "__main__":
    main()
