import argparse
import numpy as np
from pprint import pprint

def get_paths(graph):
    paths = list()
    for line in open(graph):
        line = line.strip()
        if line.startswith('P'):
            _, pid, nodes, _ = line.split('\t')
            nodes = nodes[:-1].split('+,')
            paths.append(list(map(int, nodes)))
    return paths

def get_gaf_path(gaf):
    for line in open(gaf):
        line = line.strip()
        p = line.split('\t')[5]
        p = p[1:].split('>')
        return list(map(int, p))

def get_min_recs(graph, path, debug=False):
    graphpaths = get_paths(graph)
    inpath = get_gaf_path(path)

    D = np.zeros((len(inpath), len(graphpaths)))

    for ix, j in enumerate(graphpaths):
        if inpath[0] in j:
            D[0, ix] = 0
        else:
            D[0, ix] = np.inf
    if debug:
        print(inpath)
        print(D[0,:])
    for i in range(1, len(inpath)):
        for j in range(len(graphpaths)):
            if not inpath[i] in graphpaths[j]:
                D[i,j] = np.inf
                continue
            opt1 = np.inf
            if inpath[i-1] in graphpaths[j]:
                opt1 = D[i-1, j]

            minh = np.inf
            for h in range(0, len(graphpaths)):
                x = np.inf
                if inpath[i-1] in graphpaths[h]:
                    x = D[i-1, h] + 1
                if x < minh:
                    minh = x
            D[i,j] = min(opt1, minh)
        if debug:
            print(D[i,:])
    # print(D)
    recmin = int(min(D[-1,:]))
    return recmin
            
def main(args):
    recmin = get_min_recs(args.GRAPH, args.PATH, debug=True)
    print(recmin)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter')
    parser.add_argument('GRAPH', type=str, help="GFA file")
    parser.add_argument('PATH', type=str, help="GAF file")

    args = parser.parse_args()
    main(args)