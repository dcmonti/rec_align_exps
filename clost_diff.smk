import os.path
import re


rule all:
    input:
        expand(
        "output/cdifficile/full.csv"
        )
    
    conda: "envs/csvkit.yaml"
    shell:
        """
        MYFI=""
        MYGR=""
        for f in {input:q}; do MYFI="$MYFI $f"; MYGR="$MYGR,${{f#output\/cdifficile\/}}"; done
        csvstack -g ${{MYGR#,}} $MYFI > {output}
        """

rule build_recgraph_a_star:
    output:
        "bin/recgraph_a_star"
    shadow: "shallow"
    conda: "envs/rust.yaml"
    threads: 4
    shell:
        """
        mkdir a_star
        cd a_star
        git clone https://github.com/AlgoLab/RecGraph.git
        cd RecGraph/
        git checkout a_star
        cargo build --release --jobs {threads}
        cd ../..
        cp a_star/RecGraph/target/release/recalign {output}
        """

rule make_graph:
    output:
        gfa_folder= "output/cdifficile/gfa",
        gfa_file = "output/cdifficile/gfa/graph.gfa",
    conda: "envs/make_graph.yaml"
    threads: 4
    shadow: "shallow"
    shell:
        """
        pggb -i data/cdifficile/slpa-basis.fa -o {output.gfa_folder} -n 11 -s 100 -p 50 -K 8 -k 1 -B 10 -G 2750,4500 -P asm5 -d 12 -O 0.5
        cp output/cdifficile/gfa/slpa-basis.fa.3f0681a.b682f59.ebb3f5c.smooth.fix.gfa {output.gfa_file}
        """

rule run_recgraph_a_star:
    input:
        recgraph="bin/recgraph_a_star",
        fa = "output/cdifficile/split/simulated/split.{seq}.fa",
        gfa = "output/cdifficile/gfa/graph.gfa",
    output:
        gaf = "output/cdifficile/recgraph/split.{seq}.gaf",
    threads: 1
    resources:
        mem="12GB"
    shell:
        """
        {input.recgraph} \
            -s 8 -r 4 -k 1 \
            -q {input.fa} -g {input.gfa} > {output.gaf}
        """  

def aggregate_recgraph_input(wildcards):
    checkpoint_output = checkpoints.split_seqs.get(**wildcards).output[0]
    ret = expand(
        "output/cdifficile/recgraph/split.{seq}.gaf",
        **wildcards,
        seq=glob_wildcards(os.path.join(checkpoint_output, "split.{seq}.fa")).seq
    )
    return sorted(ret)

rule merge_recgraph:
    input:
        aggregate_recgraph_input
    output:
        "output/cdifficile/recgraph/full.gaf",
    threads: 1
    shell:
        """
        cat {input} > {output}
        """

rule make_csv:
    input:
        gaf = "output/cdifficile/recgraph/full.gaf",
    output:
        csv = "output/cdifficile/full.csv",
    run:
        out = open(output.csv, 'w+')
        print("ReadName,RecPaths,RecPos,Score,Displacement", file=out)
        orig_lines = iter(open(input.gaf).readlines()) 
        for line,comment in zip(orig_lines, orig_lines):
            line = line.strip().split('\t')
            comment = comment.strip().split('\t')
            rid=line[0]

            paths = comment[3].split(',')
            if len(paths) == 1:
                path = paths[0].split(':')
                recpath = f"{int(path[0])+1}"
                bp = ''
            else:
                path1 = paths[0].split(':')
                path2 = paths[1].split(':')
                recpath = f"{int(path1[0])+1}>{int(path2[0])+1}"
                bp = path1[2]
            
            score=f"{comment[1]}"
            
            
            displacement=''

            print(rid, recpath, bp, score, displacement, sep=',', file=out)
        out.close()

checkpoint split_seqs:
    input:
        fa = "data/cdifficile/simulated.fa",
    output:
        dir = directory("output/cdifficile/split/simulated"),
    shell:
        """
        mkdir -p {output.dir}
        awk 'BEGIN {{ NUM = -1 }} /^>/ {{NUM = NUM+1; F=sprintf("{output.dir}/split.%.4d.fa", NUM); print > F; next; }} {{print >> F;}}' < {input.fa}
        """