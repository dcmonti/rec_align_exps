import os

split_read_files = []
output_dir = "output/sars-cov-2/reads/"
if os.path.exists(output_dir):
    for read_file in os.listdir(output_dir):
        if read_file.endswith(".fasta"):
            path = read_file.split(".")[0]
            split_read_files.append(path)

split_rec_read_files = []
output_dir = "output/sars-cov-2/rec_reads/"
if os.path.exists(output_dir):
    for read_file in os.listdir(output_dir):
        if read_file.endswith(".fasta"):
            path = read_file.split(".")[0]
            split_rec_read_files.append(path)


rule all:
    input:
        "bin/recgraph_a_star",
        "bin/minichain",
        expand("output/sars-cov-2/ori/{mode}/{read_path}.gaf",
                mode=["ga", "ra", "mc"], 
                read_path=split_read_files
        ),
        expand("output/sars-cov-2/rec/{mode}/{rec_read_path}.gaf",
                mode=["ga", "ra", "mc"], 
                rec_read_path=split_rec_read_files
        )

rule build_minichain:
    output:
        "bin/minichain"
    shadow: "shallow"
    #conda: "envs/rust.yaml"
    threads: 4
    shell:
        """
        mkdir minichain
        cd minichain
        git clone https://github.com/at-cg/minichain
        cd minichain && make
        cd ../..
        cp minichain/minichain/minichain {output}
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
rule generate_mc_graph:
    input:
        gfa="output/sars-cov-2/graph.gfa"
    output:
        w_gfa="output/sars-cov-2/w_graph.gfa"
    shell:
        "python scripts/mc_fix.py {input.gfa} > {output.w_gfa} "

rule run_recgraph_a_star:
    input:
        fa="output/sars-cov-2/reads/{read_file}.fasta",
        gfa="output/sars-cov-2/graph.gfa",
        rg = "bin/recgraph_a_star"
    log:
        "output/sars-cov-2/ori/ra/{read_file}.log"
    output:
        "output/sars-cov-2/ori/ra/{read_file}.gaf"
    shell:
        """
        touch {output}
        touch {log}
        timeout 2m /usr/bin/time -v {input.rg} -q {input.fa} -g {input.gfa} -s 14 -k 0 -m -e fast > {output} 2> {log}
        """

rule run_minichain:
    input:
        fa="output/sars-cov-2/reads/{read_file}.fasta",
        gfa="output/sars-cov-2/w_graph.gfa",
        mc = "bin/minichain"
    log:
        "output/sars-cov-2/ori/mc/{read_file}.log"
    output:
        "output/sars-cov-2/ori/mc/{read_file}.gaf"
    shell:
        "/usr/bin/time -v bin/minichain -cx lr {input.gfa} {input.fa} -R 4 -k 12 -w 8 > {output} 2> {log}"

rule run_graphaligner:
    input:
        fa="output/sars-cov-2/reads/{read_file}.fasta",
        gfa="output/sars-cov-2/graph.gfa"
    log:
        "output/sars-cov-2/ori/ga/{read_file}.log"
    output:
        "output/sars-cov-2/ori/ga/{read_file}.gaf"
    shell:
        "/usr/bin/time -v GraphAligner -f {input.fa}  -g {input.gfa} -a {output} -x vg 2> {log}"


rule run_recgraph_a_star_r:
    input:
        fa="output/sars-cov-2/rec_reads/{read_file}.fasta",
        gfa="output/sars-cov-2/graph.gfa",
        rg = "bin/recgraph_a_star"
    log:
        "output/sars-cov-2/rec/ra/{read_file}.log"
    output:
        "output/sars-cov-2/rec/ra/{read_file}.gaf"
    shell:
        "timeout 2m /usr/bin/time -v {input.rg} -q {input.fa} -g {input.gfa} -s 10 -r 4 -k 2 -m -e fast > {output} 2> {log}"

rule run_minichain_r:
    input:
        fa="output/sars-cov-2/rec_reads/{read_file}.fasta",
        gfa="output/sars-cov-2/w_graph.gfa",
        mc = "bin/minichain"
    log:
        "output/sars-cov-2/rec/mc/{read_file}.log"
    output:
        "output/sars-cov-2/rec/mc/{read_file}.gaf"
    shell:
        """
        touch {output}
        touch {log}
        /usr/bin/time -v bin/minichain -cx lr {input.gfa} {input.fa} -R 4 -k 12 -w 8 > {output} 2> {log}
        """

rule run_graphaligner_r:
    input:
        fa="output/sars-cov-2/rec_reads/{read_file}.fasta",
        gfa="output/sars-cov-2/graph.gfa"
    log:
        "output/sars-cov-2/rec/ga/{read_file}.log"
    output:
        "output/sars-cov-2/rec/ga/{read_file}.gaf"
    shell:
        "/usr/bin/time -v GraphAligner -f {input.fa}  -g {input.gfa} -a {output} -x vg 2> {log}"