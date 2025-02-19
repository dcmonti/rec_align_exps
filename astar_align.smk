import os

genes = []
with open("genes_HLA_short.txt", "r") as f:
    genes = [gene.strip() for gene in f.readlines()]

error_levels = ["0", "3", "5"]
recs= ["0", "1", "2"]

split_read_files = []
for gene in genes:
    for rec in recs:
        for err in error_levels:
            output_dir = f"output/HLA/genes/{gene}/{rec}/reads_{err}_split"
            if os.path.exists(output_dir):
                for read_file in os.listdir(output_dir):
                    if read_file.startswith("read_") and read_file.endswith(".fa"):
                        path = os.path.join(f"{gene}/{rec}/reads_{err}_split", read_file.split(".")[0])
                        split_read_files.append(path)
print(split_read_files)

rule all:
    input:
        "bin/recgraph",
        "bin/recgraph_a_star",
        expand("output/HLA/{mode}_recgraph/{read_path}.gaf",
                mode=["old", "new"], 
                read_path=split_read_files
        )

rule build_recgraph:
    output:
        "bin/recgraph"
    shadow: "shallow"
    conda: "envs/rust.yaml"
    threads: 4
    shell:
        """
        git clone https://github.com/AlgoLab/RecGraph.git
        cd RecGraph/
        git checkout ef684e8fe6a9cd1218f8e418aacf195f3702e3b7
        cargo build --release --jobs {threads}
        cd ..
        cp RecGraph/target/release/recgraph {output}
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
        cp a_star/RecGraph/target/release/recgraph {output}
        """

rule run_recgraph:
    input:
        fa="output/HLA/genes/{gene}/{rec}/reads_{err}_split/{read_file}.fa",
        gfa="output/HLA/genes/{gene}/graph.gfa",
        rg = "bin/recgraph"
    log:
        "output/HLA/old_recgraph/{gene}/{rec}/reads_{err}_split/{read_file}.log"
    output:
        "output/HLA/old_recgraph/{gene}/{rec}/reads_{err}_split/{read_file}.gaf"
    shell:
        "/usr/bin/time -v {input.rg} {input.fa} {input.gfa} -m 8 > {output} 2> {log}"

rule run_recgraph_a_star:
    input:
        fa="output/HLA/genes/{gene}/{rec}/reads_{err}_split/{read_file}.fa",
        gfa="output/HLA/genes/{gene}/graph.gfa",
        rg = "bin/recgraph_a_star"
    log:
        "output/HLA/new_recgraph/{gene}/{rec}/reads_{err}_split/{read_file}.log"
    output:
        "output/HLA/new_recgraph/{gene}/{rec}/reads_{err}_split/{read_file}.gaf"
    shell:
        "/usr/bin/time -v {input.rg} -q {input.fa} -g {input.gfa} -s 8 -r 4 -k 2 > {output} 2> {log}"
