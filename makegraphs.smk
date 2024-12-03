from os.path import join as pjoin
import glob

SEQSDIR = config["seqsdir"]
RECGRAPH_BIN = (
    config["recgraph"]
    if "recgraph" in config
    else "/data/RecGraph/target/release/recgraph"
)
MAKEPRG_BIN = (
    config["mkprg"] if "mkprg" in config else "make_prg"
)

SEQS = []
for fa in glob.glob(pjoin(SEQSDIR, "*.fa")):
    SEQS.append(fa.split("/")[-1][:-3])


rule run:
    input:
        expand(pjoin(SEQSDIR, "{seq}", "MPCSIM", "mosaics.txt"), seq=SEQS),


rule make_graph:
    input:
        msa=pjoin(SEQSDIR, "{seq}.fa"),
    output:
        gfa=pjoin(SEQSDIR, "{seq}", "graph.prg.gfa"),
    log:
        pjoin(SEQSDIR, "{seq}", "graph.prg.log"),
    params:
        prefix=pjoin(SEQSDIR, "{seq}", "graph"),
    shell:
        """
        {MAKEPRG_BIN} from_msa -F -i {input.msa} -o {params.prefix} --output-type g -t {threads} -v --log {log}
        """


rule clean_prg_graph:
    input:
        gfa=pjoin(SEQSDIR, "{seq}", "graph.prg.gfa"),
    output:
        gfa=pjoin(SEQSDIR, "{seq}", "graph.gfa"),
    shell:
        """
        python3 ./scripts/clean_gfa_from_ast.py {input.gfa} > {output.gfa}
        """


rule tsort_graph:
    input:
        gfa=pjoin(SEQSDIR, "{seq}", "graph.gfa"),
    output:
        gfa=pjoin(SEQSDIR, "{seq}", "graph.tsorted.gfa"),
    params:
        prefix=pjoin(SEQSDIR, "{seq}", "graph"),
    conda:
        "envs/vg.yaml"
    shell:
        """
        python3 ./scripts/increment_gfa_idx.py {input.gfa} > {params.prefix}.incr.gfa
        odgi build -g {params.prefix}.incr.gfa -o {params.prefix}.incr.gfa.og
        odgi sort -i {params.prefix}.incr.gfa.og -o {params.prefix}.incr.tsorted.gfa.og
        odgi view -i {params.prefix}.incr.tsorted.gfa.og -g > {output.gfa}
        """


# rule check_acyclicity:
#     input:
#         gfa=pjoin(SEQSDIR, "{seq}", "graph.tsorted.gfa"),
#     output:
#         log=pjoin(SEQSDIR, "{seq}", "graph.gfa.ACYFLAG"),
#     conda:
#         "envs/vg.yaml"
#     shell:
#         """
#         vg stats -A {input.gfa} > {output.log}
#         """


rule remove_indels:
    input:
        msa=pjoin(SEQSDIR, "{seq}.fa"),
    output:
        fa=pjoin(SEQSDIR, "{seq}", "seqs.noindels.fa"),
    shell:
        """
        python3 ./scripts/rm_indel_msa.py {input.msa} > {output.fa}
        """


rule rspoa_align_1:
    input:
        fa=pjoin(SEQSDIR, "{seq}", "seqs.noindels.fa"),
        gfa=pjoin(SEQSDIR, "{seq}", "graph.tsorted.gfa"),
    output:
        gaf=pjoin(SEQSDIR, "{seq}", "rspoa.gaf"),
    log:
        log=pjoin(SEQSDIR, "{seq}", "rspoa.log"),
    threads: 1
    shell:
        """
        {RECGRAPH_BIN} -m0 -f1 {input.fa} {input.gfa} > {output.gaf} 2> {log.log}
        """


rule add_paths:
    input:
        gfa=pjoin(SEQSDIR, "{seq}", "graph.tsorted.gfa"),
        gaf=pjoin(SEQSDIR, "{seq}", "rspoa.gaf"),
    output:
        gfa=pjoin(SEQSDIR, "{seq}", "graph.tsorted.wpaths.gfa"),
    shell:
        """
        cp {input.gfa} {output.gfa}
        cut -f 1,6 {input.gaf} | while read idx p ; do echo -e "P\\t$idx\\t$(echo $p | cut -c 2- | sed "s/>/+,/g")+\\t*" ; done >> {output.gfa}
        """


rule remove_nodes_not_in_a_path:
    input:
        gfa=pjoin(SEQSDIR, "{seq}", "graph.tsorted.wpaths.gfa"),
    output:
        gfa=pjoin(SEQSDIR, "{seq}", "graph.tsorted.wpaths.clean.gfa"),
    shell:
        """
        python3 ./scripts/remove_nodes_not_in_a_path.py {input.gfa} > {output.gfa}
        """


rule minimal_path_cover:
    input:
        gfa=pjoin(SEQSDIR, "{seq}", "graph.tsorted.wpaths.clean.gfa"),
    output:
        gfa=pjoin(SEQSDIR, "{seq}", "MPCSIM", "graph.gfa"),
        recombinations=pjoin(SEQSDIR, "{seq}", "MPCSIM", "recombinations.log"),
    shell:
        """
        cat {input.gfa} | python3 ./scripts/minimal_path_cover.py > {output.gfa} 2> {output.recombinations}
        """


checkpoint get_mosaics:
    input:
        log=pjoin(SEQSDIR, "{seq}", "MPCSIM", "recombinations.log"),
    output:
        txt=pjoin(SEQSDIR, "{seq}", "MPCSIM", "mosaics.txt"),
    shell:
        """
        grep "is a mosaic" {input.log} | cut -f 3 -d' ' | sort > {output.txt} || true
        """
