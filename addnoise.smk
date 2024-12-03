from os.path import join as pjoin
import glob

SEQSDIR = config["seqsdir"]

seqs_to_generate = []
for fa in glob.glob(pjoin(SEQSDIR, "*", "MPCSIM", "*", "sequence.fa")):
    for p in ["0", "1", "3", "5", "10", "20"]:
        seqs_to_generate.append(
            pjoin("/".join(fa.split("/")[:-1]), f"p{p}", "sequence.fa")
        )


rule run:
    input:
        seqs_to_generate,


rule add_noise:
    input:
        fa=pjoin(SEQSDIR, "{gene}", "MPCSIM", "{path}", "sequence.fa"),
    output:
        fa=pjoin(SEQSDIR, "{gene}", "MPCSIM", "{path}", "p{p}", "sequence.fa"),
    shell:
        """
        python3 scripts/randsnp.py {input.fa} {wildcards.p} > {output.fa}
        """
