rule all:
    input:
        "output/sars-cov-2/graph.gfa",
        "bin/pbsim",
        "bin/QSHMM-ONT-HQ.model",
        "bin/ERRHMM-SEQUEL.model",
        expand("output/sars-cov-2/reads/{haplo}/reads.fasta", haplo=glob_wildcards("data/sars-cov-2/haplotypes/sequence_{haplo}.fasta").haplo),
        expand("output/sars-cov-2/rec_reads/{rec_haplo}/reads.fasta", rec_haplo=glob_wildcards("data/sars-cov-2/rec_haplos/sequence_r{rec_haplo}.fasta").rec_haplo),
    shell:
        """
        bash split_sars_fasta.sh
        echo "Reads simulated and graph built"
        """

rule make_graph:
    params:
        gfa_folder = "output/sars-cov-2/gfa",
    output:
        gfa_file = "output/sars-cov-2/graph.gfa",
        haplos = "output/sars-cov-2/haplos.fa",
        rec_haplos = "output/sars-cov-2/rec_haplos.fa"
    conda: "envs/make_graph.yaml"
    threads: 4
    shell:
        """
        rm -rf output/sars-cov-2/haplos.fa
        rm -rf output/sars-cov-2/rec_haplos.fa
        cat data/sars-cov-2/haplotypes/*.fasta >> {output.haplos}
        cat data/sars-cov-2/rec_haplos/*.fasta >> {output.rec_haplos}
        pggb -i output/sars-cov-2/haplos.fa -o {params.gfa_folder} -n 11 -s 100 -p 50 -K 8 -k 1 -B 10 -G 2750,4500 -P asm5 -d 12 -O 0.5
        cp output/sars-cov-2/gfa/haplos.fa.3f0681a.b682f59.ebb3f5c.smooth.fix.gfa {output.gfa_file}
        rm -rf output/sars-cov-2/gfa/
        """

rule install_pbsim:
    output:
        pbsim = "bin/pbsim",
        method ="bin/QSHMM-ONT-HQ.model",
        err = "bin/ERRHMM-SEQUEL.model"
    shell:
        """
        rm -rf pbsim3
        git clone https://github.com/yukiteruono/pbsim3
        cd pbsim3
         autoreconf -i
        ./configure
        make
        cd ..
        cp pbsim3/src/pbsim {output.pbsim}
        cp pbsim3/data/QSHMM-ONT-HQ.model {output.method}
        cp pbsim3/data/ERRHMM-SEQUEL.model {output.err}
        rm -rf pbsim3
        """

rule generate_reads:
    input:
        haplos="output/sars-cov-2/haplos.fa",
        pbsim="bin/pbsim"
    output:
        "output/sars-cov-2/sim/sd_00{haplo}.fastq",
    shell:
        """
        {input.pbsim} --strategy wgs --method qshmm --qshmm bin/QSHMM-ONT-HQ.model --length-min 26000 --length-max 28500 --length-mean 27000  --accuracy-mean .99 --hp-del-bias 10 --depth 100 --genome {input.haplos} --prefix output/sars-cov-2/sim/sd --difference-ratio 2:1:1
        """

rule generate_rec_reads:
    input:
        haplos="output/sars-cov-2/rec_haplos.fa",
        pbsim="bin/pbsim"
    output:
        "output/sars-cov-2/sim_rec/sd_00{rec_haplo}.fastq",
    shell:
        """
        {input.pbsim} --strategy wgs --method qshmm --qshmm bin/QSHMM-ONT-HQ.model --length-min 26000 --length-max 28500 --length-mean 27000  --accuracy-mean .99 --hp-del-bias 10 --depth 100 --genome {input.haplos} --prefix output/sars-cov-2/sim_rec/sd --difference-ratio 2:1:1
        """

rule filter_reads:
    input:
        reads="output/sars-cov-2/sim/sd_00{haplo}.fastq",
        pbsim="bin/pbsim",
        haplos="data/sars-cov-2/haplotypes/sequence_{haplo}.fasta"
    output:
        "output/sars-cov-2/sim/filtered_sd_{haplo}_0001.fastq",
    shell:
        """ 
        {input.pbsim} --strategy wgs --method sample --sample {input.reads} --genome {input.haplos} --accuracy-min 0.96 --depth 80 --prefix output/sars-cov-2/sim/filtered_sd_{wildcards.haplo} --id-prefix S{wildcards.haplo}_
        """

rule filter_reads_rec:
    input:
        reads="output/sars-cov-2/sim_rec/sd_00{rec_haplo}.fastq",
        pbsim="bin/pbsim",
        haplos="data/sars-cov-2/rec_haplos/sequence_r{rec_haplo}.fasta"
    output:
        "output/sars-cov-2/sim_rec/filtered_sd_{rec_haplo}_0001.fastq",
    shell:
        """ 
        {input.pbsim} --strategy wgs --method sample --sample {input.reads} --genome {input.haplos} --accuracy-min 0.96 --depth 80 --prefix output/sars-cov-2/sim_rec/filtered_sd_{wildcards.rec_haplo} --id-prefix S{wildcards.rec_haplo}_
        """

rule convert_reads:
    input:
        "output/sars-cov-2/sim/filtered_sd_{haplo}_0001.fastq",
    output:
        "output/sars-cov-2/reads/{haplo}/reads.fasta"
    shell:
        """
        sed -n '1~4s/^@/>/p;2~4p' {input} > {output}
        """

rule convert_reads_rec:
    input:
        "output/sars-cov-2/sim_rec/filtered_sd_{rec_haplo}_0001.fastq",
    output:
        "output/sars-cov-2/rec_reads/{rec_haplo}/reads.fasta"
    shell:
        """
        sed -n '1~4s/^@/>/p;2~4p' {input} > {output}
        """
