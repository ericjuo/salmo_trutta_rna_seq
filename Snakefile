# RNA-Seq analysis pipeline
#
# Written by Eric Juo
#
# Last modified April 9, 2022
# First written April 7, 2022


configfile: "config.yaml"

rule all:
    input:
        "report/initial_multiqc/raw_read_multiqc_report.html",
        "report/post_trim_multiqc/trimmomatic_trim_multiqc_report.html",
        expand("report/pre_fastq_screen/{sample}_{replicate}_paired_trimmomatic_screen.png", sample=config["samples"], replicate=[1,2]),
        expand("report/post_fastq_screen/{sample}_unclassified_{replicate}_screen.png", sample=config["samples"], replicate=[1,2]),
        "data/03_processed/trinity/brown_trout_trinity_rna_seq.fasta"

# Convert SRA into paired-end fastq files
rule fastq_dump:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        temp("data/02_intermediate/{sample}_1.fastq"),
        temp("data/02_intermediate/{sample}_2.fastq")
    params:
        outdir="data/02_intermediate/"
    shell:
        "fastq-dump -O {params.outdir} --split-3 {input}"

# Compress fastq files in gzip format to save storage space
rule gzip_compress:
    input:
        "data/02_intermediate/{sample}_{direction}.fastq"
    output:
        "data/02_intermediate/gzip_compress/{sample}_{direction}.fastq.gz"
    shell:
        "gzip -c {input} > {output}"

# Report quality of raw reads using FastQC software
rule initial_fastqc:
    input:
        expand("data/02_intermediate/gzip_compress/{sample}_{replicate}.fastq.gz", sample=config["samples"], replicate=[1,2])
    output:
        "report/initial_multiqc/raw_read_multiqc_report.html"
    params:
        init_qc_dir='data/02_intermediate/inital_fastqc/'
    threads: 8
    shell:
        "mkdir -p {params.init_qc_dir} && fastqc -t 8 -o {params.init_qc_dir} -f fastq {input} && multiqc -f -n {output} {params.init_qc_dir}"

# Trim off bases with quality scores lower than 20 using Trimmomatic software. Residual adapters are removed in this step as well.
rule trimmomatic:
    input:
        "data/02_intermediate/gzip_compress/{sample}_1.fastq.gz",
        "data/02_intermediate/gzip_compress/{sample}_2.fastq.gz"
    output:
        "data/02_intermediate/trimmomatic_trimmed/{sample}_1_paired_trimmomatic.fastq.gz",
        "data/02_intermediate/trimmomatic_trimmed/{sample}_1_unpaired_trimmomatic.fastq.gz",
        "data/02_intermediate/trimmomatic_trimmed/{sample}_2_paired_trimmomatic.fastq.gz",
        "data/02_intermediate/trimmomatic_trimmed/{sample}_2_unpaired_trimmomatic.fastq.gz",
    params:
        adapters=config['adapter']
    log: 
        "logs/trimmomatic/{sample}.log",
        "logs/trimmomatic/{sample}_summary.log",
    threads: 10
    shell:
        "trimmomatic PE -threads {threads} -phred33 -trimlog {log[0]} {input[0]} {input[1]} " \
        "{output[0]} {output[1]} {output[2]} {output[3]} " \
        "ILLUMINACLIP:{params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:35 2>{log[1]}"

# Report quality of trimmed reads from trimmomatic software using FastQC
rule post_trim_fastqc:
    input:
        expand("data/02_intermediate/trimmomatic_trimmed/{sample}_{replicate}_paired_trimmomatic.fastq.gz", sample=config["samples"], replicate=[1,2])
    output:
        "report/post_trim_multiqc/trimmomatic_trim_multiqc_report.html"
    params:
        post_trim_qc_dir='data/02_intermediate/trimmomatic_trimmed_fastqc/'
    threads: 8
    shell:
        "mkdir -p {params.post_trim_qc_dir} && fastqc -t 8 -o {params.post_trim_qc_dir} -f fastq {input} && multiqc -f -n {output} {params.post_trim_qc_dir}"

# Report foreign contaminants (Human, Mouse, Rat, Fly, Worm, Yeast, Plant, E.Coli, rRNA
# Mictocondria, PhiX, Lambda, Vectors and adapters) in the trimmed reads using fastq_screen software.
rule pre_fastq_screen:
    input:
        expand("data/02_intermediate/trimmomatic_trimmed/{sample}_{replicate}_paired_trimmomatic.fastq.gz", sample=config["samples"], replicate=[1,2])
    output:
        expand("report/pre_fastq_screen/{sample}_{replicate}_paired_trimmomatic_screen.png", sample=config["samples"], replicate=[1,2])
    params:
        db_conf = "contaminants/FastQ_Screen_Genomes/fastq_screen.conf",
        pre_dir = "report/pre_fastq_screen/"
    threads: 10
    run:
        for f in input:
            shell("fastq_screen -threads {threads} -conf {params.db_conf} --aligner 'bowtie2' --force --outdir {params.pre_dir} {f}")

# Remove foreign contaminants (Human, Archea, Bacteria, Viral) in reads using kraken2
rule kraken2:
    input:
        "data/02_intermediate/trimmomatic_trimmed/{sample}_1_paired_trimmomatic.fastq.gz",
        "data/02_intermediate/trimmomatic_trimmed/{sample}_2_paired_trimmomatic.fastq.gz"
    output:
        "data/02_intermediate/filter_foreign_contaminants/{sample}_classified_1.fastq",
        "data/02_intermediate/filter_foreign_contaminants/{sample}_classified_2.fastq",
        "data/02_intermediate/filter_foreign_contaminants/{sample}_unclassified_1.fastq",
        "data/02_intermediate/filter_foreign_contaminants/{sample}_unclassified_2.fastq"
    params:
        minikraken="contaminants/minikraken2",
        classified="data/02_intermediate/filter_foreign_contaminants/{sample}_classified#.fastq",
        unclassified="data/02_intermediate/filter_foreign_contaminants/{sample}_unclassified#.fastq"
    log: 
        "logs/kraken2/{sample}_kraken2_summary.log",
        "logs/kraken2/{sample}_kraken2_stdout.log"
    threads: 10
    shell:
        "kraken2 -db {params.minikraken} --threads {threads} --classified-out {params.classified} --unclassified-out {params.unclassified} --paired --use-names --gzip-compressed --report {log[0]} " \ 
        "{input[0]} {input[1]} > {log[1]}"

# Report residul contaminants using fastq_screen software.
rule post_fastq_screen:
    input:
        expand("data/02_intermediate/filter_foreign_contaminants/{sample}_unclassified_{replicate}.fastq", sample=config["samples"], replicate=[1,2])
    output:
        expand("report/post_fastq_screen/{sample}_unclassified_{replicate}_screen.png", sample=config["samples"], replicate=[1,2])
    params:
        db_conf = "contaminants/FastQ_Screen_Genomes/fastq_screen.conf",
        post_dir = "report/post_fastq_screen/"
    threads: 10
    run:
        for f in input:
            shell("fastq_screen -threads {threads} -conf {params.db_conf} --aligner 'bowtie2' --force --outdir {params.post_dir} {f}")

# Reformatting fastq file to be feed into Trinity de novo assembler using bioawk.
# The sequence name (e.g. @SRR799769.1) in the fastq file is not recognized by Trinity.
# The /1 or /2 is append to the end of original sequence name (e.g. @SRR799769.1/1) using bioawk
rule pre_trinity_formatting:
    input:
        "data/02_intermediate/filter_foreign_contaminants/{sample}_unclassified_{replicate}.fastq"
    output:
        "data/02_intermediate/pre_trinity_formatting/{sample}_pre_trinity_{replicate}.fastq"
    shell:
        "bioawk -c fastx '{{print \"@\"$name\"/{wildcards.replicate}\\n\"$seq\"\\n+\\n\"$qual}}' {input} > {output}"



# Trinity is performed using snakemake wrapper downloaded from 
# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/trinity.html
rule trinity:
    input:
        left=[
            "data/02_intermediate/pre_trinity_formatting/SRR799769_pre_trinity_1.fastq",
            "data/02_intermediate/pre_trinity_formatting/SRR799770_pre_trinity_1.fastq",
            "data/02_intermediate/pre_trinity_formatting/SRR799771_pre_trinity_1.fastq",
            "data/02_intermediate/pre_trinity_formatting/SRR799772_pre_trinity_1.fastq",
            "data/02_intermediate/pre_trinity_formatting/SRR799773_pre_trinity_1.fastq",
            "data/02_intermediate/pre_trinity_formatting/SRR799774_pre_trinity_1.fastq",
            "data/02_intermediate/pre_trinity_formatting/SRR799775_pre_trinity_1.fastq",
            "data/02_intermediate/pre_trinity_formatting/SRR799776_pre_trinity_1.fastq"
        ],
        right=[
            "data/02_intermediate/pre_trinity_formatting/SRR799769_pre_trinity_2.fastq",
            "data/02_intermediate/pre_trinity_formatting/SRR799770_pre_trinity_2.fastq",
            "data/02_intermediate/pre_trinity_formatting/SRR799771_pre_trinity_2.fastq",
            "data/02_intermediate/pre_trinity_formatting/SRR799772_pre_trinity_2.fastq",
            "data/02_intermediate/pre_trinity_formatting/SRR799773_pre_trinity_2.fastq",
            "data/02_intermediate/pre_trinity_formatting/SRR799774_pre_trinity_2.fastq",
            "data/02_intermediate/pre_trinity_formatting/SRR799775_pre_trinity_2.fastq",
            "data/02_intermediate/pre_trinity_formatting/SRR799776_pre_trinity_2.fastq"          
        ]
    output:
        "data/03_processed/trinity/brown_trout_trinity_rna_seq.fasta"
    log:
        'logs/trinity/trinity_assembly.log'
    params:
        extra="",
        seqtype="fq"
    threads: 10
    resources:
        mem_gb=10
    wrapper:
        "v1.3.2/bio/trinity"