# RNA-Seq analysis pipeline
#
# Written by Eric Juo
#
# Last modified April 14, 2022
# First written April 7, 2022


configfile: "config.yaml"

rule all:
    input:
        "report/initial_multiqc/raw_read_multiqc_report.html",
        "report/post_trim_multiqc/trimmomatic_trim_multiqc_report.html",
        expand("report/pre_fastq_screen/{sample}_{replicate}_paired_trimmomatic_screen.png", sample=config["samples"], replicate=[1,2]),
        expand("report/post_fastq_screen/{sample}_unclassified_{replicate}_screen.png", sample=config["samples"], replicate=[1,2]),
        # "trinity_Trinity.fasta"

# Convert SRA into paired-end fastq files
rule fastq_dump:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "data/02_intermediate/{sample}_1.fastq",
        "data/02_intermediate/{sample}_2.fastq"
    params:
        outdir="data/02_intermediate/"
    resources:
        mem_gb=16
    conda:
        "env/sratools.yaml"
    shell:
        "fastq-dump -O {params.outdir} --split-3 {input}"

# Compress fastq files in gzip format to save storage space
rule gzip_compress:
    input:
        "data/02_intermediate/{sample}_{direction}.fastq"
    output:
        "data/02_intermediate/gzip_compress/{sample}_{direction}.fastq.gz"
    resources:
        mem_gb=16
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
    resources:
        mem_gb=16
    conda:
        "env/fastqc.yaml"
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
    resources:
        mem_gb=16
    conda:
        "env/trimmomatic.yaml"
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
    resources:
        mem_gb=16
    conda:
        "env/fastqc.yaml"
    shell:
        "mkdir -p {params.post_trim_qc_dir} && fastqc -t 8 -o {params.post_trim_qc_dir} -f fastq {input} && multiqc -f -n {output} {params.post_trim_qc_dir}"

# Report foreign contaminants (Human, Mouse, Rat, Fly, Worm, Yeast, Plant, E.Coli, rRNA
# Mictocondria, PhiX, Lambda, Vectors and adapters) in the trimmed reads using fastq_screen software.
rule pre_fastq_screen:
    input:
        "data/02_intermediate/trimmomatic_trimmed/{sample}_{replicate}_paired_trimmomatic.fastq.gz"
    output:
        "report/pre_fastq_screen/{sample}_{replicate}_paired_trimmomatic_screen.png"
    params:
        db_conf = "contaminants/FastQ_Screen_Genomes/fastq_screen.conf",
        pre_dir = "report/pre_fastq_screen/"
    threads: 10
    resources:
        mem_gb=16
    conda:
        "env/fastqscreen.yaml"
    shell:
        "fastq_screen -threads {threads} -conf {params.db_conf} --aligner 'bowtie2' --force --outdir {params.pre_dir} {input}"

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
    resources:
        mem_gb=16
    conda:
        "env/kraken2.yaml"
    shell:
        "kraken2 -db {params.minikraken} --threads {threads} --classified-out {params.classified} --unclassified-out {params.unclassified} --paired --use-names --gzip-compressed --report {log[0]} " \ 
        "{input[0]} {input[1]} > {log[1]}"

# Report residul contaminants using fastq_screen software.
rule post_fastq_screen:
    input:
        "data/02_intermediate/filter_foreign_contaminants/{sample}_unclassified_{replicate}.fastq"
    output:
        "report/post_fastq_screen/{sample}_unclassified_{replicate}_screen.png"
    params:
        db_conf = "contaminants/FastQ_Screen_Genomes/fastq_screen.conf",
        post_dir = "report/post_fastq_screen/"
    threads: 10
    resources:
        mem_gb=16
    conda:
        "env/fastqscreen.yaml"
    shell:
        "fastq_screen -threads {threads} -conf {params.db_conf} --aligner 'bowtie2' --force --outdir {params.post_dir} {input}"

# Reformatting fastq file to be feed into Trinity de novo assembler using bioawk.
# The sequence name (e.g. @SRR799769.1) in the fastq file is not recognized by Trinity.
# The /1 or /2 is append to the end of original sequence name (e.g. @SRR799769.1/1) using bioawk
rule pre_trinity_formatting:
    input:
        "data/02_intermediate/filter_foreign_contaminants/{sample}_unclassified_{replicate}.fastq"
    output:
        "data/02_intermediate/pre_trinity_formatting/{sample}_pre_trinity_{replicate}.fastq"
    conda:
        "env/bioawk.yaml"
    resources:
        mem_gb=16
    shell:
        "bioawk -c fastx '{{print \"@\"$name\"/{wildcards.replicate}\\n\"$seq\"\\n+\\n\"$qual}}' {input} > {output}"

rule concatenate_pre_trinity:
    input:
        "data/02_intermediate/pre_trinity_formatting/SRR799769_pre_trinity_{replicate}.fastq",
        "data/02_intermediate/pre_trinity_formatting/SRR799770_pre_trinity_{replicate}.fastq",
        "data/02_intermediate/pre_trinity_formatting/SRR799771_pre_trinity_{replicate}.fastq",
        "data/02_intermediate/pre_trinity_formatting/SRR799772_pre_trinity_{replicate}.fastq",
        "data/02_intermediate/pre_trinity_formatting/SRR799773_pre_trinity_{replicate}.fastq",
        "data/02_intermediate/pre_trinity_formatting/SRR799774_pre_trinity_{replicate}.fastq",
        "data/02_intermediate/pre_trinity_formatting/SRR799775_pre_trinity_{replicate}.fastq",
        "data/02_intermediate/pre_trinity_formatting/SRR799776_pre_trinity_{replicate}.fastq"

    output:
        "data/02_intermediate/concatenated_pre_trinity/pre_trinity_{replicate}.fastq"
    shell:
        "cat {input} > {output}"
# =====================================================================================================
# The following rule is run on google cloud virtual machine e2-highmem-8 (8 vCPU, 64 GB RAM)
# If your machine has at least 8 CPU and 60 GB RAM, you can uncommand the following rule and run it on you local machine.
# singularity was installed using the bash script 'install_singularity.sh'
# Trinity docker was download and converted into singularity image by singularity
# Trinity command is as follow
# singularity exec -e docker://registry.hub.docker.com/trinityrnaseq/trinityrnaseq Trinity --seqType fq --left "concatenated_pre_trinity/pre_trinity_1.fastq" --right "concatenated_pre_trinity/pre_trinity_2.fastq" --CPU 8 --output data/03_processed/trinity/ --max_memory 60G > logs/trinity/run.log 2> logs/trinity/error.log
# ========================================================================================================
# rule trinity:
#     input:
#         left="data/02_intermediate/concatenated_pre_trinity/pre_trinity_1.fastq",
#         right="data/02_intermediate/concatenated_pre_trinity/pre_trinity_2.fastq"       
#     output:
#         "trinity_Trinity.fasta"
#     log:
#         'logs/trinity/trinity_assembly.log'
#     singularity:
#         "docker://registry.hub.docker.com/trinityrnaseq/trinityrnaseq"
#     threads: 8
#     params:
#         seqtype='fq',
#         outdir='data/03_processed/trinity/'
#     resources:
#         mem_gb=60
#     shell:
#         "Trinity --seqType {params.seqtype} --left {input.left} --right {input.right} --CPU {threads} --max_memory {resources.mem_gb}G --output {params.outdir} > logs/trinity/run.log 2> logs/trinity/error.log"
# ============================================================================================================
# End of section
# ============================================================================================================

rule read_representation:
    input:
        "data/03_processed/trinity/trinity.Trinity.fasta",
        "data/02_intermediate/concatenated_pre_trinity/pre_trinity_1.fastq",
        "data/02_intermediate/concatenated_pre_trinity/pre_trinity_2.fastq"
    output:
        "report/read_representation/align_stats.txt",
        "data/02_intermediate/bowtie2/aligned_reads.bam"
    params:
        indexed_basename = "data/03_processed/trinity/brown_trout_trinity.fai",
        max_aligned = 20
    threads: 10
    shell:
        "bowtie2-build {input[0]} {params.indexed_basename} && bowtie2 -q -p {threads} "\
        " --no-unal -k {params.max_aligned} -x {params.indexed_basename} " \
        "-1 {input[1]} -2 {input[2]} 2> {output[0]} | samtools view -@ {threads} -b -o {output[1]}"
