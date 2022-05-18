# De novo assembling using Trinity
After clean up foreign species contaminants, forward reads and reverse reads preprocessed to conform format requirements of Trinity software and separately pooled together. Pooled forward reads and pooled reverse reads were uploaded to google cloud virtual machine for de novo assebmly using trinity.

## Steps
1. Forward and reverse reads were separately pooled together:  
```
$ cat *_Unclassified_1.fastq > Pooled_reads_1.fastq
$ cat *_Unclassified_2.fastq > Pooled_reads_2.fastq
```  


2. marker /1 and /2 were added to the end of sequence name in order to distinguish forward or reverse read. This step is conducted using bioawk software:  
```
$ bioawk -c fastx '{print "@"$name"/1\n"$seq"\n+\n"$qual"}' Pooled_reads_1.fastq > pre_trinity_1.fasta
$ bioawk -c fastx '{print "@"$name"/2\n"$seq"\n+\n"$qual"}' Pooled_reads_2.fastq > pre_trinity_2.fasta
```  

3. Reads are assembled using Trinity container on google cloud virtual machine with singularity software: 
```
singularity exec -e docker://registry.hub.docker.com/trinityrnaseq/trinityrnaseq Trinity --seqType fq --left pre_trinity_1.fasta --right pre_trinity_2.fasta --CPU 8 --max_memory 60G -O trinity/ > run.log 2> error.log
```

## Output
The assembled transcripts are sotred in trinity.Trinity.fasta file.

## Post assemble quality control
### Basic statistic
Basic statistics was caculated using trinity's perl script:
```
singularity exec -e docker://registry.hub.docker.com/trinityrnaseq/trinityrnaseq /usr/local/bin/util/TrinityStats.pl data/03_processed/trinity/trinity.Trinity.fasta > report/trinitystats/basic_stats.txt
```

### Representation of assembled transcripts
To assess the representation of assembled transcripts to input reads, input reads were aligned back to the assembly using bowtie2.
```
# assembled transcripts are indexed first
$ bowtie2-build trinity.Trinity.fasta brown_trout_trinity.fai
# The indexed transcripts are used in the input reads realignment
$ bowtie2 -p 10 -q --no-unal -k 20 -x data/03_processed/trinity/brown_trout_trinity.fai -1 data/03_processed/trinity/brown_trout_trinity.fai -2 data/02_intermediate/concatenated_pre_trinity/pre_trinity_2.fastq 2> report/read_representation/align_stats.txt | samtools view -@ 10 -b data/02_intermediate/bowtie2/aligned_reads.bam
```
The statistic of realignment showed that 89.91% of reads can align to the assembled transcripts, indicating good quality of this assembly.  
```
$ cat align-stats.txt
23453284 reads; of these:
  23453284 (100.00%) were paired; of these:
    3880194 (16.54%) aligned concordantly 0 times
    4959861 (21.15%) aligned concordantly exactly 1 time
    14613229 (62.31%) aligned concordantly >1 times
    ----
    3880194 pairs aligned concordantly 0 times; of these:
      161179 (4.15%) aligned discordantly 1 time
    ----
    3719015 pairs aligned 0 times concordantly or discordantly; of these:
      7438030 mates make up the pairs; of these:
        4731395 (63.61%) aligned 0 times
        869484 (11.69%) aligned exactly 1 time
        1837151 (24.70%) aligned >1 times
89.91% overall alignment rate
```
### Blast against brown_trout cDNA database
Salmo trutta cDNA database is downloaded from ensembl:
```
$ wget http://ftp.ensembl.org/pub/release-106/fasta/salmo_trutta/cdna/Salmo_trutta.fSalTru1.1.cdna.all.fa.gz
```
Perform blastn against borwn trout cDNA database:
```
singularity exec -e docker://registry.hub.docker.com/trinityrnaseq/trinityrnaseq blastn -query data/03_processed/trinity/trinity.Trinity.fasta -db contaminants/brown_trout_cDNA/brown_trout -out data/02_intermediate/blastn/blastn_1e-20.outfmt6 -evalue 1e-20 -num_threads 10 -max_target_seqs 1 -outfmt 6
```
Caculate percent coverage of transcripts to cDNA database using trinity's perl script
```
gunzip -c contaminants/brown_trout_cDNA/Salmo_trutta.fSalTru1.1.cdna.all.fa.gz | singularity exec -e docker://registry.hub.docker.com/trinityrnaseq/trinityrnaseq /usr/local/bin/util/analyze_blastPlus_topHit_coverage.pl data/02_intermediate/blastn/blastn_1e-20.outfmt6 data/03_processed/trinity/trinity.Trinity.fasta -
```
Draw brown trout cDNA coverage using R
```
trout <- read.table("data/02_intermediate/blastn/blastn_1e-20.outfmt6.hist")
names(trout) <- c("pct_cov_bin", "count_in_bin", "accumulated_counts")
jpeg(file="report/trinitystats/full_len_brown_trout_cDNA.jpg")
barplot(trout$count_in_bin, names.arg=trout$pct_cov_bin, xlab="Percent coverage", ylab="Transcript counts", main="Brown Trout cDNA coverage")
dev.off()
```

### Assess completness of transcript assembly using BUSCO
Install
```
mamba create -n busco -c bioconda busco=5.3.2
```
Use BUSCO genes from actinopterygii to search against trinity assembly
```
busco -c 10 -i data/03_processed/trinity/trinity.Trinity.fasta -l actinopterygii_odb10 -o data/02_intermediate/busco -m transcriptome
```
Result in the short_summary.specific.actinopterygii_odb10.busco.txt file showed that 54.3% of highly conserved genes are reconstructed to full length. The commonly seen result of completness for non-model animal ranges from 50% up to 90% depending on the complexity of the species' biology (genome size, repetitive sequences).
```
# BUSCO version is: 5.3.2 
# The lineage dataset is: actinopterygii_odb10 (Creation date: 2021-02-19, number of genomes: 26, number of BUSCOs: 3640)
# Summarized benchmarking in BUSCO notation for file /home/ericjuo/Projects/salmo_trutta_rna_seq/data/03_processed/trinity/trinity.Trinity.fasta
# BUSCO was run in mode: transcriptome

	***** Results: *****

	C:54.3%[S:25.2%,D:29.1%],F:9.3%,M:36.4%,n:3640	   
	1974	Complete BUSCOs (C)			   
	916	Complete and single-copy BUSCOs (S)	   
	1058	Complete and duplicated BUSCOs (D)	   
	337	Fragmented BUSCOs (F)			   
	1329	Missing BUSCOs (M)			   
	3640	Total BUSCO groups searched		   

Dependencies and versions:
	hmmsearch: 3.1
	metaeuk: 5.34c21f2
```
### Estimate abunance of each transcrpt
BAM file is sorted before input to abundance estimation
```
$ samtools sort -@ 10 ata/02_intermediate/bowtie2/aligned_reads.bam > data/02_intermediate/bowtie2/aligned_reads_sorted.bam
```

concatenated pre_trinity paired end reads were used to esitmate abundance
```
$ singularity exec -e docker://registry.docker.com/trinityrnaseq/trinityrnaseq /usr/local/bin/util/align_and_estimate_abundance.pl --transcripts data/03_processed/trinity/trinity.Trinity.fasta --seqType fq --prep_reference --left data/02_intermediate/concatenated_pre_trinity/pre_trinity_1.fastq --right data/02_intermediate/concatenated_pre_trinity/pre_trinity_2.fastq --est_method RSEM --output_dir data/02_intermediate/RSEM_abundance --aln_method bowtie
```
The output is stored in data/02_intermdeiate/RSEM_abundance/RSEM.isoforms.results
```
transcript_id	gene_id	length	effective_length	expected_count	TPM	FPKM	IsoPct
TRINITY_DN0_c0_g1_i11	TRINITY_DN0_c0_g1_i11	1071	894.55	6557.95	390.32	453.17	100.00
TRINITY_DN0_c0_g1_i17	TRINITY_DN0_c0_g1_i17	1756	1579.55	1929.61	65.04	75.52	100.00
TRINITY_DN0_c0_g1_i2	TRINITY_DN0_c0_g1_i2	1727	1550.55	1520.22	52.20	60.61	100.00
TRINITY_DN0_c0_g1_i3	TRINITY_DN0_c0_g1_i3	1533	1356.55	469.59	18.43	21.40	100.00
TRINITY_DN0_c0_g1_i5	TRINITY_DN0_c0_g1_i5	1684	1507.55	539.02	19.04	22.10	100.00
TRINITY_DN0_c0_g1_i7	TRINITY_DN0_c0_g1_i7	1686	1509.55	1346.52	47.49	55.14	100.00
TRINITY_DN0_c0_g1_i9	TRINITY_DN0_c0_g1_i9	1755	1578.55	4198.08	141.60	164.40	100.00
TRINITY_DN10000_c0_g1_i1	TRINITY_DN10000_c0_g1_i1	337	161.66	4.00	1.32	1.53	100.00
TRINITY_DN10000_c1_g1_i1	TRINITY_DN10000_c1_g1_i1	506	329.61	7.00	1.13	1.31	100.00
```
Create transcript and gene expression matrix
```
$ singularity exec -e docker://registry.docker.com/trinityrnaseq/trinityrnaseq /usr/local/bin/util/abundance_estimates_to_matrix.pl --est_method RSEM --gene_trans_map data/03_processed/trinity/trinity.Trinity.fasta.gene_trans_map --out_prefix RSEM --name_sample_by_basedir data/02_intermediate/RSEM_abundance/RSEM.isoforms.results
```
The output directory is the current working directory. 


### Calculate E90N50 using trinity perl script
Generate table group by transcirpt expression level percentil
```
$ singularity exec -e docker://registry.docker.com/trinityrnaseq/trinityrnaseq /usr/local/bin/util/misc/contig_ExN50_statistic.pl RSEM.isoform.counts.matrix data/03_processed/trinity/trinity.Trinity.fasta | tee ExN50.stats
```
Plot transcript expression percentil against N50
```
$ singularity exec -e docker://registry.docker.com/trinityrnaseq/trinityrnaseq /usr/local/bin/util/misc/plot_ExN50_statistic.Rscript  ExN50.stats
```
The output is ExN50_plot.pdf

