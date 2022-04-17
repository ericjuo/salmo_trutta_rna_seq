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
$ bowtie2 -p 10 -q --no-unal -k 20 -x brown_trout_trinity.fai -1 pre_trinity_1.fasta -2 pre_trinity_2.fasta 2> align-stats.txt | samtools view -@ 10 -b brown_trout.bam
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

