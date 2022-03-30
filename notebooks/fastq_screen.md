# Fastq Screen
Fastq Screen is a program to find potential contamination in reads.

## Installation of FastQ screen.
Fastq screen software was downloaded from the [homepage](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/). My install version is 0.15.2.

Unpack the downloaded file:
```
tar -xzf fastq_screen_v0.15..tar.gz
```

Install dependency, bowtie2
```
sudo apt update
sudo apt install bowtie2
```

fastq_screen requires configuration before use, the configuration is in the following path:
```
/home/ericjuo/bin/FastQ-Screen-0.15.2/fastq_screen.conf
```
In the configuration file, I edited the path of bowtie2 executable and databases (Human, Mouse).
The Human genome was downloaded from ensembl, and indexed using bowite2.
Human genome: http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz 
bowtie2 parameters:
```
bowtie2-build Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz Homo_sapiens.GRCh38
```

The indexed mouse genome was downloaded directly from bowtie2 website:
Mouse genome: https://genome-idx.s3.amazonaws.com/bt/GRCm39.zip
Genome was downloaded using wget:
```
wget https://genome-idx.s3.amazonaws.com/bt/GRCm39.zip
```
