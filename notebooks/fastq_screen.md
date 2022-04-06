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
In the configuration file, I edited the path of bowtie2 executable and databases (Human, Mouse, and E coli).

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

The E coli genome was downloaded from EMBL, and indexed using bowtie2.
E coli genome: https://www.ebi.ac.uk/ena/browser/api/fasta/U00096.2?download=true
bowtie2 parameters:
```
bowtie2-build U00096.2.fasta ecoli
```

## fastq_screen analysis
Trimmomatic_trimmed files were subjected to fastq_screen using ptyhon code:
```
ptyhon ../../scripts/run_fastq_screen.py
```

The foreign species contaminants depleted files from Kraken2 were analysed with FastQ screen again:
```
for i in $(ls *_Unclassified_[1,2].fastq);do ~/bin/FastQ-Screen-0.15.2/fastq_screen --aligner bowtie2 $i; done
```