# Intermediate Data
This intermdiate directory stores data processed by softwares and only for intermiate use.

## Content
```
$ ls -1
README.md
SRR799769_1.fastq
SRR799769_2.fastq
SRR799770_1.fastq
SRR799770_2.fastq
SRR799771_1.fastq
SRR799771_2.fastq
SRR799772_1.fastq
SRR799772_2.fastq
SRR799773_1.fastq
SRR799773_2.fastq
SRR799774_1.fastq
SRR799774_2.fastq
SRR799775_1.fastq
SRR799775_2.fastq
SRR799776_1.fastq
SRR799776_2.fastq
```
The SRR799769_1-SRR799776_2.fastaq files were created by `fasterq-dump` with input SRA files in the 01_raw folder.
```
$ for i in data/01_raw/*/*.sra; do fasterq-dump $i -O data/02_intermediate/; done
```

## FastQC analysis
Raw reads were inspected for their quality by FastQC application.
```
$ for i in *.fastq; do fastqc $i; done
Started analysis of SRR799769_1.fastq
Approx 5% complete for SRR799769_1.fastq
Approx 10% complete for SRR799769_1.fastq
Approx 15% complete for SRR799769_1.fastq
Approx 20% complete for SRR799769_1.fastq
Approx 25% complete for SRR799769_1.fastq
Approx 30% complete for SRR799769_1.fastq
Approx 35% complete for SRR799769_1.fastq
Approx 40% complete for SRR799769_1.fastq
Approx 45% complete for SRR799769_1.fastq
Approx 50% complete for SRR799769_1.fastq
Approx 55% complete for SRR799769_1.fastq
Approx 60% complete for SRR799769_1.fastq
Approx 65% complete for SRR799769_1.fastq
Approx 70% complete for SRR799769_1.fastq
Approx 75% complete for SRR799769_1.fastq
Approx 80% complete for SRR799769_1.fastq
Approx 85% complete for SRR799769_1.fastq
Approx 90% complete for SRR799769_1.fastq
Approx 95% complete for SRR799769_1.fastq
Analysis complete for SRR799769_1.fastq
...
```
