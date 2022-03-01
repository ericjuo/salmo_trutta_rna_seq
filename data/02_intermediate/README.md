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
