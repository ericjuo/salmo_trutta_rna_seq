# Quality trimming using trimmomatic.md
trimmomatic is a quality trimming tool specific for use in fastq data from illumina sequencer.

##  Download trimmomatic
trimmomatic was download from http://www.usadellab.org/cms/?page=trimmomatic.  
Zipped file was unzipped:
```
$ unzip Trimmomatic-0.39.zip
```
Check the content of Trimmomatic-0.39 folder
```
$ tree Trimmomatic-0.39/
Trimmomatic-0.39/
├── LICENSE
├── adapters
│   ├── NexteraPE-PE.fa
│   ├── TruSeq2-PE.fa
│   ├── TruSeq2-SE.fa
│   ├── TruSeq3-PE-2.fa
│   ├── TruSeq3-PE.fa
│   └── TruSeq3-SE.fa
└── trimmomatic-0.39.jar

1 directory, 8 files
```
Remove zipped file
```
rm Trimmomatic-0.39.zip
```

## Quality trimming
fastq files were trimmed by trimmomatic:
```
$ java -jar  Trimmomatic-0.39/trimmomatic-0.39.jar PE SRR799770_1.fastq  SRR799770_2.fastq SRR799770_1_paired_trimmomatic.fastq SRR799770_2_paired_trimmomatic.fastq SRR799770_1_unpaired_trimmomatic.fastq SRR799770_2_unpaired_trimmomatic.fastq ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:35
TrimmomaticPE: Started with arguments:
 SRR799770_1.fastq SRR799770_2.fastq SRR799770_1_paired_trimmomatic.fastq SRR799770_2_paired_trimmomatic.fastq SRR799770_1_unpaired_trimmomatic.fastq SRR799770_2_unpaired_trimmomatic.fastq ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:35
Multiple cores found: Using 4 threads
Using PrefixPair: 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT'
Using Long Clipping Sequence: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
Using Long Clipping Sequence: 'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG'
Using Long Clipping Sequence: 'TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC'
Using Long Clipping Sequence: 'TTTTTTTTTTCAAGCAGAAGACGGCATACGA'
Using Long Clipping Sequence: 'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT'
Using Long Clipping Sequence: 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
ILLUMINACLIP: Using 1 prefix pairs, 6 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Quality encoding detected as phred33
Input Read Pairs: 3198880 Both Surviving: 2750056 (85.97%) Forward Only Surviving: 132911 (4.15%) Reverse Only Surviving: 216323 (6.76%) Dropped: 99590 (3.11%)
TrimmomaticPE: Completed successfully
```
parameters:
- ILLUMINACLIP:2:30:10, clip adapter based on the input adapter sequence file `Trimmomatic-0.39/adapters/TruSeq2-PE.fa`. The allowable mismatch is 2; the score of palindrome score between two ligated adapter is 30; the minimum threshold for matching any adapter is 10.
- SLIDINGWINDOW:4:20, Scan each 4 bases from 5' end and trim from where the average quality score below 20.
- MINLEN: Discard sequences less than 35 bp after trimming.

