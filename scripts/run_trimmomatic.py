import os
import sys
import subprocess

Sample = ['SRR7997' + str(x) for x in range(69, 77)]

for i in Sample:
    subprocess.run(f"java -jar  ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -trimlog {i}_log.txt {i}_1.fastq  {i}_2.fastq {i}_1_paired_trimmomatic.fastq {i}_1_unpaired_trimmomatic.fastq {i}_2_paired_trimmomatic.fastq {i}_2_unpaired_trimmomatic.fastq ILLUMINACLIP:/home/ericjuo/bin/Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:35", shell=True)