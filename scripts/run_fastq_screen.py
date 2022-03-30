import subprocess
import os

for file in os.listdir("./"):
    print(file)
    if file.endswith("paired_trimmomatic.fastq.gz"):
        print(file)
        subprocess.run(f"~/bin/FastQ-Screen-0.15.2/fastq_screen --aligner bowtie2 {file}", shell=True)