# Remove foreign species contamination from reads using Kraken2
Kraken2 is an classification software for classifying species in reads. We can use it to classify human genome and get clean brown trout reads in the unclassified output.

## System requirement
The size of default databasein the kraken is 29G, and it will all load into memory while processing. So, the minimum requirement for RAM size would be about 30G.  

Since I don't have that much RAM space, I used [galaxy platform](https://usegalaxy.org/) to run kraken2 software.

## Steps
### Upload data
Data were uploaded to galaxy platform:
-   SRR799769_1_paired_trimmomatic.fastq.gz


### Kraken parameters
-   --paired
-   --use-name
-   --minimum hit groups 2
-   --unclassified-out and --classified-out
-   --report
-   --db Standard (Created: 2020-11-26T0217062, kmer-len=35, minimizer-len=31, minimizer-space=6) load-factor 0.7
