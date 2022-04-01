# Rerun de novo RNA-Seq assembly
## Initiative
This project is to reproduce the results of brown trout RNA-Seq analysis published by Uren Webster et al. in 2013.

Publication:
Uren Webster, T. M., Bury, N., van Aerle, R., & Santos, E. M. (2013). Global transcriptome profiling reveals molecular mechanisms of metal tolerance in a chronically exposed wild population of brown trout. Environmental science & technology, 47(15), 8869–8877. https://doi.org/10.1021/es401380p  


##  Data source
GSE45637 dataset downloaded from GEO website:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45637

Steps to download raw reads: [[Link](notebooks/obtain_raw_reads_from_GEO.md)].

## Analysis outline
- [x] Remove adapter sequences (Base on FASTQC reports, no adapter contamination was found)  
- [x] Filter out sequences that did not meet the quality threshold  
- [x] Remove seuqences less than 35 bp  
- [ ] Pool paried reads of the adult tissues and embryo  
- [ ] De novo assemble using Velvet (version 1.2.08) and Oases by kmers  
- [ ] Annotate transcripts using Blastn and Blastx and fish and mammalian nt and protein databases with a cut off e-value < 1E-15.  
- [ ] Analyse gene expression using RSEM with `--no_polyA` parameter and default settings.  
- [ ] Statistical analysis differential expression using edgeR with a cutoff FDR <0.1.  
- [ ] Draw differential expressed genes using VennDiagram in R, Bioconductor.  

##  Analysis
-   The matadata analysis show that there are 2 sample groups (fish group from control River Teign and high metal content River Hayle ). For each group, fish's gill, stomach and intestine, trunk kidney and liver were taken to prepare RNA-Seq.  
Steps to perform metadata analysis: [[Link](notebooks/metadata_analysis.ipynb)]  

    |    | Run       |       Gb | tissue                | treatment              |
    |---:|:----------|---------:|:----------------------|:-----------------------|
    |  0 | SRR799769 | 0.255852 | gill                  | chronic metal exposure |
    |  1 | SRR799770 | 0.280547 | stomach and intestine | chronic metal exposure |
    |  2 | SRR799771 | 0.312282 | trunk kidney          | chronic metal exposure |
    |  3 | SRR799772 | 0.321309 | liver                 | chronic metal exposure |
    |  4 | SRR799773 | 0.305214 | gill                  | control                |
    |  5 | SRR799774 | 0.299001 | stomach and intestine | control                |
    |  6 | SRR799775 | 0.272784 | trunk kidney          | control                |
    |  7 | SRR799776 | 0.338982 | liver                 | control                |  

-   The initial qualtiy of these raw reads were inspected using FastQC software. 

    [SRR799769_1](https://htmlpreview.github.io/?https://github.com/ericjuo/salmo_trutta_rna_seq/blob/master/data/02_intermediate/SRR799769_1_fastqc.html),
    [SRR799769_2](https://htmlpreview.github.io/?https://github.com/ericjuo/salmo_trutta_rna_seq/blob/master/data/02_intermediate/SRR799769_2_fastqc.html),
    [SRR799770_1](https://htmlpreview.github.io/?https://github.com/ericjuo/salmo_trutta_rna_seq/blob/master/data/02_intermediate/SRR799770_1_fastqc.html),
    [SRR799770_2](https://htmlpreview.github.io/?https://github.com/ericjuo/salmo_trutta_rna_seq/blob/master/data/02_intermediate/SRR799770_2_fastqc.html),
    [SRR799771_1](https://htmlpreview.github.io/?https://github.com/ericjuo/salmo_trutta_rna_seq/blob/master/data/02_intermediate/SRR799771_1_fastqc.html),
    [SRR799771_2](https://htmlpreview.github.io/?https://github.com/ericjuo/salmo_trutta_rna_seq/blob/master/data/02_intermediate/SRR799771_2_fastqc.html),
    [SRR799772_1](https://htmlpreview.github.io/?https://github.com/ericjuo/salmo_trutta_rna_seq/blob/master/data/02_intermediate/SRR799772_1_fastqc.html),
    [SRR799772_2](https://htmlpreview.github.io/?https://github.com/ericjuo/salmo_trutta_rna_seq/blob/master/data/02_intermediate/SRR799772_2_fastqc.html),
    [SRR799773_1](https://htmlpreview.github.io/?https://github.com/ericjuo/salmo_trutta_rna_seq/blob/master/data/02_intermediate/SRR799773_1_fastqc.html),
    [SRR799773_2](https://htmlpreview.github.io/?https://github.com/ericjuo/salmo_trutta_rna_seq/blob/master/data/02_intermediate/SRR799773_2_fastqc.html),
    [SRR799774_1](https://htmlpreview.github.io/?https://github.com/ericjuo/salmo_trutta_rna_seq/blob/master/data/02_intermediate/SRR799774_1_fastqc.html),
    [SRR799774_2](https://htmlpreview.github.io/?https://github.com/ericjuo/salmo_trutta_rna_seq/blob/master/data/02_intermediate/SRR799774_2_fastqc.html),
    [SRR799775_1](https://htmlpreview.github.io/?https://github.com/ericjuo/salmo_trutta_rna_seq/blob/master/data/02_intermediate/SRR799775_1_fastqc.html),
    [SRR799775_2](https://htmlpreview.github.io/?https://github.com/ericjuo/salmo_trutta_rna_seq/blob/master/data/02_intermediate/SRR799775_2_fastqc.html),
    [SRR799776_1](https://htmlpreview.github.io/?https://github.com/ericjuo/salmo_trutta_rna_seq/blob/master/data/02_intermediate/SRR799776_1_fastqc.html),
    [SRR799776_2](https://htmlpreview.github.io/?https://github.com/ericjuo/salmo_trutta_rna_seq/blob/master/data/02_intermediate/SRR799776_2_fastqc.html) 

- Shift in GC content peak indicates potential contamination with other speices.  
![gc](./data/02_intermediate/SRR799769_1_fastqc_gc.png)  

-   Biased nucleotide content at 5'end of reads is commonly seen in RNA-Seq dataset, which is stemmed from priming on fragmented RNA with hexamer or fragmentation with transposases ([Hansen et al., 2010](https://academic.oup.com/nar/article/38/12/e131/2409775)). But it doesn't affect the downstream de novo asseblmly.  
![bnc](./data/02_intermediate/SRR799769_1_fastqc.png)

- Raw reads were subjected to quality trimming using trimmomatic. The leading and trailing 3 bases of reads were trimmed. Those bases with phred scores below 20 were trimmed, and those read length shorter than 35 bp after trimming were discarded. In addition, the remaining illumina adapter were also removed.    
Quality report of trimmomatic-trimmed reads from FastQC software: [[Link](https://htmlpreview.github.io/?https://github.com/ericjuo/salmo_trutta_rna_seq/blob/master/data/02_intermediate/SRR799770_1_paired_trimmomatic_fastqc.html)]

- Raw reads were also trimmed with sickle software for comparison.   
Quality report of sickle-trimmed reads from FastQC software: [[Link](https://htmlpreview.github.io/?https://github.com/ericjuo/salmo_trutta_rna_seq/blob/master/data/02_intermediate/SRR799770_1_sickle_fastqc.html)]

- There is no much difference between trimming with trimmomatic or sickle program. Both trimming program sucessfully trimmed low quality score bases in raw reads.  
![BvA](./data/02_intermediate/SRR799770_1_Before_After_trimmed.jpg)  

- MultiQC report summarized the quality report of before and after trimming either using trimmomatic or sickle: [[Link](https://htmlpreview.github.io/?https://github.com/ericjuo/salmo_trutta_rna_seq/blob/master/data/02_intermediate/multiqc_report_1.html)]

- Trommatic-trimmed reads were searched against human, mouse and E. Coli genome using Fastq Screen to identify foreign species contamination.
![Contamination](./data/02_intermediate/SRR799769_1_paired_trimmomatic_screen.png)

- Foreign species contamination can be removed by Kraken2 or Centrifuge software. However, both softwares require at least 30G RAM, which totally beyond my desktop's RAM size (I only have 16G RAM) ([CCB website](http://ccb.jhu.edu/software/choosing-a-metagenomics-classifier/)). Despite I  really want to remove foreign species contamination de novo assebmly, I can't perform it at this moment.
