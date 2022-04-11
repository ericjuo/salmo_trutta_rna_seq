# De novo assembling using Trinity
After clean up foreign species contaminants, forward reads and reverse reads were separately pooled together. Pooled forward reads and pooled reverse reads were uploaded to galaxy for de nova assembling.

## Steps to pool reads
`cat` command was used to concatenate gzipped forward and reverse reads:
```
$ cat cat $(ls *_Unclassified_1.fastq.gz) > Pooled_reads_1.fastq.gz
$ cat cat $(ls *_Unclassified_2.fastq.gz) > Pooled_reads_2.fastq.gz
```

## Parameters to run Trinity
-- left Pooled_reads_1.fastq.gz  
-- right Pooled_reads_2.fastq.gz  
-- min_contig_length 200  
-- min_kmer_cov 1  

The default kmer length is 25bp, and this can't be changed in the Trinity in the galaxy platform.  
The minimum kmer coverage was setted to 1, which will extend contigs if any k-1 mer is matched.
