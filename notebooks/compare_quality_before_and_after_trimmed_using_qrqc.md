# Compare quality before and after trimmed using qrqc module

The qualities of before and after-trimmed files were compared using bioconductor module, qrqc.

##  Install bioconductor in R
In R console:
```
> if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")
```
Check Bioconductor version:
```
> BiocManager::version()
[1] '3.14'
```
## Install qrqc module
Install qrqc module through BiocManager
```
> BiocManager::install("qrqc")
```

## Load in qrqc module
```
> library(qrqc)
```

## Produce quality plot
```
> fqfiles <- c(none="untreated1_chr4.fq", sickle="untreated1_chr4_sickle.fq", trimfq="untreated1_chr4_trimfq.fq")
> seq_info <- lapply(fqfiles, function(file){readSeqFile(file, hash = FALSE, kmer = FALSE)})
> p <- qualPlot(seq_info, quartile.color = NULL, mean.color = NULL) + theme_bw() + scale_y_continuous("qualty (Phred33)")
> jpeg(file="SRR799770_1_Before_After_trimmed.jpg")
> print(p)
> dev.off()
```
![BvA](../data/02_intermediate/SRR799770_1_Before_After_trimmed.jpg)