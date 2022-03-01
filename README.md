# Rerun de novo RNA-Seq assembly
## Initiative
This project was to practice de novo RNA-Seq assembly through rerunning the experiment of brown trout RNA-Seq published by Uren Webster et al. in 2013.

Publication:
Uren Webster, T. M., Bury, N., van Aerle, R., & Santos, E. M. (2013). Global transcriptome profiling reveals molecular mechanisms of metal tolerance in a chronically exposed wild population of brown trout. Environmental science & technology, 47(15), 8869–8877. https://doi.org/10.1021/es401380p  


##  Data source
GSE45637 dataset downloaded from GEO website:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45637


## Analysis outline
-   Remove adapter sequences
-   Filter out sequences that did not meet the quality shreshold
-   Remove seuqences less than 30 bp
-   Pool paried reads of the adult tissues and embryo
-   De novo assemble using Velvet (version 1.2.08) and Oases by kmers
-   Annotate transcripts using Blastn and Blastx and fish and mammalian nt and protein databases with a cut off e-value < 1E-15.
-   Analyse gene expression using RSEM with `--no_polyA` parameter and default settings.
-   Statistical analysis differential expression using edgeR with a cutoff FDR <0.1.
-   Draw differential expressed genes using VennDiagram in R, Bioconductor.

