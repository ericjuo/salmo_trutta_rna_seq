# Rerun de novo RNA-Seq assembly
## Initiative
This project was to practice de novo RNA-Seq assembly through rerunning the experiment of brown trout RNA-Seq published by Uren Webster et al. in 2013.

Publication:
Uren Webster, T. M., Bury, N., van Aerle, R., & Santos, E. M. (2013). Global transcriptome profiling reveals molecular mechanisms of metal tolerance in a chronically exposed wild population of brown trout. Environmental science & technology, 47(15), 8869–8877. https://doi.org/10.1021/es401380p  


## Download raw reads from Gene Expression Omnibus (GEO), NCBI
 
-   The raw reads associated with Uren Webste's publication were downloaded from Series GSE45637 on the GEO website on Feb 28, 2022:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45637

-   Steps to download dataset:
    1. Get list of accession numbers through `SRA RUN Selector`:  
    List of accession numbers was downloaded in a file named `SRR_Acc_List.txt`.
    2. SRA files was fetched using `prefetch` function of SRA toolkit:
    ```
    prefetch --option-file SRR_Acc_List.txt
    ```
    3. Downloaded SRA files were stored separately in folders in its name:
    ```
    ls -1 data/01_raw/
    SRR799769
    SRR799770
    SRR799771
    SRR799772
    SRR799773
    SRR799774
    SRR799775
    SRR799776
    SRR_Acc_List.txt
    ```
    



### Download SRA Toolkit
The SRA Toolkit was downloaded following the instruction describe in https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit  

**For installing on Ubuntu**  
-   Fetch the tar file from the canonical lacation at NCBI:
```
wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
```

-   Extract the contents of the tar file:
```
tar -vxzf sratoolkit.tar.gz
```

-   Append the path to the `PATH` environment variables:
```
export PATH=$PATH:$PWD/sratoolkit.2.11.2-ubuntu64/bin
```

-   Verify the binaries will be found by the shell:
```
which fastq-dump
```
The output was:
```
/home/ericjuo/sratoolkit.2.11.2-ubuntu64/bin/fastq-dump
```

