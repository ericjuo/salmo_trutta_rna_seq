# Raw Data
This raw data direcotry mainly stores sequencing data downloaded directly from the internet.

## Content
--- 
```
$ ls -1
README.md
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
- SRR799769-SRR799776 are folders containing SRA files of RNA-Seq from embryos, gill, gut, kidney and liver of brown trout exposed to metals or in control river.
-   SRR_Acc_List.txt was downloaded from GEO entry of Uren Webster's publication using SRA toolkit online. (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45637)

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

## Download raw reads from Gene Expression Omnibus (GEO), NCBI
--- 
 
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
    
##  Convert SRA file into FastQ file
The SRA files were converted into FastQ files using `fasterq-dump` function of SRA toolkits. The converted files were stored in `data/02_intermediate` folder.
```
$ for i in data/01_raw/*/*.sra; do fasterq-dump $i -O data/02_intermediate/; done
spots read      : 2,936,116
reads read      : 5,872,232
reads written   : 5,872,232
spots read      : 3,198,880
reads read      : 6,397,760
reads written   : 6,397,760
spots read      : 3,587,826
reads read      : 7,175,652
reads written   : 7,175,652
spots read      : 3,648,593
reads read      : 7,297,186
reads written   : 7,297,186
spots read      : 3,489,607
reads read      : 6,979,214
reads written   : 6,979,214
spots read      : 3,418,757
reads read      : 6,837,514
reads written   : 6,837,514
spots read      : 3,147,358
reads read      : 6,294,716
reads written   : 6,294,716
spots read      : 3,855,726
reads read      : 7,711,452
reads written   : 7,711,452
```


