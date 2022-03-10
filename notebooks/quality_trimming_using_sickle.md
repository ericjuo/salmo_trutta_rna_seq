# Quality trimming using sickle
Sickle is a quality trimming tool for fastq files.

## Download sickle
sickle was downloaded on Ubuntu:
-   Update package list in apt
```
$ sudo apt update
```
-   Check if sickle is in the package list
```
$ apt search sickle
Full Text Search... Done 
sickle/focal 1.33+git20150314.f3d6ae3-1 amd64 
windowed adaptive trimming tool for FASTQ files using quality
```
-   Install sickle
```
$ sudo apt install sickle
```
-   Sanity check
```
$ sickle
Usage: sickle <command> [options] 
Command: 
pe      paired-end sequence trimming 
se      single-end sequence trimming 
--help, display this help and exit 
--version, output version information and exit
```

##  Quality trimming using sickle
fastq files were quality trimmed using sickle
```
$ sickle pe -f SRR799770_1.fastq -r SRR799770_2.fa
stq -t sanger -o SRR799770_1_sickle.fastq -p SRR799770_2_sickle.fastq -s SRR799770_singles_sickle.fastq -q 20 -l 35 -x

PE forward file: SRR799770_1.fastq
PE reverse file: SRR799770_2.fastq

Total input FastQ records: 6397760 (3198880 pairs)

FastQ paired records kept: 5709724 (2854862 pairs)
FastQ single records kept: 278970 (from PE1: 111575, from PE2: 167395)
FastQ paired records discarded: 130096 (65048 pairs)
FastQ single records discarded: 278970 (from PE1: 167395, from PE2: 111575)
```

