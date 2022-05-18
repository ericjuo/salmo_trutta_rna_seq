# Count reads to exons on brown trout genome
## Build reference genome index
```
$ PREFIX=contaminants/brown_trout_genome
$ ls -d --format=commas $PREFIX/Salmo_trutta.fSalTru1.1.dna.primary_assembly*.fa | tr -d "\n" > ref
$ hisat2-build $(cat ref) $PREFIX/brown_trout
```
## Perform gene level read couting
```
$ mkdir -p data/02_intermediate/bam
$ parallel -j 1 echo SRR7997{1}_pre_trinity ::: 69 70 71 71 72 73 74 75 76 > ids
$ cat ids | parallel -j 1 "hisat2 -x $PREFIX/brown_trout -1 data/02_intermediate/pre_trinity_formatting/{}_1.fastq -2 data/02_intermediate/pre_trinity_formatting/{}_2.fastq | samtools sort > data/02_intermediate/bam/{}.bam"
```

## Perform featureCounts to obtain quantification matrix
```
$ mkdir -p report/featureCounts
$ wget http://ftp.ensembl.org/pub/release-106/gff3/salmo_trutta/Salmo_trutta.fSalTru1.1.106.chr.gff3.gz
$ featureCounts -T 8 -p -F GFF3 -a contaminants/brown_trout_genome/Salmo_trutta.fSalTru1.1.106.chr.gff3.gz -t exon -g Name -o report/featureCounts/counts_against_genome.txt data/02_intermediate/bam/*.bam
```

## Calculate differential expression using edgeR
```
counts = read.delim("report/featureCounts/counts_against_genome.txt", header=TRUE)
