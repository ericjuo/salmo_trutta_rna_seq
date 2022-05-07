
Download brown trout genome from ensembl
```
rsync -avh --exclude '*.toplevel.fa.gz' --exclude '*dna_rm.*' --exclude '*dna_sm.*' --exclude 'README' --exclude 'CHECKSUMS' rsync://ftp.ensembl.org/ensembl/pub/release-106/fasta/salmo_trutta/dna/ .
skipping directory .
```
Gunzipped files
```
gunzip contaminants/brown_trout_genome/*
```

Change directory 
```
cd contaminants/brown_trout_genome/ 
```


Build indexed genome with HISAT2
```
hisat2-build Salmo_trutta.fSalTru1.1.dna.nonchromosomal.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.10.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.11.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.12.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.13.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.14.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.15.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.16.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.17.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.18.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.19.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.1.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.20.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.21.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.22.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.23.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.24.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.25.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.26.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.27.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.28.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.29.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.2.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.30.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.31.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.32.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.33.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.34.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.35.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.36.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.37.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.38.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.39.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.3.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.40.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.4.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.5.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.6.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.7.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.8.fa,Salmo_trutta.fSalTru1.1.dna.primary_assembly.9.fa brown_trout
```

Change to top directory
```
cd ../..
```

Map reads to genome
```
mkdir data/02_intermediate/hisat2
mkdir logs/hisat2

hisat2 -p 10 -x contaminants/brown_trout_genome/brown_trout -1 data/02_intermediate/pre_trinity_formatting/SRR799769_pre_trinity_1.fastq -2 data/02_intermediate/pre_trinity_formatting/SRR799769_pre_trinity_2.fastq -S data/02_intermediate/hisat2/SRR799769_map_brown_trout_genome.sam --summary-file logs/hisat2/SRR799769_map_brown_trout_genome.log

hisat2 -p 10 -x contaminants/brown_trout_genome/brown_trout -1 data/02_intermediate/pre_trinity_formatting/SRR799770_pre_trinity_1.fastq -2 data/02_intermediate/pre_trinity_formatting/SRR799770_pre_trinity_2.fastq -S data/02_intermediate/hisat2/SRR799770_map_brown_trout_genome.sam --summary-file logs/hisat2/SRR799770_map_brown_trout_genome.log

hisat2 -p 10 -x contaminants/brown_trout_genome/brown_trout -1 data/02_intermediate/pre_trinity_formatting/SRR799771_pre_trinity_1.fastq -2 data/02_intermediate/pre_trinity_formatting/SRR799771_pre_trinity_2.fastq -S data/02_intermediate/hisat2/SRR799771_map_brown_trout_genome.sam --summary-file logs/hisat2/SRR799771_map_brown_trout_genome.log

hisat2 -p 10 -x contaminants/brown_trout_genome/brown_trout -1 data/02_intermediate/pre_trinity_formatting/SRR799772_pre_trinity_1.fastq -2 data/02_intermediate/pre_trinity_formatting/SRR799772_pre_trinity_2.fastq -S data/02_intermediate/hisat2/SRR799772_map_brown_trout_genome.sam --summary-file logs/hisat2/SRR799772_map_brown_trout_genome.log

hisat2 -p 10 -x contaminants/brown_trout_genome/brown_trout -1 data/02_intermediate/pre_trinity_formatting/SRR799773_pre_trinity_1.fastq -2 data/02_intermediate/pre_trinity_formatting/SRR799773_pre_trinity_2.fastq -S data/02_intermediate/hisat2/SRR799773_map_brown_trout_genome.sam --summary-file logs/hisat2/SRR799773_map_brown_trout_genome.log

hisat2 -p 10 -x contaminants/brown_trout_genome/brown_trout -1 data/02_intermediate/pre_trinity_formatting/SRR799774_pre_trinity_1.fastq -2 data/02_intermediate/pre_trinity_formatting/SRR799774_pre_trinity_2.fastq -S data/02_intermediate/hisat2/SRR799774_map_brown_trout_genome.sam --summary-file logs/hisat2/SRR799774_map_brown_trout_genome.log

hisat2 -p 10 -x contaminants/brown_trout_genome/brown_trout -1 data/02_intermediate/pre_trinity_formatting/SRR799775_pre_trinity_1.fastq -2 data/02_intermediate/pre_trinity_formatting/SRR799775_pre_trinity_2.fastq -S data/02_intermediate/hisat2/SRR799775_map_brown_trout_genome.sam --summary-file logs/hisat2/SRR799775_map_brown_trout_genome.log

hisat2 -p 10 -x contaminants/brown_trout_genome/brown_trout -1 data/02_intermediate/pre_trinity_formatting/SRR799776_pre_trinity_1.fastq -2 data/02_intermediate/pre_trinity_formatting/SRR799776_pre_trinity_2.fastq -S data/02_intermediate/hisat2/SRR799776_map_brown_trout_genome.sam --summary-file logs/hisat2/SRR799776_map_brown_trout_genome.log
```

Download annotation, general feature format (GFF3), file from esembl
http://ftp.ensembl.org/pub/release-106/gff3/salmo_trutta/

```
wget http://ftp.ensembl.org/pub/release-106/gff3/salmo_trutta/Salmo_trutta.fSalTru1.1.106.gff3.gz
```