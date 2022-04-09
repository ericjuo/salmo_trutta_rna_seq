#!/bin/sh

prefix=/home/ericjuo/Projects/salmo_trutta_rna_seq/.snakemake/conda/9ad3273d11a628166183e9b102ca99a6
exec_prefix=/home/ericjuo/Projects/salmo_trutta_rna_seq/.snakemake/conda/9ad3273d11a628166183e9b102ca99a6
libdir=${exec_prefix}/lib

LD_PRELOAD=${libdir}/libjemalloc.so.2
export LD_PRELOAD
exec "$@"
