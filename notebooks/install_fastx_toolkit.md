# Install fastx_toolkit
fastx_toolkit was installed following the instruction from the [website](http://hannonlab.cshl.edu/fastx_toolkit/install_ubuntu.txt) with serveral modifcations.

## Install libgtextutils (dependency of fastx-toolkit)
```
$ wget --timeout=10 http://hannonlab.cshl.edu/fastx_toolkit/libgtextutils-0.6.1.tar.bz2
```
Check integrity of downloaded file with md5sum
```
$ md5sum -c <<<"40e7df4e5a72efe50aa789af8caeb935 libgtextutils-0.6.1.tar.bz2"
libgtextutils-0.6.1.tar.bz2: OK
```
Compile libgtextutils.  
I ran into error while complining with make. I checked stackoverflow and found an [answer](https://stackoverflow.com/questions/38659115/make-fails-with-error-cannot-convert-stdistream-aka-stdbasic-istreamchar). The error is stemmed from newer version of GCC defaulting to C++14 mode instead of C++03 mode. The solution is override the C++ flags to `-std=c++03 -01`
```
$ tar -xjf libgtextutils-0.6.tar.bz2
$ cd libgtextutils-0.6.1/
$ ./configure
$ make CXXFLAGS='-std=c++03 -O1'
$ sudo make install
```
##  Install fastx_toolkit
```
$ wget --timeout=10 http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit-0.0.13.2.tar.bz2
```
Check integrity of downloaded file with md5sum
```
$ md5sum -c <<<"1d1238cb3029cb1c5d1a3eea7d7d2ca8 fastx_toolkit-0.0.13.2.tar.bz2"
fastx_toolkit-0.0.13.2.tar.bz2: OK
```
Compile fastx_toolkit
```
$ tar -xjf fastx_toolkit-0.0.13.2.tar.bz2
$ cd fastx_toolkit-0.0.13.2/
$ ./configure
$ make CXXFLAGS=-O1
$ sudo make install
```

##  add /usr/local/lib to LD_LIBRARY_PATH environment variable
```
$ export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
```

## Sanity check
```
$ fastx_uncollapser -h
usage: fasta_uncollapser [-c N] [-h] [-v] [-i INFILE] [-o OUTFILE]
Part of FASTX Toolkit 0.0.13.2 by A. Gordon (gordon@cshl.edu)

   [-h]         = This helpful help screen.
   [-v]         = verbose: print short summary of input/output counts
   [-c N]       = Assume input is a tabular file (not FASTA file),
                  And the collapsed identifier (e.g. '1-1000') is on column N.
   [-i INFILE]  = FASTA/Tabular input file. default is STDIN.
   [-o OUTFILE] = FASTA/Tabular output file. default is STDOUT.
```
