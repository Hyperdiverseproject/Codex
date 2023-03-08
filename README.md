# Codexpipeline takes fastq files of paired-end reads and a tag file of tags of each fragments
Usage: $0 [-r|--R1 <R1 fastq file, required>] [-f|--R2 <R2 fastq file, required>] [-t|--tag <Barcode file, tabular, required>] 
[-p|-prj <project name directory>]

-h, --help              This help
-f, --R1 <file.fastq>   Fastq file unzipped with R1 reads; required
-r, --R2 <file.fastq>   Fastq file unzipped with R2 reads; required
-t, --tag <file.txt>    Barcode file; one line per well with the well name, the fragment 1 forward tag, then the fragment 1 reverse tag then the fragment 2 forward tag, then the fragm
ent 2 reverse tag separated by a tabulation; required
-p, --prj <directory>   Project name; Creates a directory were computing will be done; default BTC1
-m, --min <number>      Minimum number of reads for a fragment; default 10
-a, --minfr1 <number>   Minimum length for fragment 1; default 331
-b, --maxfr1 <number>   Maximum length for fragment 1; default 451
-c, --minfr2 <number>   Minimum length for fragment 2; default 313
-d, --maxfr2 <number>   Maximum length for fragment 2; default 434
-n, --minreads <number> Minimum reads matching a contig fragment; default 10
-e, --mincov <number>   Minimum contig fragment converage; default 5

First step - meerging R1 and R2 - need program PEAR under licence: https://www.h-its.org/downloads/pear-academic/
        Install Pear locally:
        Once downloaded, unzip the file pear-0.9.11-linux-x86_64.tar.gz by running:
                tar xzf pear-0.9.11-linux-x86_64.tar.gz
        then export it to your PATH:
                export PATH=\$PATH:\$PWD/pear-0.9.11-linux-x86_64/bin/

Second step - demultiplexing per well - no external program

Third step - demultiplexing per fragment type - need Blast - remove reads lower than 200 nucleotides

Four step - First step assembly per fragment - need trinity 

Fifth step - Second step assembly (assembly of both fragments) - need Blast - Contig filtering (see options)
