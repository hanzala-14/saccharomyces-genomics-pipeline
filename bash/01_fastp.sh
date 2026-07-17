#!/usr/bin/env bash

for infile in *R1.fastq.gz; do
   infile2=${infile/R1.fastq.gz/R2.fastq.gz}
   outfile1=${infile/R1.fastq.gz/.filt.R1.fastq.gz}
   outfile2=${infile/R1.fastq.gz/.filt.R2.fastq.gz}
   htmlrep=${infile/R1.fastq.gz/.report.html}
   jsonrep=${infile/R1.fastq.gz/.report.json}

fastp \
   -i $infile \
   -I $infile2 \
   -o $outfile1 \
   -O $outfile2 \
   --dont_overwrite \
   -G \
   -l 30 \
   -h $htmlrep \
   -j $jsonrep

done
