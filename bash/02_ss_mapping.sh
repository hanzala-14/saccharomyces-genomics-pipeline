#!/usr/bin/env bash

for inname in ./*R1.filt.f*q.gz; do
ifname2=${inname/R1.filt.f*q.gz/R2.f*q.gz}
  sam1name=${inname/R1.filt.f*q.gz/sam}
    samname=$(basename "$sam1name")
     REFERENCE=/home/user1/Desktop/Hanzala/References/sensu_stricto_formatted.fasta
     THREADS=6

bwa mem -t $THREADS $REFERENCE $inname $ifname2 > $samname &&

# 1. Convert SAM to BAM

  ofname1=${samname/.sam/.tempfile.raw.bam}

samtools view -Sb -F4 -@6 $samname > $ofname1 &&

# 2. sort BAM

  ofname2=${samname/.sam/.tempfile.sort.bam}

samtools sort -@$THREADS $ofname1 > $ofname2 &&

# 3. mapping stats 

  ofname3=${samname/.sam/.mstat}

samtools flagstat $samname > $ofname3 &&

# 4. Rename readgroups for GATK

  ofname4=${samname/.sam/.sort.tempfile.rg.bam}

java -jar /home/user1/picard.jar AddOrReplaceReadGroups \
 RGID=sthg RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$samname \
 I=$ofname2 \
 O=$ofname4 &&

# 5. Mark duplicates for GATK

  ofname5=${samname/.sam/.mdup.bam}
  tmpname1=${samname/.sam/.metrics}

java -Xmx18G -jar /home/user1/picard.jar MarkDuplicates I=$ofname4 O=$ofname5 M=$tmpname1 &&

# 6. Index bam for GATK

samtools index $ofname5 &&

rm ./*tempfile* ;
rm $samname

done
