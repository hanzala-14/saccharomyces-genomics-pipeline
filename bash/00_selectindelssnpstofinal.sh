#!/bin/bash

REFERENCE=/home/racz_hanna/nemeth_balint/02_bioinformatics/references/sensu_stricto_formatted.fasta  # BWA-indexed reference genome

####################################################################################################################################

for infile in *called.vcf.gz; do

   outfile1=${infile/.vcf.gz/.snp.vcf.gz}

gatk SelectVariants \
 -V $infile \
 -select-type SNP \
 -XL /home/racz_hanna/nemeth_balint/02_bioinformatics/genomics/01_scripts/phylogeny/cerevisiae_repetitive.intervals \
 -O $outfile1 &&

#filter indel

outfile2=${infile/.vcf.gz/.indel.vcf.gz}

gatk SelectVariants \
 -V $infile \
 -select-type INDEL \
 -XL /home/racz_hanna/nemeth_balint/02_bioinformatics/genomics/01_scripts/phylogeny/cerevisiae_repetitive.intervals \
 -O $outfile2
 
    outfile3=${infile/.vcf.gz/.indelvarfiltered.vcf.gz}

gatk VariantFiltration \
 -R $REFERENCE \
 -V $outfile2 \
 -filter "QD < 5.0" --filter-name "QD5" \
 -filter "QUAL < 30.0" --filter-name "QUAL30" \
 -filter "FS > 60.0" --filter-name "FS60" \
 -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
 -O $outfile3
 
    outfile4=${infile/.vcf.gz/indelvarfilteredleftal.vcf.gz}

gatk LeftAlignAndTrimVariants \
 -R $REFERENCE \
 -V $outfile3 \
 --dont-trim-alleles \
 -O $outfile4
 
    outfile5=${infile/.vcf.gz/.snpvarfiltered.vcf.gz}

gatk VariantFiltration \
 -R $REFERENCE \
 -V $outfile1 \
 -filter "QD < 5.0" --filter-name "QD2" \
 -filter "QUAL < 30.0" --filter-name "QUAL30" \
 -filter "SOR > 3.0" --filter-name "SOR3" \
 -filter "FS > 60.0" --filter-name "FS60" \
 -filter "MQ < 40.0" --filter-name "MQ40" \
 -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
 -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
 -O $outfile5
 
   outfile6=${infile/.vcf.gz/.final.vcf.gz}

java -jar /home/racz_hanna/nemeth_balint/02_bioinformatics/softwares/Picard/picard.jar MergeVcfs \
 I=$outfile4 \
 I=$outfile5 \
 O=$outfile6
 
    outfile7=${infile/.vcf.gz/.final.onlyvariant.vcf.gz}

gatk SelectVariants \
 -V $outfile6 \
 --exclude-non-variants \
 --exclude-filtered \
 -O $outfile7

done
