REFERENCE=/home/racz_hanna/nemeth_balint/02_bioinformatics/references/sensu_stricto_formatted.fasta

for inname in ./*mdup.bam; do 

ofname1=${inname/.bam/.ss.cer.chr01.g.vcf.gz}
ofname2=${inname/.bam/.ss.cer.chr02.g.vcf.gz}
ofname3=${inname/.bam/.ss.cer.chr03.g.vcf.gz}
ofname4=${inname/.bam/.ss.cer.chr04.g.vcf.gz}
ofname5=${inname/.bam/.ss.cer.chr05.g.vcf.gz}
ofname6=${inname/.bam/.ss.cer.chr06.g.vcf.gz}
ofname7=${inname/.bam/.ss.cer.chr07.g.vcf.gz}
ofname8=${inname/.bam/.ss.cer.chr08.g.vcf.gz}
ofname9=${inname/.bam/.ss.cer.chr09.g.vcf.gz}
ofname10=${inname/.bam/.ss.cer.chr10.g.vcf.gz}
ofname11=${inname/.bam/.ss.cer.chr11.g.vcf.gz}
ofname12=${inname/.bam/.ss.cer.chr12.g.vcf.gz}
ofname13=${inname/.bam/.ss.cer.chr13.g.vcf.gz}
ofname14=${inname/.bam/.ss.cer.chr14.g.vcf.gz}
ofname15=${inname/.bam/.ss.cer.chr15.g.vcf.gz}
ofname16=${inname/.bam/.ss.cer.chr16.g.vcf.gz}

gatk --java-options "-Xmx4g" HaplotypeCaller \
 -R $REFERENCE \
 -I $inname \
 -O $ofname1 \
 -G StandardAnnotation \
 -G StandardHCAnnotation \
 --intervals c_I \
 -ploidy 2 \
 -ERC GVCF &&

gatk --java-options "-Xmx4g" HaplotypeCaller \
 -R $REFERENCE \
 -I $inname \
 -O $ofname2 \
 -G StandardAnnotation \
 -G StandardHCAnnotation \
 --intervals c_II \
 -ploidy 2 \
 -ERC GVCF &&

gatk --java-options "-Xmx4g" HaplotypeCaller \
 -R $REFERENCE \
 -I $inname \
 -O $ofname3 \
 -G StandardAnnotation \
 -G StandardHCAnnotation \
 --intervals c_III \
 -ploidy 2 \
 -ERC GVCF &&

gatk --java-options "-Xmx4g" HaplotypeCaller \
 -R $REFERENCE \
 -I $inname \
 -O $ofname4 \
 -G StandardAnnotation \
 -G StandardHCAnnotation \
 --intervals c_IV \
 -ploidy 2 \
 -ERC GVCF &&

gatk --java-options "-Xmx4g" HaplotypeCaller \
 -R $REFERENCE \
 -I $inname \
 -O $ofname5 \
 -G StandardAnnotation \
 -G StandardHCAnnotation \
 --intervals c_V \
 -ploidy 2 \
 -ERC GVCF &&

gatk --java-options "-Xmx4g" HaplotypeCaller \
 -R $REFERENCE \
 -I $inname \
 -O $ofname6 \
 -G StandardAnnotation \
 -G StandardHCAnnotation \
 --intervals c_VI \
 -ploidy 2 \
 -ERC GVCF &&

gatk --java-options "-Xmx4g" HaplotypeCaller \
 -R $REFERENCE \
 -I $inname \
 -O $ofname7 \
 -G StandardAnnotation \
 -G StandardHCAnnotation \
 --intervals c_VII \
 -ploidy 2 \
 -ERC GVCF &&

gatk --java-options "-Xmx4g" HaplotypeCaller \
 -R $REFERENCE \
 -I $inname \
 -O $ofname8 \
 -G StandardAnnotation \
 -G StandardHCAnnotation \
 --intervals c_VIII \
 -ploidy 2 \
 -ERC GVCF &&

gatk --java-options "-Xmx4g" HaplotypeCaller \
 -R $REFERENCE \
 -I $inname \
 -O $ofname9 \
 -G StandardAnnotation \
 -G StandardHCAnnotation \
 --intervals c_IX \
 -ploidy 2 \
 -ERC GVCF &&

gatk --java-options "-Xmx4g" HaplotypeCaller \
 -R $REFERENCE \
 -I $inname \
 -O $ofname10 \
 -G StandardAnnotation \
 -G StandardHCAnnotation \
 --intervals c_X \
 -ploidy 2 \
 -ERC GVCF &&

gatk --java-options "-Xmx4g" HaplotypeCaller \
 -R $REFERENCE \
 -I $inname \
 -O $ofname11 \
 -G StandardAnnotation \
 -G StandardHCAnnotation \
 --intervals c_XI \
 -ploidy 2 \
 -ERC GVCF &&

gatk --java-options "-Xmx4g" HaplotypeCaller \
 -R $REFERENCE \
 -I $inname \
 -O $ofname12 \
 -G StandardAnnotation \
 -G StandardHCAnnotation \
 --intervals c_XII \
 -ploidy 2 \
 -ERC GVCF &&

gatk --java-options "-Xmx4g" HaplotypeCaller \
 -R $REFERENCE \
 -I $inname \
 -O $ofname13 \
 -G StandardAnnotation \
 -G StandardHCAnnotation \
 --intervals c_XIII \
 -ploidy 2 \
 -ERC GVCF &&

gatk --java-options "-Xmx4g" HaplotypeCaller \
 -R $REFERENCE \
 -I $inname \
 -O $ofname14 \
 -G StandardAnnotation \
 -G StandardHCAnnotation \
 --intervals c_XIV \
 -ploidy 2 \
 -ERC GVCF &&

gatk --java-options "-Xmx4g" HaplotypeCaller \
 -R $REFERENCE \
 -I $inname \
 -O $ofname15 \
 -G StandardAnnotation \
 -G StandardHCAnnotation \
 --intervals c_XV \
 -ploidy 2 \
 -ERC GVCF &&

gatk --java-options "-Xmx4g" HaplotypeCaller \
 -R $REFERENCE \
 -I $inname \
 -O $ofname16 \
 -G StandardAnnotation \
 -G StandardHCAnnotation \
 --intervals c_XVI \
 -ploidy 2 \
 -ERC GVCF

done
