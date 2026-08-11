gatk --java-options "-Xmx100g -Xms100g" GenomicsDBImport \
--sample-name-map samplename_map_comp2025_chr01.txt \
--genomicsdb-workspace-path my_database20260209_Compendium2025_chr01 \
--intervals c_I

gatk --java-options "-Xmx100g" GenotypeGVCFs \
 -R /home/racz_hanna/nemeth_balint/02_bioinformatics/references/sensu_stricto_formatted.fasta \
 -G StandardAnnotation \
 -V gendb://my_database20260209_Compendium2025_chr01 \
 -O 0cohort.20260209.Compendium2025_chr01.called.vcf.gz

