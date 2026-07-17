for inname in ./*mdup.bam; do
 ofname1=${inname/.bam/.temporary.tab}

 bedtools genomecov -ibam $inname -g /home/racz_hanna/nemeth_balint/02_bioinformatics/references/sensustricto_genomefile.tab -d > $ofname1 &&

# perbase bedgraph

 ofname2=${inname/.bam/.per_base_temporary.bedGraph}

awk '{print $1"\t"$2-1"\t"$2"\t"$3}' $ofname1 > $ofname2 &&

# sorted bedgraph
 ofname3=${inname/.bam/.sort_perbase_temporary.bedGraph}

bedtools sort -g /home/racz_hanna/nemeth_balint/02_bioinformatics/references/sensustricto_genomefile.tab -i $ofname2 > $ofname3 &&

# sliding cov

 ofname5=${inname/.bam/.slidingwindow.tab}

bedtools map -a /home/racz_hanna/nemeth_balint/02_bioinformatics/references/sensustrictoslidingwindows.bed -b $ofname3 -c 4 -o median -g /home/racz_hanna/nemeth_balint/02_bioinformatics/references/sensustricto_genomefile.tab > $ofname5 &&

rm ./*temporary.*

done



