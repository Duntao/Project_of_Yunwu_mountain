######qiime2 pipeline for analysis########
######project: Rare prokaryotic sub-communities dominate the complexity of ecological networks and soil multinutrient cycling during long-term secondary #######succession in China¡¯s Loess Plateau
#######author: DUNTAO SHU#################
#######Update 12,6th,2020#################


#import data
#https://docs.qiime2.org/2019.1/tutorials/importing/#sequence-data-with-sequence-quality-information-i-e-fastq
#https://view.qiime2.org/
#https://www.jianshu.com/p/009963ac7393
#https://blog.csdn.net/weixin_34304013/article/details/86888256
#https://www.jianshu.com/p/d9e056108629
#https://cloud.tencent.com/developer/news/218445
#https://blog.csdn.net/woodcorpse/article/details/75204871
#https://blog.csdn.net/woodcorpse/article/details/86773312
#https://benjjneb.github.io/dada2/tutorial.html
#https://www.jianshu.com/p/d9e056108629
#https://forum.qiime2.org/t/qiime2-chinese-manual/838
#https://docs.qiime2.org/2019.4/plugins/available/dada2/
#https://www.jianshu.com/p/1235f6024b63
#https://forum.qiime2.org/t/self-consistency-loop-terminated-before-convergence-an-error-was-encountered-while-running-dada2-in-r-return-code-1/3532
#filt.result <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, compress = TRUE, truncQ = 2, truncLen = c(290, 260), trimLeft = c(26, 26), minLen = 20, maxN = #0, minQ = 0, maxEE = c(2, 2), rm.phix = TRUE, multithread = TRUE)
##https://github.com/YongxinLiu/QIIME2ChineseManual/blob/master/docs/QIIME%202%E6%95%99%E7%A8%8B.%2012%E6%95%B0%E6%8D%AE%E7%AD%9B%E9%80%89Filtering%20data(2020.2).md
  
#1.Import data with quality control fastq
##qiime and qiime2 
1. activate qiime2
  source activate qiime2-2020.8

2. Deactivate
   source deactivate qiime2-2020.8
#input
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path pe-33-manifest.txt \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33

# summary
qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization demux.qzv
  
  
#2. filter 
qiime dada2 denoise-paired \
--i-demultiplexed-seqs paired-end-demux.qza \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 215 \
--p-trunc-len-r 215 \
--o-representative-sequences rep-seqs-dada2.qza \
--o-table table-dada2.qza \
--o-denoising-stats dada2-stats.qza \
--p-n-threads 40

 # Visulization
qiime metadata tabulate \
  --m-input-file dada2-stats.qza \
  --o-visualization denoising-stats.qzv

#3.1 filter samples
###Prior to subsmapling, your must filter samples
###filter feature should be conducted in the qiime2
###filter sample should be performed in the R

qiime feature-table filter-features \
  --i-table table-dada2.qza \
  --p-min-samples 2 \
  --o-filtered-table sample-filtered-table.qza
  
 qiime feature-table filter-features \
  --i-table rarefied_table.qza \
  --p-min-samples 2 \
  --o-filtered-table sample-filtered-table.qza
  
#remove mitochondria
qiime taxa filter-table \
  --i-table table-dada2.qza  \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-no-mitochondria-no-chloroplast.qza  
  
#3.2 Feature
cp sample-filtered-table.qza table_filter.qza

qiime feature-table summarize \
  --i-table table_filter.qza \
  --o-visualization table_filter.qzv \
  --m-sample-metadata-file sample-metadata.txt
  
#4. represent fasta
qiime feature-table tabulate-seqs \
  --i-data rep-seqs-dada2.qza \
  --o-visualization rep-seqs.qzv

#5. tree
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-dada2.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
  
#6. calcuate diversity
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table_filter.qza \
  --p-sampling-depth 5874 \
  --m-metadata-file sample-metadata.txt \
  --output-dir core-metrics-results
  
# statistics analysis
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file sample-metadata.txt \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv
# eveness
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file sample-metadata.txt \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

#7.Beta diversity
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata.txt \
  --m-metadata-column Group \
  --o-visualization core-metrics-results/unweighted-unifrac-Group-significance.qzv \
  --p-pairwise
# unweighted-unifrac
qiime emperor plot \
  --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file sample-metadata.txt \
  --p-custom-axes timegradient \
  --o-visualization core-metrics-results/unweighted-unifrac-emperor-timegradient.qzv
# bray-curtis
qiime emperor plot \
  --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file sample-metadata.txt \
  --p-custom-axes timegradient \
  --o-visualization core-metrics-results/unweighted-unifrac-emperor-timegradient.qzv

#8. rarefraction curves
# --p-max-depth should be determined by reviewing the ¡°Frequency per sample¡± information presented in the table.qzv file
  that was created above. In general, choosing a value that is somewhere around the median frequency seems to work well. 
  but you may want to increase that value if the lines in the resulting rarefaction plot don¡¯t appear to be leveling out, 
  or decrease that value if you seem to be losing many of your samples due to low total frequencies closer to the minimum
  sampling depth than the maximum sampling depth.
# --p-max-depth middle of Frequency per sample table.qzv 
qiime diversity alpha-rarefaction \
  --i-table table_filter.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 5874 \
  --m-metadata-file sample-metadata.txt \
  --o-visualization alpha-rarefaction.qzv

#9. assign  

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path silva_132_99_16S.fna \
  --output-path 99_otus.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path majority_taxonomy_7_levels.txt \
  --output-path ref-taxonomy.qza

# Extract reference reads

qiime feature-classifier extract-reads \
  --i-sequences 99_otus.qza \
  --p-f-primer CCTAYGGGRBGCASCAG \      #341F
  --p-r-primer GGACTACNNGGGTATCTAAT \  #806R
  --o-reads ref-seqs.qza

# Train the classifier

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier Silva_132_99%_OTUs_341F-806R_classifier.qza \


#10. 
qiime feature-classifier classify-sklearn \
  --i-classifier Silva_132_99%_OTUs_341F-806R_classifier.qza \
  --i-reads rep-seqs-dada2.qza \
  --o-classification taxonomy.qza
  
qiime feature-classifier classify-sklearn \
  --i-classifier Silva_132_99%_OTUs_341F-806R_classifier.qza \
  --i-reads rep-seqs-dada2.qza \
  --o-classification taxonomy_rep-seqs.qza
  --p-threads 0

#alternative methods
  ##########
  qiime feature-classifier classify-sklearn \
  --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

  
# 
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
  
# barplot
qiime taxa barplot \
  --i-table table-dada2.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata.txt \
  --o-visualization taxa-bar-plots.qzv
  
qiime taxa barplot \
  --i-table table-dada2.qza \
  --i-taxonomy taxonomy_rep-seqs.qza \
  --m-metadata-file sample-metadata.txt \
  --o-visualization taxa-bar-plots_rep-seqs.qzv 
  
#11.filter

qiime feature-table filter-samples \
  --i-table table-dada2.qza \
  --p-min-frequency 1500 \
  --o-filtered-table sample-frequency-filtered-table.qza
  
# filter features 
qiime feature-table filter-features \
  --i-table table-dada2.qza \
  --p-min-frequency 10 \
  --o-filtered-table feature-frequency-filtered-table.qza
  
##
qiime feature-table filter-features \
  --i-table table-dada2.qza \
  --p-min-samples 2 \
  --o-filtered-table sample-filtered-table.qza

#
qiime feature-table filter-samples \
  --i-table table-dada2.qza \
  --p-min-features 10 \
  --o-filtered-table feature-contingency-filtered-table.qza
  
    
#12.export
##exporting a feature table
qiime tools export \
  --input-path table_filter.qza \
  --output-path exported-feature-table
  
  
qiime tools export \
  --input-path sample-filtered-table.qza \
  --output-path exported-feature-rarified-table
  
cd exported-feature-table
biom convert -i feature-table.biom -o otu_table.txt --to-tsv --header-key taxonomy --table-type="OTU table"

#filter 1
qiime tools export \
  --input-path feature-frequency-filtered-table.qza \
  --output-path feature-frequency-filtered-table
  
cd feature-frequency-filtered-table
biom convert -i feature-table.biom -o otu_table.txt --to-tsv --header-key taxonomy --table-type="OTU table"

#filter 2
qiime tools export \
  --input-path sample-contingency-filtered-table.qza \
  --output-path sample-contingency-filtered-table
  
cd sample-contingency-filtered-table
biom convert -i feature-table.biom -o otu_table.txt --to-tsv --header-key taxonomy --table-type="OTU table"

#filter 3
qiime tools export \
  --input-path feature-contingency-filtered-table.qza \
  --output-path feature-contingency-filtered-table
  
cd feature-contingency-filtered-table
biom convert -i feature-table.biom -o otu_table.txt --to-tsv --header-key taxonomy --table-type="OTU table"

#export tree 
qiime tools export   --input-path unrooted-tree.qza   --output-path exported-unroot-tree

qiime tools export   --input-path rooted-tree.qza   --output-path exported-root-tree

qiime tools export \
  --input-path taxonomy_rep-seqs.qza \
  --output-path exported-taxonomy-table 
  
qiime tools export \
  --input-path unrooted-tree.qza \
  --output-path exported-tree

mkdir extracted-feature-table
qiime tools extract \
  --input-path feature-table.qza \
  --output-path extracted-feature-table
  
cd core-metrics-results
qiime tools export \
  --input-path rarefied_table.qza \
  --output-path exported-rarefied-table 


#convert biom to otu table
biom convert -i feature-table.biom -o otu_table.txt --to-tsv --header-key taxonomy --table-type="OTU table"
