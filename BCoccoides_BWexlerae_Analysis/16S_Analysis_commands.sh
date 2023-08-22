qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path Manifest.csv --output-path paired-end-demux.qza --input-format PairedEndFastqManifestPhred33

qiime demux summarize --i-data paired-end-demux.qza --o-visualization paired-end-demux.qzv

qiime dada2 denoise-paired --verbose --i-demultiplexed-seqs paired-end-demux.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 240  --p-trunc-len-r 200 --o-table table.qza --o-representative-sequences rep-seqs.qza --p-n-threads 120

qiime tools export --input-path table.qza --output-path exported-feature-table
#Exported table.qza as BIOMV210DirFmt to directory exported-feature-table

cd exported-feature-table

biom convert -i feature-table.biom -o feature-table.tsv --to-tsv

qiime metadata tabulate --m-input-file Metadata.tsv --o-visualization tabulated-metadata.qzv

qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file Metadata.tsv --verbose

qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv

qiime alignment mafft --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza --p-n-threads 25

qiime alignment mask --i-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza

qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza

qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza


#-------- Taxonomic Assignment
qiime tools import --type 'FeatureData[Sequence]' --input-path /home/vineetks/WORK/BIOAEROSOL/Bacterial/06_QIIME2_ANALYSIS/PE_READS/gg_13_5_otus/rep_set/99_otus.fasta --output-path 99_otus.qza

qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path /home/vineetks/WORK/BIOAEROSOL/Bacterial/06_QIIME2_ANALYSIS/PE_READS/gg_13_5_otus/taxonomy/99_otu_taxonomy.txt --output-path ref-taxonomy.qza

#--Extract sequences

qiime feature-classifier extract-reads --i-sequences 99_otus.qza --p-f-primer CCTACGGGNGGCWGCAG --p-r-primer WTTACCGCGGCTGCTGG --o-reads ref-seqs.qza

#--- Train the classifier

qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ref-seqs.qza --i-reference-taxonomy ref-taxonomy.qza --o-classifier classifier.qza 

#-- Test the classifier 
qiime feature-classifier classify-sklearn --i-classifier classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza 
unzip taxonomy.qza


qiime taxa collapse --i-table table.qza  --i-taxonomy taxonomy.qza --p-level 7 --o-collapsed-table table-l7.qza
qiime taxa collapse --i-table table.qza  --i-taxonomy taxonomy.qza --p-level 7 --o-collapsed-table table-l7.qza
