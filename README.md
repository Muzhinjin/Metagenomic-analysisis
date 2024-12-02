# Metagenomic-analysisis

qiime tools inspect-metadata metadata1.tsv
qiime tools import --input-path metadata1.tsv --type 'SampleData[SequencesWithQuality]'  --input-format SingleEndFastqManifestPhred33V2 --output-path demux.qza
qiime demux summarize --i-data demux.qza
gunzip *.fastq.gz (all with fastq.qz)
qiime tools import --input-path metadata1.tsv --type 'SampleData[SequencesWithQuality]'  --input-format SingleEndFastqManifestPhred33V2 --output-path demux.qza
qiime dada2 denoise-single --i-demultiplexed-seqs demux-filtered.qza --p-trim-left 0 --p-trunc-len 0 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza

