qstat -a 448956 (16RA straberry dadda)
module load cluster/hpc
qstat -a
qdel 
ern cluster create --name=hpc --head-node=login.hpc.ufs.ac.za --node-select=ern --email-address=MuzhinjiN@ufs.ac.za --shared-fs


module load qiime distribution=amplicon
#Import
ern jobs submit --name=16S_preplanting2024.qiime --threads=8 --memory=16gb --hours=1 --inplace --module='qiime/3.0_47df43c distribution=amplicon' --command=qiime -- demux summarize --i-data paired-end-demux.qza --o-visualization demux-paired-end.qzv
qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path Metadata.tsv --output-path pacbio-demux.qza --input-format SingleEndFastqManifestPhred33V2

ern jobs submit --name=petrus_nust.qiime --threads=16 --memory=64gb --hours=100 --inplace --module='qiime/4.0_bca2b43 distribution=amplicon' --command=qiime -- tools import --type 'SampleData[SequencesWithQuality]' --input-path Metadata.tsv --output-path pacbio-demux.qza --input-format SingleEndFastqManifestPhred33V2

ern jobs submit --name=petrus_nust.qiime --threads=16 --memory=64gb --hours=100 --inplace --module='qiime/4.0_bca2b43 distribution=amplicon' --command=qiime -- demux summarize --i-data pacbio-demux.qza --o-visualization pacbio-demux.qzv

ern jobs submit --name=petrus_nust.qiime --threads=16 --memory=64gb --hours=100 --inplace --module='qiime/4.0_bca2b43 distribution=amplicon' --command=qiime -- quality-filter q-score --i-demux pacbio-demux.qza --o-filtered-sequences filtered-pacbio-demux.qza --o-filter-stats filter-stats.qza 

ern jobs submit --name=petrus_nust.qiime --threads=16 --memory=64gb --hours=100 --inplace --module='qiime/4.0_bca2b43 distribution=amplicon' --command=qiime vsearch dereplicate-sequences --i-sequences filtered-pacbio-demux.qza --o-dereplicated-table table.qza --o-dereplicated-sequences rep-seqs.qza

ern jobs submit --name=petrus_nust.qiime --threads=16 --memory=64gb --hours=100 --inplace --module='qiime/4.0_bca2b43 distribution=amplicon' --command=qiime -- feature-classifier classify-sklearn  --i-classifier classifier-V34.qza  --i-reads rep-seqs.qza  --o-classification taxonomy.qza

ern jobs submit --name=16SRNAstraberryclassfier --threads=16 --memory=128gb --hours=100 --inplace --module='qiime/4.0_bca2b43 distribution=amplicon' --command=qiime -- feature-classifier classify-sklearn  --i-reads rep-seqs.qza --i-classifier unite_ver10_dynamic.qza --p-n-jobs  16 --o-classification vtaxonomy.qza  --verbose

ern jobs submit --name=petrus_nust.qiime --threads=16 --memory=64gb --hours=100 --inplace --module='qiime/4.0_bca2b43 distribution=amplicon' --command=qiime -- alignment mafft --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza --p-n-threads 4

ern jobs submit --name=petrus_nust.qiime --threads=16 --memory=64gb --hours=100 --inplace --module='qiime/4.0_bca2b43 distribution=amplicon' --command=qiime -- feature-classifier fit-classifier-naive-bayes --i-reference-reads silva-138-99-nb-classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza

#In linux (Command line)

#Check the quality of sequenves
fastqc sample_R1.fastq.gz sample_R2.fastq.gz -o fastqc_results/

#Trim reads with fastp
fastp -i input_R1.fastq.gz -I input_R2.fastq.gz -o trimmed_R1.fastq.gz -O trimmed_R2.fastq.gz \
      --trim_front1 10 --trim_tail1 50 --trim_front2 10 --trim_tail2 50 \
      --qualified_quality_phred 30 --length_required 100

#Removes short reads
java -jar trimmomatic PE -phred33 input_R1.fastq.gz input_R2.fastq.gz \
     trimmed_R1_paired.fastq.gz trimmed_R1_unpaired.fastq.gz \
     trimmed_R2_paired.fastq.gz trimmed_R2_unpaired.fastq.gz \
     SLIDINGWINDOW:4:30 MINLEN:100

 #Creating an environment
conda activate
conda update conda conda create -n qiime2-2024.8 --file https://data.qiime2.org/distro/core/qiime2-2024.8-py38-linux-conda.yml conda activate qiime2-2023
conda update conda
conda activate qiime2-2024.8

#Installing plug in
conda install -c qiime2 -c conda-forge -c bioconda -c defaults q2-dada2
conda install -c qiime2 -c conda-forge -c bioconda -c defaults q2-phylogeny

# Importing Data into QIIME 2
qiime tools import  --type 'SampleData[PairedEndSequencesWithQuality]'  --input-path 16SRNATomato.tsv --output-path paired-end-demux.qza  --input-format PairedEndFastqManifestPhred33V2



#Visualizing Demultiplexed Sequences 
qiime demux summarize --i-data paired-end-demux.qza --o-visualization demux-paired-end.qzv
rm stdout stderr script.sh (remoce 

#visualize on qiime tools view (view.qiime2.org )
qiime tools view demux-paired-end.qzv 

#Installing dada2 plugin
conda install -c qiime2 -c conda-forge -c bioconda -c defaults q2-dada2
conda install -c qiime2 -c conda-forge -c bioconda -c defaults q2-deblur,

#Check quality and what to use in dada3

https://view.qiime2.org/visualization/?src=f710f4e8-0e52-4255-a8c0-e818a0a2dd2c

qiime dada2 denoise-paired --i-demultiplexed-seqs paired-end-demux.qza  --p-trunc-len-f 250  --p-trunc-len-r 200 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza
ern jobs submit --name=16SRNATika_strawberry.qiime --threads=16 --memory=64gb --hours=100 --input=paired-end-demux.qza --module='qiime/4.0_bca2b43 distribution=amplicon' --command=qiime -- dada2 denoise-paired --i-demultiplexed-seqs paired-end-demux.qza  --p-trunc-len-f 250  --p-trunc-len-r 200 --p-n-threads 16  --o-table table.qza --o-representative-sequences vrep-seqs.qza --o-denoising-stats vdenoising-stats.qza


qiime feature-table summarize --i-table feature-table-0.qz.qza  --m-sample-metadata-file metadata.tsv --o-visualization feature-table-0-summ.qzv

qiime feature-table tabulate-seqs --i-data asv-sequences-0.qza  --o-visualization asv-sequences-0-summ.qzv

#Visualise DADA results

qiime metadata tabulate  --m-input-file denoising-stats.qza  --o-visualization denoising-stats.qzv

#Rooted tree
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.qza --o-alignment alligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seq.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

#Training feature classifiers with q2-feature-classifier¶
mkdir training-feature-classifiers
cd training-feature-classifiers

#Two elements are required for training the classifier: the reference sequences and the corresponding taxonomic classifications.

For ITS download from UNITE and for 16SRNA use SILVA or Greengenes

Import the data into QIIME 2 Artifacts. The reference taxonomy file (.txt) is a tab-separated (TSV) file without a header, so specify  HeaderlessTSVTaxonomyFormat as the source format since the default source format requires a header.

qiime tools import --type 'FeatureData[Sequence]' --input-path sh_refs_qiime_ver9_99_25.07.2023.fasta  --output-path UNITE.REF.qza
qiime tools import   --type 'FeatureData[Taxonomy]'  --input-format HeaderlessTSVTaxonomyFormat --input-path sh_taxonomy_qiime_ver9_99_25.07.2023.txt --output-path UNITEref-taxonomy.qza

#Extract reference reads (for 16SRNA include primers and for ITS skip)

#16SRNA

qiime feature-classifier extract-reads \
  --i-sequences 85_otus.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-trunc-len 120 \
  --p-min-length 100 \
  --p-max-length 400 \
  --o-reads ref-seqs.qza

#ITS

qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads UNITE.REF.qza --i-reference-taxonomy UNITEref-taxonomy.qza --o-classifier unite-classifier.qza

#Training the classifier using reference read and taxonomy

qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads UNITE.REF.qza --i-reference-taxonomy UNITEref-taxonomy.qza --o-classifier unite-classifier.qza

# Taxonomy asignment 
qiime feature-classifier classify-sklearn  --i-classifier unite-classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza

#Taxonomy Assignment
Use a pre-trained classifier based on the 16S rRNA region (e.g., Silva, Greengenes):
wget https://data.qiime2.org/2024.8/common/silva-138-99-nb-classifier.qza

#Assigning Taxonomy
qiime feature-classifier classify-sklearn  --i-classifier silva-138-99-nb-classifier.qza  --i-reads rep-seqs.qza  --o-classification taxonomy.qza

#Visualisation
qiime metadata tabulate  --m-input-file taxonomy.qza  --o-visualization taxonomy.qzv

#Remove features that contain mitochondria or chloroplast
qiime taxa filter-table --i-table filtered-sequences/feature-frequency-filtered-table.qza --i-taxonomy taxonomy.qza --p0-include p__ --p-exclude mitochondria,chloroplast --o-filtered-table filtered-sequences/table-with-phyla-no-mitochondria-chloroplast.qza
#Remove Archaea
qiime taxa filter-table --i-table filtered-sequences/table-with-phyla-no-mitochondria-chloroplast.qza --i-taxonomy taxonomy.qza --p-exclude "k__Archaea" \
--o-filtered-table filtered-sequences/table-with-phyla-no-mitochondria-chloroplasts-archaea.qza
#Filter Eukaryota
qiime taxa filter-table --i-table table-with-phyla-no-mitochondria-chloroplasts-archaea.qza --i-taxonomy taxonomy.qza --p-exclude "k__Eukaryota" --o-filtered-table table-with-phyla-no-mitochondria-chloroplasts-archaea-eukaryota.qza

#Filter from sequences
Remove features that contain mitochondria or chloroplast 
qiime taxa filter-seqs --i-sequences rep-seqs.qza --i-taxonomy taxonomy.qza --p-include p__ --p-exclude mitochondria,chloroplast --o-filtered-sequences rep-seqs-with-phyla-no-mitochondria-chloroplast.qza


#Remove Archaea
qiime taxa filter-seqs --i-sequences rep-seqs-with-phyla-no-mitochondria-chloroplast.qza --i-taxonomy taxonomy.qza --p-exclude "k__Archaea" --o-filtered-sequences filtered-sequences/rep-seqs-with-phyla-no-mitochondria-chloroplasts-archaea.qza

#Let’s rename the final filtered file to proceed
mv table-with-phyla-no-mitochondria-chloroplasts-archaea-eukaryota.qza filtered-table2.qza
mv rep-seqs-with-phyla-no-mitochondria-chloroplasts-archaea-eukaryota.qza filtered-rep-seqs2.qza

#Visualize taxonomic classifications

qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy.qza --m-metadata-file ITSTomato.tsv --o-visualization taxa-bar-plots.qzv

#Visualise

qiime tools view taxa-bar-plots.qzv

#Generate Diversity Metrics

#Create a phylogenetic tree:

qiime phylogeny align-to-tree-mafft-fasttree  --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza  --o-masked-alignment masked-aligned-rep-seqs.qza  --o-tree unrooted-tree.qza  --o-rooted-tree rooted-tree.qza

#Calculate diversity metrics

qiime diversity core-metrics-phylogenetic  --i-phylogeny rooted-tree.qza  --i-table table.qza --p-sampling-depth 1109  --m-metadata-file sample-metadata.tsv --output-dir core-metrics-results

#FeatureTable and FeatureData summaries
qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file manifestfilefinalnnnnf1.tsv

qiime feature-table tabulate-seqs  --i-data 16SRNArep-seqs.qza --o-visualization rep-seqs.qzv

qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.qza --p-sampling-depth 3920 --m-metadata-file manifestfilefinalnnnnf1.tsv --output-dir core-metrics-results


#Visualize diversity results Alpha and beta diversity analysis¶

qiime diversity alpha-rarefaction  --i-table table.qza  --i-phylogeny rooted-tree.qza  --p-max-depth 4000  --m-metadata-file sample-metadata.tsv  --o-visualization alpha-rarefaction.qzv

qiime diversity alpha-group-significance  --i-alpha-diversity diversity-core-metrics-phylogenetic/faith_pd_vector.qza  --m-metadata-file sample-metadata.tsv --o-visualization faith-pd-group-significance.qzv
qiime diversity alpha-group-significance \
  --i-alpha-diversity diversity-core-metrics-phylogenetic/evenness_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization evenness-group-significance.qzv

Differential Abundance
qiime composition ancom  --i-table table.qza  --m-metadata-file sample-metadata.tsv  --m-metadata-column group_column  --o-visualization ancom.qzv

# Convert taxonomy.qza to taxonomy.txt
qiime tools export --input-path taxonomy.qza --output-path exported-taxonomy
mv exported-taxonomy/taxonomy.tsv taxonomy.txt

# Convert table.qza to table.txt
qiime tools export --input-path table.qza --output-path exported-table
biom convert -i exported-table/feature-table.biom -o table.txt --to-tsv

# Convert rooted-tree.qza to rooted-tree.txt
qiime tools export --input-path rooted-tree.qza --output-path exported-tree
mv exported-tree/tree.nwk rooted-tree.txt
#In Rstudio
x<-c(1,2,3,4,5,6,7,8,9,10) 
y<-c(3,6,9,10,13,30,20,100,220,100) 
plot(x,y) 
plot(x,y,col="red") # many functions have additional parameters 
boxplot(x,y,col="red")
boxplot(x,y,col=c("hotpink", "yellow")) 
boxplot(x,y,col=c("hotpink", "yellow"),main="Lec2") # can use the ? sign to find out more about function ?boxplot

ggsave(filename = "foo3.png",width = 90, height = 180, unit = "mm", dpi = 300)

# Set library path (if needed)
# .libPaths(new = "C:/Users/MuzhinjiN/AppData/Local/Programs/R/R-4.4.1/library")
.libPaths()

# Install CRAN packages
install.packages(c(
  "tidytree", "RNeXML", "remotes", "biomformat", "adespatial", "ape", "devtools", 
  "ggdendro", "gridExtra", "knitr", "MicrobeR", "pander", "plotly", "png", 
  "tidyverse", "vegan", "reshape", "picante", "rms", "nlme", "ggnewscale", 
  "patchwork", "pals", "ggstar", "ggraph", "ggpubr", "ggtree", "treeio", 
  "tidytree", "forcats", "ggtreeExtra", "ggClusterNet", "ggraph", "grid", 
  "gridExtra", "knitr", "png", "ggdendro", "ggpubr", "ggtree", "ggtreeExtra", 
  "treeio", "tidytree", "forcats", "ggstar", "patchwork", "pals", "ggraph"
), dependencies = TRUE, force = TRUE)

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2", "phyloseq", "microbiomeSeq", "microbiome", "ALDEx2", "metagenomeSeq", 
  "HMP", "dendextend", "adespatial", "MicrobiotaProcess", "GO.db", "impute", 
  "Biostrings", "DECIPHER", "philr", "microbiomeMarker", "curatedMetagenomicData", 
  "clusterProfiler", "enrichplot", "MicrobiomeProfiler"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ShortRead")

# Install GitHub packages
devtools::install_github(repo = "malucalle/selbal")
devtools::install_github("adw96/breakaway")
devtools::install_github("taowenmicro/EasyStat", upgrade = "never")
devtools::install_github("taowenmicro/ggClusterNet", upgrade = "never")
devtools::install_github("zdk123/SpiecEasi")
devtools::install_github("fjossandon/Tax4Fun2", force = TRUE)
devtools::install_github("kylebittinger/qiimer", force = TRUE)
devtools::install_github("joey711/biom", upgrade = "never", force = TRUE)
devtools::install_github("gauravsk/ranacapa")
devtools::install_github("umerijaz/microbiomeSeq")
devtools::install_github("microsud/microbiomeutilities")
remotes::install_github("jbisanz/qiime2R")
1
# Install local packages (if needed)
install.packages("./install.packages.for.local/Tax4Fun/", repos = NULL, type = "source")   

install.packages("extrafont")
library(extrafont)
font_import()  # Import system fonts
loadfonts(device = "win")  # Load fonts for Windows (use "pdf" or "postscript" for other devices)

# Load CRAN packages
library(ShortRead)
library(rms)
library(biomformat)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyverse)
library(vegan)    # vegetation analysis
library(RColorBrewer)
library(reshape2)
library(ggtreeExtra)
library(forcats)
library(ggstar)
library(patchwork)
library(pals)
library(ggpubr)
library(ggtree)
library(treeio)
library(tidytree)
library(ggraph)
library(grid)
library(gridExtra)
library(knitr)
library(png)
library(ggdendro)
library(ggpubr)
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(tidytree)
library(forcats)
library(ggstar)
library(patchwork)
library(pals)
library(ggraph)

# Load Bioconductor packages
library(DESeq2)
library(phyloseq)
library(microbiome)
library(ALDEx2)
library(metagenomeSeq)
library(HMP)
library(dendextend)
library(MicrobiotaProcess)
library(Biostrings)
library(DECIPHER)
library(philr)
library(curatedMetagenomicData)
library(clusterProfiler)
library(enrichplot)
library(MicrobiomeProfiler)

# Load GitHub packages
library(qiime2R)
library(ranacapa)
library(microbiomeutilities)
                     
## Set Working Directory                       
setwd("C:/Users/MuzhinjiN/OneDrive - University of the Free State/Documents/Metagenomic2016097426/Training/microbiome_data")
getwd()

# GUnzipping fastq files

# Install packages
install.packages("R.utils")
library(R.utils)

# Define file paths
fastq_gzR1 <- "C:/Users/MuzhinjiN/OneDrive - University of the Free State/Documents/Metagenomic2016097426/Training/microbiome_data/303-soil-V3V4_S125_L001_R1_001.fastq.gz"
fastq_gzR2 <- "C:/Users/MuzhinjiN/OneDrive - University of the Free State/Documents/Metagenomic2016097426/Training/microbiome_data/303-soil-V3V4_S125_L001_R2_001.fastq.gz"

#Prepare file for unzipping
fastq_unzippedR1 <- "C:/Users/MuzhinjiN/OneDrive - University of the Free State/Documents/Metagenomic2016097426/Training/microbiome_data/303-soil-V3V4_S125_L001_R1_001.fastq"
fastq_unzippedR2 <- "C:/Users/MuzhinjiN/OneDrive - University of the Free State/Documents/Metagenomic2016097426/Training/microbiome_data/303-soil-V3V4_S125_L001_R2_001.fastq"

# Unzip the file
gunzip(fastq_gzR1, destname = fastq_unzippedR1, overwrite = TRUE)
gunzip(fastq_gzR2, destname = fastq_unzippedR2, overwrite = TRUE)

#Remove files from the folder if you don't need it
file_path <- "C:/Users/MuzhinjiN/OneDrive - University of the Free State/Documents/Metagenomic2016097426/Training/microbiome_data/303-root-V3V4_S150_L001_R2_001.fastq.gz"
file.remove(file_path)
list.files()

# Load necessary libraries
library(ShortRead)
library(vctrs)

# Define the path to your FASTQ file
fastq_file <- "C:/Users/MuzhinjiN/OneDrive - University of the Free State/Documents/Metagenomic2016097426/Training/microbiome_data/303-soil-V3V4_S125_L001_R1_001.fastq"

# Check if the file exists
if (!file.exists(fastq_file)) {
  stop("FASTQ file not found!")
}

# Read the FASTQ file
fq <- tryCatch({
  readFastq(fastq_file)
}, error = function(e) {
  stop("Error reading FASTQ file: ", e$message)
})

# Extract sequences and convert to character vector
sequences <- sread(fq)
seq_char <- as.character(sequences)

# Check the structure of the FASTQ object
print(fq)

# Group identical sequences
group_ids <- vctrs::vec_group_id(seq_char)
print("Group IDs for the first 5 sequences:")
print(group_ids[1:5])

# View the first few sequences and quality scores
print("First 5 sequences:")
print(sread(fq)[1:5])
print("Quality scores of the first 5 reads:")
print(quality(fq)[1:5])

# Extract and display components of the first read
identifier <- as.character(vctrs::vec_group_id(fq)[1])
sequence <- as.character(sread(fq)[1])
plus_line <- "+"
quality_scores <- as.character(PhredQuality(quality(fq)[1]))

cat("First read components:\n")
cat("Identifier: ", identifier, "\n")
cat("Sequence: ", sequence, "\n")
cat("Plus line: ", plus_line, "\n")
cat("Quality scores: ", quality_scores, "\n")

# Importing QIIME 2 artifacts
ps <- qza_to_phyloseq(
  features = "table.qza",            # ASV table
  tree = "rooted-tree.qza",          # Phylogenetic tree
  taxonomy = "taxonomy.qza",         # Taxonomy assignments
  metadata = "metadata.tsv"          # Sample metadata
)

ps
#view structure
str(ps)

#Check the metadata

metadata <- read.table("metadata.tsv", , sep='\t', header=T, row.names=1, comment="")
metadata <- metadata[-1,] # remove the second line that specifies the data type
metadata
smin <- min(sample_sums(ps))
smean <- mean(sample_sums(ps))
smax <- max(sample_sums(ps))
smin
smean
smax
table(taxa_sums(ps))[1:10]
tail(table(taxa_sums(ps)))
range(taxa_sums(ps))

#source("microbiome_data/rarefaction.R")
rarefaction <- ggrare(ps, step = 1000, color = "SampleName", label = NULL, se = FALSE)

ggsave("rarefaction.png", p, width = 5, height = 3)

set.seed(111) # keep result reproductive
ps.rarefied = rarefy_even_depth(ps, rngseed=1, replace=F)
ps.rarefied

alpha_div <- plot_richness(ps.rarefied, x="body.site", measures=c("Observed", "Shannon", "Chao1")) +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

alpha_div
ggsave("alpha_div.png", alpha_div, width = 8, height = 8)

# Define custom colors for body.site
body_site_colors <- c("gut" = "#1f77b4",  # Blue
                      "tongue" = "#ff7f0e",  # Orange
                      "right palm" = "#2ca02c", # Green
                      "left palm" = "#d62728") # Red

colors()

# Plot richness with custom colors and font settings
alpha_div2 <- plot_richness(ps.rarefied, x = "body.site", measures = c("Observed", "Shannon", "Chao1")) +
  geom_boxplot(aes(fill = body.site), alpha = 0.7) +  # Fill boxplots by body.site
  theme_classic() +
  theme(
    strip.background = element_blank(),  # Remove strip background
    axis.text.x.bottom = element_text(angle = 45, hjust = 0.7, vjust = 0.5, color = "black", size = 12, family = "Times New Roman"),  # Customize x-axis text
    axis.text.y = element_text(color = "black", size = 12, family = "Times New Roman"),  # Customize y-axis text
    axis.title = element_text(color = "black", size = 12, family = "Times New Roman"),  # Customize axis titles
    legend.text = element_text(color = "black", size = 12, family = "Times New Roman"),  # Customize legend text
    legend.title = element_text(color = "black", size = 12, family = "Times New Roman")  # Customize legend title
  ) +
  scale_fill_manual(values = body_site_colors) +  # Use custom colors for body.site
  labs(x = "Body Site", y = "Alpha Diversity Measure")  # Add axis labels

# Display the plot with the legend
print(alpha_div2)

# Display the plot without the legend
alpha_div2_no_legend <- alpha_div2 + theme(legend.position = "none")
print(alpha_div2_no_legend)

# Save the plot with high resolution
ggsave("alpha_div2_no_legend.png", alpha_div2, width = 8, height = 4, dpi = 360)                 



# Plot Shannon only

# Plot richness with custom colors and font settings
alpha_div3 <- plot_richness(ps.rarefied, x = "body.site", measures = c("Shannon")) +
  geom_boxplot(aes(fill = body.site), alpha = 0.7) +  # Fill boxplots by body.site
  theme_classic() +
  theme(
    strip.background = element_blank(),  # Remove strip background
    axis.text.x.bottom = element_text(angle = 45, hjust = 0.7, vjust = 0.8, color = "black", size = 12, family = "Times New Roman"),  # Customize x-axis text
    axis.text.y = element_text(color = "black", size = 12, family = "Times New Roman"),  # Customize y-axis text
    axis.title = element_text(color = "black", size = 12, family = "Times New Roman"),  # Customize axis titles
    legend.text = element_text(color = "black", size = 12, family = "Times New Roman"),  # Customize legend text
    legend.title = element_text(color = "black", size = 12, family = "Times New Roman")  # Customize legend title
  ) +
  scale_fill_manual(values = body_site_colors) +  # Use custom colors for body.site
  labs(x = "Body Site", y = "Alpha Diversity Measure")  # Add axis labels

# Display the plot with the legend
print(alpha_div3)

# Display the plot without the legend
alpha_div3_no_legend <- alpha_div3 + theme(legend.position = "none")
print(alpha_div3_no_legend)

# Save the plot with high resolution
ggsave("alpha_div2_no_legend.png", alpha_div2, width = 8, height = 6, dpi = 360)  

# Estimate alpha diversity

rich = estimate_richness(ps.rarefied, measures = c("Observed", "Shannon", "Chao"))

# Ensure body.site is a factor

wilcox.observed <- pairwise.wilcox.test(rich$Observed,
                                        sample_data(ps.rarefied)$body.site,
                                        p.adjust.method = "BH")

tab.observed <- wilcox.observed$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()

# Print result
tab.observed

write.csv(tab.observed, file = "tab.observed.csv")

#Shannon
wilcox.shannon <- pairwise.wilcox.test(rich$Shannon,
                                       sample_data(ps.rarefied)$body.site,
                                       p.adjust.method = "BH")
tab.shannon <- wilcox.shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()

#Print table shannon
tab.shannon

write.csv(tab.shannon, file = "tab.shannon.csv")

#Beta diversity

ord1 <- ordinate(ps.rarefied, method="PCoA", distance="bray")
colourCount_ord = length(unique(metadata$body.site))

p1 <- plot_ordination(ps.rarefied, ord1) +
  geom_point(aes(color = as.factor(body.site)), size = 8, alpha=0.75) +
  scale_fill_manual(values = body_site_colors) +  # Use custom colors for body.site 
  labs(color = "Body Site") +
  theme_bw() +
  theme(legend.position = "right")+
  theme(
    strip.background = element_blank(),  # Remove strip background
    axis.text.y = element_text(color = "black", size = 12, family = "Times New Roman"),  # Customize y-axis text
    axis.text.x = element_text(color = "black", size = 12, family = "Times New Roman"),  # Customize y-axis text
    axis.title = element_text(color = "black", size = 12, family = "Times New Roman"),  # Customize axis titles
    legend.text = element_text(color = "black", size = 12, family = "Times New Roman"),  # Customize legend text
    legend.title = element_text(color = "black", size = 12, family = "Times New Roman")  # Customize legend title
  )
p1

#customise colors
p2 <- plot_ordination(ps.rarefied, ord1) +
  geom_point(aes(color = as.factor(body.site)), size = 8, alpha=0.75) +
  scale_color_manual(values = body_site_colors) +  # Use custom colors for body.site 
  labs(color = "Body Site") +
  theme_bw() +
  theme(legend.position = "right")+
  theme(
    strip.background = element_blank(),  # Remove strip background
    axis.text.y = element_text(color = "black", size = 12, family = "Times New Roman"),  # Customize y-axis text
    axis.text.x = element_text(color = "black", size = 12, family = "Times New Roman"),  # Customize y-axis text
    axis.title = element_text(color = "black", size = 12, family = "Times New Roman"),  # Customize axis titles
    legend.text = element_text(color = "black", size = 12, family = "Times New Roman"),  # Customize legend text
    legend.title = element_text(color = "black", size = 12, family = "Times New Roman")  # Customize legend title
  )
p2

ggsave("betadiversity.png", p2, width = 12, height = 10, dpi =360)


# Calculate distance matrix (Bray-Curtis)
dist_matrix <- phyloseq::distance(ps, method = "bray")

# Perform PCoA
pcoa <- ordinate(ps, method = "PCoA", distance = dist_matrix)

# Plot PCoA
p2 <- plot_ordination(ps, pcoa, color = "body.site") + 
  geom_point(size = 3) + 
  theme_minimal()
p2

# Assuming 'ps' is your phyloseq object and 'pcoa' is your PCoA result
p3 <- plot_ordination(ps, pcoa, color = "body.site") + 
  geom_point(size = 3, shape = 17) +  # Use triangles (shape = 17)
  scale_color_manual(values = body_site_colors) +
  theme_minimal() +
  theme(
    text = element_text(size = 12, color = "black"),  # Set font size and color for all text
    axis.title = element_text(color = "black"),       # Axis titles in black
    axis.text = element_text(color = "black"),        # Axis labels in black
    legend.text = element_text(color = "black"),      # Legend text in black
    legend.title = element_text(color = "black")      # Legend title in black
  )
p3

# Perform RDA (Redundancy Analysis)
rda_result <- ordinate(ps, method = "RDA", distance = "bray", formula = ~ year)

# Plot RDA
p4 <- plot_ordination(ps, rda_result, color = "body.site") + 
  geom_point(size = 3, shape = 17) +  # Use triangles (shape = 17)
  scale_color_manual(values = body_site_colors) +
  theme_minimal() +
theme(
    text = element_text(size = 12, color = "black"),  # Set font size and color for all text
    axis.title = element_text(color = "black"),       # Axis titles in black
    axis.text = element_text(color = "black"),        # Axis labels in black
    legend.text = element_text(color = "black"),      # Legend text in black
    legend.title = element_text(color = "black")      # Legend title in black
p4

sample_variables(ps)

# Perform RDA (Redundancy Analysis)
rda_result <- ordinate(ps, method = "RDA", distance = "bray", formula = ~ year)

# Plot RDA
p4 <- plot_ordination(ps, rda_result, color = "body.site") + 
  geom_point(size = 3, shape = 17) +  # Use triangles (shape = 17)
  scale_color_manual(values = body_site_colors) +
  theme_minimal() +
  theme(
    text = element_text(size = 12, color = "black"),  # Set font size and color for all text
    axis.title = element_text(color = "black"),       # Axis titles in black
    axis.text = element_text(color = "black"),        # Axis labels in black
    legend.text = element_text(color = "black"),      # Legend text in black
    legend.title = element_text(color = "black")      # Legend title in black
  ) +
  labs(
    x = "RDA Axis 1",
    y = "RDA Axis 2",
    title = "Redundancy Analysis (RDA) of Microbial Communities",
    color = "Body Site"
  )

# Display the plot
print(p4)

# Save the plot
ggsave("rda_plot.png", p4, width = 8, height = 6, dpi = 300)

# Perform RDA with conditioning variable of antibiotic resistance
rda_result <- ordinate(
  ps, 
  method = "RDA", 
  distance = "bray", 
  formula = ~ year + Condition(reported.antibiotic.usage)
)
# Plot RDA
p5 <- plot_ordination(ps, rda_result, color = "body.site") + 
  geom_point(size = 3, shape = 17) +  # Use triangles (shape = 17)
  scale_color_manual(values = body_site_colors) +
  theme_minimal() +
  theme(
    text = element_text(size = 12, color = "black"),  # Set font size and color for all text
    axis.title = element_text(color = "black"),       # Axis titles in black
    axis.text = element_text(color = "black"),        # Axis labels in black
    legend.text = element_text(color = "black"),      # Legend text in black
    legend.title = element_text(color = "black")      # Legend title in black
  ) +
  labs(
    x = "RDA Axis 1",
    y = "RDA Axis 2",
    title = "Redundancy Analysis (RDA) of Microbial Communities",
    color = "Body Site"
  )

# Display the plot
print(p5)

# Save the plot
ggsave("rda_plot_with_conditioning.png", p4, width = 8, height = 6, dpi = 300)

# Relative abundance
physeq_rel <- transform_sample_counts(ps, function(x) x / sum(x))


# Aggregate data to Kingdom level
physeq_kingdom <- tax_glom(physeq_rel, taxrank = "Kingdom")

# Aggregate data to Phylum level
physeq_phyla <- tax_glom(physeq_rel, taxrank = "Phylum")

# Aggregate data to Family level
physeq_family <- tax_glom(physeq_rel, taxrank = "Family")

# Aggregate data to Order level
physeq_order <- tax_glom(physeq_rel, taxrank = "Order")

# Aggregate data to Class level
physeq_class <- tax_glom(physeq_rel, taxrank = "Class")

# Aggregate data to Genus level
physeq_genus <- tax_glom(physeq_rel, taxrank = "Genus")

# Aggregate data to Species level
physeq_species <- tax_glom(physeq_rel, taxrank = "Species")

# Top 30 Phyla
top30_phyla <- names(sort(taxa_sums(physeq_phyla), decreasing = TRUE))[1:30]
physeq_phyla_top30 <- prune_taxa(top30_phyla, physeq_phyla)

# Top 30 Families
top30_family <- names(sort(taxa_sums(physeq_family), decreasing = TRUE))[1:30]
physeq_family_top30 <- prune_taxa(top30_family, physeq_family)

# Top 30 Genera
top30_genus <- names(sort(taxa_sums(physeq_genus), decreasing = TRUE))[1:30]
physeq_genus_top30 <- prune_taxa(top30_genus, physeq_genus)

# Top 30 Species
top30_species <- names(sort(taxa_sums(physeq_species), decreasing = TRUE))[1:30]
physeq_species_top30 <- prune_taxa(top30_species, physeq_species)

# Plot relative abundance at Phyla level
plot_phyla <- plot_bar(physeq_phyla_top30, fill = "Phylum") +
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Top 30 Phylum - Relative Abundance")
plot_phyla

# Plot relative abundance at Family level
plot_family <- plot_bar(physeq_family_top30, fill = "Family") +
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Top 30 Families - Relative Abundance")
plot_family

# Ensure you have the extrafont package for Times New Roman
if (!requireNamespace("extrafont", quietly = TRUE)) {
  install.packages("extrafont")
}
library(extrafont)
custom_colors <- c(
  # Basic colors
  "red", "blue", "green", "purple", "orange", "pink", "cyan", "yellow", "brown", "black",
  
  # Dark shades
  "darkblue", "darkgreen", "darkorange", "darkcyan", "darkorchid4", "violetred4",
  
  # Light shades
  "lightcoral", "lightblue", "lightgreen", "lightpink", "lightcyan", "lightyellow",
  
  # Pastel colors
  "thistle1", "seagreen1", "rosybrown3", "peachpuff3", "mistyrose", "lemonchiffon",
  
  # Miscellaneous
  "gold", "mediumaquamarine", "tan4", "orchid3", "wheat2"
)

# Create the bar plot
plot_family <- plot_bar(physeq_family_top30, fill = "Family") +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  
  # Customize axis text and titles
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12, family = "Times New Roman", color = "black"),
    axis.text.y = element_text(size = 12, family = "Times New Roman", color = "black"),
    axis.title.x = element_text(size = 12, family = "Times New Roman", color = "black"),
    axis.title.y = element_text(size = 12, family = "Times New Roman", color = "black"),
    plot.title = element_text(size = 14, family = "Times New Roman", color = "black", hjust = 0.5),
    legend.text = element_text(size = 10, family = "Times New Roman", color = "black"),
    legend.title = element_text(size = 12, family = "Times New Roman", color = "black")
  ) +
  
  # Add a custom color palette
  scale_fill_manual(values = custom_colors) +
  
  # Add a title
  labs(title = "Top 20 Families - Relative Abundance", x = "Samples", y = "Relative Abundance") +
  
  # Split the legend into two panels
  guides(fill = guide_legend(ncol = 2, title.position = "top"))


# Print the plots
print(plot_phyla)
print(plot_family)
print(plot_genus)
print(plot_species)

# Convert phyloseq object to a data frame for plotting at Family level
family_df <- psmelt(physeq_family_top30)  # Melt phyloseq object to a long format data frame

# Convert phyloseq object to a data frame for plotting at phyla level
phyla_df <- psmelt(physeq_phyla_top30)  # Melt phyloseq object to a long format data frame

# Convert phyloseq object to a data frame for plotting at Family level
genus_df <- psmelt(physeq_genus_top30)  # Melt phyloseq object to a long format data frame

# Convert phyloseq object to a data frame for plotting at Family level
species_df <- psmelt(physeq_species_top30)  # Melt phyloseq object to a long format data frame

# Plotting the relative abundance at Phylum level, faceted by 'body.site'
pbacteriaphylumA <- ggplot(phyla_df, aes(x = body.site, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ year, scales = "free_x") +  # Facet by 'year'
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Increase X-axis font size
    axis.text.y = element_text(size = 12),  # Increase Y-axis font size
    strip.text = element_text(size = 12),  # Increase facet label font size
    legend.title = element_text(size = 12),  # Increase legend title size
    legend.text = element_text(size = 12)  # Increase legend text size
  ) +
  labs(
    x = "body.site",  # Label for X-axis
    y = "Relative Abundance",  # Label for Y-axis
    title = "A"
  ) +
  scale_fill_manual(values = custom_colors) + # Use a color palette for Family
  guides(fill = guide_legend(ncol = 2))

# Print the plot
print(pbacteriaphylumA)

# Plotting the relative abundance at Family level, faceted by 'Sampling_part'
phylumA <- ggplot(family_df, aes(x = treatment, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Sampling_part, scales = "free_x") +  # Facet by 'SamplingPart'
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Increase X-axis font size
    axis.text.y = element_text(size = 12),  # Increase Y-axis font size
    strip.text = element_text(size = 12),  # Increase facet label font size
    legend.title = element_text(size = 12),  # Increase legend title size
    legend.text = element_text(size = 12)  # Increase legend text size
  ) +
  labs(
    x = "Treatment",  # Label for X-axis
    y = "Relative Abundance",  # Label for Y-axis
    title = "A"
  ) +
  scale_fill_manual(values = palletet) + # Use a color palette for Family
  guides(fill = guide_legend(ncol = 2))

# Print the plot
print(Phylum)

library(VennDiagram)

vennlist <- get_vennlist(obj=ps, factorNames="body.site")
vennp <- venn.diagram(vennlist,
                      height=5,
                      width=5, 
                      filename=NULL, 
                      fill=c("cyan","red", "grey", "#FD9347"),
                      cat.col=c("cyan", "red", "grey",  "#FD9347"),
                      alpha = 1, 
                      fontfamily = "serif",
                      fontface = "bold",
                      cex = 1.0,
                      cat.cex = 1.3,
                      cat.default.pos = "outer",
                      cat.dist=0.1,
                      margin = 0.1, 
                      lwd = 3,
                      lty ='dotted',
                      imagetype = "svg")

grid::grid.draw(vennp)


# Load necessary libraries
library(phyloseq)
library(ggtree)
library(ggplot2)
library(dplyr)

# Step 1: Aggregate data at the genus level
ps_genus <- tax_glom(ps, taxrank = "Genus")

# Step 2: Subset data for a specific body.site
body_site <- "gut"  # Replace with your desired body.site
ps_subset <- subset_samples(ps_genus, body.site == body_site)

# Step 3: Identify the top 30 genera

# Calculate relative abundance
ps_rel_abund <- transform_sample_counts(ps_subset, function(x) x / sum(x))
genus_abund <- taxa_sums(ps_rel_abund)
top_genera <- names(sort(genus_abund, decreasing = TRUE)[1:30])
                    
# Step 4: Extract the phylogenetic tree
tree <- phy_tree(ps_subset)
                    
# Step 5: Prune the tree to the top 30 genera
pruned_tree <- keep.tip(tree, top_genera)
                    
# Step 6: Plot the phylogenetic tree
  ggtree(pruned_tree) + geom_tiplab(size = 3, color = "blue") +  # Add genus labels
                      theme_tree2() +
                      ggtitle(paste("Phylogenetic Tree of Top 30 Genera in", body_site))

 
 


  
                  


               





















# Metagenomic-analysisis

qiime tools inspect-metadata metadata1.tsv
qiime tools import --input-path metadata1.tsv --type 'SampleData[SequencesWithQuality]'  --input-format SingleEndFastqManifestPhred33V2 --output-path demux.qza
qiime demux summarize --i-data demux.qza
gunzip *.fastq.gz (all with fastq.qz)
qiime tools import --input-path metadata1.tsv --type 'SampleData[SequencesWithQuality]'  --input-format SingleEndFastqManifestPhred33V2 --output-path demux.qza
qiime dada2 denoise-single --i-demultiplexed-seqs demux-filtered.qza --p-trim-left 0 --p-trunc-len 0 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza

