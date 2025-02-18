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
conda update conda conda create -n qiime2-2024.8 --file https://data.qiime2.org/distro/core/qiime2-2024.8-py38-linux-conda.yml conda activate qiime2-2024.8
conda update conda
conda activate qiime2-2024.8

#Installing plug in
conda install -c qiime2 -c conda-forge -c bioconda -c defaults q2-dada2
conda install -c qiime2 -c conda-forge -c bioconda -c defaults q2-phylogeny

#Prepare Your Data
qiime tools import  --type 'SampleData[PairedEndSequencesWithQuality]'  --input-path /yourpath/manifestfilefinalnnn.tsv --output-path paired-end-demux.qza  --input-format PairedEndFastqManifestPhred33V2

#Demultiplexing 
qiime demux summarize --i-data /yourpath/16SRNApaired-end-demux.qza --o-visualization /yourpath/16SRNAdemux-paired-end.qzv

#visualize on qiime tools view (view.qiime2.org )
qiime tools view demux.qzv

/yourpath/aligned-rep-seqs.qza --o-masked-alignment /yourpath/masked-aligned-rep-seqs.qza --o-rooted-tree /yoorpath/rooted-tree.qza



Visualise data
qiime demux summarize --i-data /yourpath/paired-end-demux.qza --o-visualization /yourpath/demux-paired-end.qzv

Installing dada2 plugin
conda install -c qiime2 -c conda-forge -c bioconda -c defaults q2-dada2
conda install -c qiime2 -c conda-forge -c bioconda -c defaults q2-deblur,

Check quality and what to use in dada3

https://view.qiime2.org/visualization/?src=f710f4e8-0e52-4255-a8c0-e818a0a2dd2c

qiime dada2 denoise-paired --i-demultiplexed-seqs /yourapth/paired-end-demux.qza  --p-trunc-len-f 203  --p-trunc-len-r 200 --o-table /yourpath/table.qza --o-representative-sequences /yourpath/rep-seqs.qza --o-denoising-stats /yourpath/denoising-stats.qza

qiime dada2 denoise-paired --i-demultiplexed-seqs tikagenomeits.qza  --p-trunc-len-f 290 --p-trunc-len-r 290 --o-representative-sequences asv-sequences-0.qza --o-table feature-table-0.qz --o-denoising-stats dada2-stats.qza

qiime feature-table summarize --i-table feature-table-0.qz.qza  --m-sample-metadata-file metadata.tsv --o-visualization feature-table-0-summ.qzv

qiime feature-table tabulate-seqs --i-data asv-sequences-0.qza  --o-visualization asv-sequences-0-summ.qzv

#Visualise DADA results

qiime metadata tabulate  --m-input-file denoising-stats.qza  --o-visualization denoising-stats.qzv

#Rooted tree
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.qza --o-tree /yourpath/unrooted-tree.qza --o-alignment

#Taxonomy Assignment
Use a pre-trained classifier based on the 16S rRNA region (e.g., Silva, Greengenes):
wget https://data.qiime2.org/2024.8/common/silva-138-99-nb-classifier.qza

#Assigning Taxonomy
qiime feature-classifier classify-sklearn  --i-classifier silva-138-99-nb-classifier.qza  --i-reads rep-seqs.qza  --o-classification taxonomy.qza

#Visualisation
qiime metadata tabulate  --m-input-file taxonomy.qza  --o-visualization taxonomy.qzv

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

qiime diversity alpha-rarefaction \ --i-table table.qza \ --i-phylogeny rooted-tree.qza \ --p-max-depth 4000 \ --m-metadata-file sample-metadata.tsv \ --o-visualization alpha-rarefaction.qzv

Differential Abundance
qiime composition ancom  --i-table table.qza  --m-metadata-file sample-metadata.tsv  --m-metadata-column group_column  --o-visualization ancom.qzv




 
 #In Rstudio
 
 Set library path (if needed)
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
devtools::install_github("xia-lab/MicrobiomeAnalystR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
remotes::install_github("jbisanz/qiime2R")

# Install local packages (if needed)
install.packages("./install.packages.for.local/Tax4Fun/", repos = NULL, type = "source")   

install.packages("extrafont")
library(extrafont)
font_import()  # Import system fonts
loadfonts(device = "win")  # Load fonts for Windows (use "pdf" or "postscript" for other devices)

# Load CRAN packages
library(rms)
library(biomformat)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyverse)
library(vegan)
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
library(microbiomeSeq)
library(ALDEx2)
library(metagenomeSeq)
library(HMP)
library(dendextend)
library(MicrobiotaProcess)
library(Biostrings)
library(DECIPHER)
library(philr)
library(microbiomeMarker)
library(curatedMetagenomicData)
library(clusterProfiler)
library(enrichplot)
library(MicrobiomeProfiler)

# Load GitHub packages
library(qiime2R)
library(ranacapa)
library(microbiomeutilities)
library(MicrobiomeAnalystR)   
                  
# Set Working Directory                       
setwd("C:/Users/MuzhinjiN/OneDrive - University of the Free State/Documents/Metagenomic2016097426/Training/microbiome_data")  (Your own path)

#Verify working directory
getwd()

#list files in the directory
list.files()

# Importing QIIME 2 artifacts (this should be in your working directory)
ps <- qza_to_phyloseq(
  features = "table.qza",            # OTU/ASV table
  tree = "rooted-tree.qza",          # Phylogenetic tree
  taxonomy = "taxonomy.qza",         # Taxonomy assignments
  metadata = "metadata.tsv"          # Sample metadata
)

# Check PS
ps
#check the structure
str(ps)

smin <- min(sample_sums(ps))
smean <- mean(sample_sums(ps))
smax <- max(sample_sums(ps))
smin
smean
smax
table(taxa_sums(ps))[1:10]
tail(table(taxa_sums(ps)))
range(taxa_sums(physeq))

#source("microbiome_data/rarefaction.R")
p <- ggrare(ps, step = 1000, color = "SampleName", label = NULL, se = FALSE)

ggsave("rarefaction.png", p, width = 5, height = 3)

set.seed(111) # keep result reproductive
ps.rarefied = rarefy_even_depth(ps, rngseed=1, replace=F)
ps.rarefied

alpha_div <- plot_richness(ps.rarefied, x="body.site", measures=c("Observed", "Shannon")) +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

alpha_div
ggsave("alpha_div.png", richness, width = 8, height = 8)

# Define custom colors for body.site
body_site_colors <- c("gut" = "#1f77b4",  # Blue
                      "tongue" = "#ff7f0e",  # Orange
                      "right palm" = "#2ca02c", # Green
                      "left palm" = "#d62728") # Red
                      
# Choose colors you want from all in R >600 colors avaiable
colors()

# Diversity Indices
# Plot richness with custom colors and font settings
alpha_div2 <- plot_richness(ps.rarefied, x = "body.site", measures = c("Observed", "Shannon")) +
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

# Save the plot with high resolution for publication pursposes (Check the equirements of your journal)
ggsave("alpha_diversity_plot.png", alpha_div2, width = 8, height = 4, dpi = 360)                 





















# Metagenomic-analysisis

qiime tools inspect-metadata metadata1.tsv
qiime tools import --input-path metadata1.tsv --type 'SampleData[SequencesWithQuality]'  --input-format SingleEndFastqManifestPhred33V2 --output-path demux.qza
qiime demux summarize --i-data demux.qza
gunzip *.fastq.gz (all with fastq.qz)
qiime tools import --input-path metadata1.tsv --type 'SampleData[SequencesWithQuality]'  --input-format SingleEndFastqManifestPhred33V2 --output-path demux.qza
qiime dada2 denoise-single --i-demultiplexed-seqs demux-filtered.qza --p-trim-left 0 --p-trunc-len 0 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza

