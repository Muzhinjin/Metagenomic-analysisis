
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

