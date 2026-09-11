# setup.R

# packages --------

# essentials
library(tidyverse)

# I/O
library(readxl)

# visualization
library(cowplot)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(viridis)
library(RColorBrewer)
library(scattermore)
library(plotly)
library(ggtext)

# tables
library(DT)

# other analysis
library(cluster)

# custom packages
library(RFunctions)

# Chapters throughout the book call calculate_proportion() unqualified. Rsc
# exports it; RFunctions also defines a same-arity version but does not export
# it. Bind the Rsc one explicitly so source order cannot decide which wins.
calculate_proportion <- Rsc::calculate_proportion

# options --------

# Seurat/sctransform ship >500 MiB of globals into future workers, exceeding the
# default future.globals.maxSize. The sctransform chapters cannot run without this.
options(future.globals.maxSize = 16 * 1024^3)

# paths --------

data_dir <- here("data")
results_dir <- here("results")
cache_dir <- here("results", "my_cache")


# ggtheme helpers --------

turn_xaxis_labels <- function() {
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}

# other functions ------

calculate_silhouette_score <- function(data, num_clusters) {
  hclust_result <- hclust(dist(data), method = "complete")
  clusters <- cutree(hclust_result, k = num_clusters)
  silhouette_score <- silhouette(clusters, dist(data))
  return(mean(silhouette_score[, "sil_width"]))
}