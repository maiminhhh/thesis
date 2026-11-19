# Optional extra R package libraries: set R_LIBS_PROJECT to one or more folders (separated by ':') to put them
# first on the library path. (On the university server this pointed at a personal and two shared package folders.)
extra_libs <- strsplit(Sys.getenv('R_LIBS_PROJECT'), ':', fixed = TRUE)[[1]]
if (length(extra_libs)) .libPaths(unique(c(extra_libs, .libPaths())))

# attach only packages whose functions are called without a pkg:: prefix somewhere
# in the projects; others (batchelor, zellkonverter, SummarizedExperiment, ...) are
# used as pkg::fun and load on demand. Order matters for masking: rliger after
# Seurat, plotly last.
library(Matrix)
library(Seurat)
library(tidyverse)      # dplyr, tidyr, ggplot2, readr, stringr, purrr, ...
library(jsonlite)
library(anndataR)
library(arrow)
library(rhdf5)
library(SingleR)
library(scCustomize)
library(data.table)
library(reticulate)
library(harmony)
library(rliger)
library(mclust)
library(sf)
library(BiocParallel)
library(pheatmap)
library(plotly)
