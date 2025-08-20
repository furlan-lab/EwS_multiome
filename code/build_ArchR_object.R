setwd("/fh/fast/furlan_s/user/owalt/ewings/multiome/all_parents")

suppressPackageStartupMessages({
  library(monocle3)
  library(m3addon)
  library(reticulate)
  library(viridis)
  library(openxlsx)  
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(xfun)
  library(pals)
  library(RColorBrewer)
  library(Signac)
  library(Seurat)
  library(ggplot2)
  library(future)
  library(GenomicRanges)
  library(EnsDb.Mmusculus.v79)
  library(ComplexHeatmap)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(JASPAR2020)
  library(TFBSTools)
  library(patchwork)
  library(stringr)
  library(ggrepel)
  library(ggsignif)
  library(ggpubr)
  library(ArchR)
  library(harmony)
  library(viridis)
  library(scCustomize)
  library(SeuratDisk)
  library(parallel)
})

set.seed(1234)
dyn.load('/app/software/ArrayFire/3.8.1/lib64/libaf.so.3')
library(RcppArrayFire)

m1<- loadArchRProject()

Sys.setenv(MODULEPATH="/app/modules/all", MODULEPATH_ROOT="/app/modules/all", MODULESHOME="/app/lmod/lmod")
source("/app/lmod/lmod/init/R", echo=FALSE, max=Inf)
module("load", "MACS2/2.2.6-foss-2019b-Python-3.7.4")

#if you're running on rhino make sure to ml MACS2 before launching R, if you're running locally you need to install macs2 with python..definitely less messy to run on rhino
pathToMacs2 <- findMacs2()

m1 <- addReproduciblePeakSet(
  ArchRProj = m1, 
  groupBy = "cell_line", 
  pathToMacs2 =pathToMacs2
)

m1<- addPeakMatrix(m1)

if(is.null(getPeakSet(m1))){
  peak_se<-getMatrixFromProject(m1, useMatrix = "PeakMatrix")
  m1<-addPeakSet(m1, peak_se@rowRanges)
}

getAvailableMatrices(m1)
getPeakSet(m1)

m1 <- addPeak2GeneLinks(
  ArchRProj = m1,
  reducedDims = "ATAC_LSI",
  useMatrix = "GeneExpressionMatrix"
)

m1 <- addBgdPeaks(m1)
m1 <- addMotifAnnotations(ArchRProj = m1, force = T, motifSet = "JASPAR2020")
m1 <- addDeviationsMatrix(
  ArchRProj = m1, 
  peakAnnotation = "Motif",
  force = TRUE
)

saveArchRProject(m1)
