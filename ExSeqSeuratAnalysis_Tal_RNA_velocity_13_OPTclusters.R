rm(list = ls())
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyverse)

# the following is according to:
# https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html


counts_raw <- read.table("conbined_counts_0_296Genes.csv", header=TRUE,row.names = 1, sep = ",") 
counts_raw<-as.data.frame(t(counts_raw))
counts_norm <-counts_raw
## normalization for RNA velocity = spliced count for a given gene per cell / total spliced counts for all genes in that given cell * mean total spliced counts per cell
counts_norm <- counts_raw/colSums(counts_raw)*(mean(colSums(counts_raw))) 
#counts_norm <- log2(counts_raw+1) 

## convert the gene expression data per cell into Seurat object
MBC_ExSeq <- CreateSeuratObject(counts = counts_norm, project = "BreastCancer", assay = "RNA", min.cells = 1, min.features = 1) # note that only cells with >1 features (expressed genes) are included -- this doesn't matter in practice
MBC_ExSeq <- AddMetaData(object = MBC_ExSeq, metadata = colSums(counts_raw), col.name = "counts_per_cell")

counts_dim <- dim(counts_raw)

numGenes <- counts_dim[1]
numCells <- counts_dim[2]


###### Choose 100 highest variable genes and plot top 10 on variance vs expression plot
MBC_ExSeq_QC <- FindVariableFeatures(MBC_ExSeq, selection.method = "vst", nfeatures = 100) 

## scale the expression of each gene so that the mean expression is 0 and the variance is 1 across all cells
all.genes <- rownames(MBC_ExSeq_QC) 
MBC_ExSeq_QC <- ScaleData(MBC_ExSeq_QC, features = all.genes)

##### Pick Gene set to use for PCA
# run PCA analysis 
# NOTE (1): this is more than just PCA - this step constructs K-nearest neighbor graph based on the euclidean distance in PCA space, then uses Louvain algorithm to find clusters in the PCA space
# NOTE (2): this is supervised as we choose the important genes for the PCA using prior knowledge and not according to variance
# the list of genes that we use for the PCA is according to this paper: 
# "A Structured Tumor-Immune Microenvironment in Triple Negative Breast Cancer Revealed by Multiplexed Ion Beam Imaging"
#MBC <- RunPCA(MBC, features = VariableFeatures(object = MBC))

## Geneset 4: extending Shahar's genelist by adding: "IGHG1", "IGHG4", "IGKC", "IGHM"
GenesTouse_ExtendedShahar <- c("CD3G" ,"CD68","FOXP3", "CD4" ,"CD8A" ,"CD3D" ,"CD3E","HLA-DRA") #non-tumor marker genes 
GenesTouse_ExtendedShahar <- c(GenesTouse_ExtendedShahar,"EGFR","GRB7","ERBB2","PGR","CD44","CD24","ALDH1A3","EPCAM","KRT19","KRT18","CDH1") #tumor marker genes
GenesTouse_ExtendedShahar <- c(GenesTouse_ExtendedShahar,"IGHG1", "IGHG4", "IGKC", "IGHM") #b-cells
#GenesTouse_ExtendedShahar <- c(GenesTouse_ExtendedShahar,"COL1A1","COL1A2","MMP11","CDH11","HSPG2") #fibroblast
GenesTouse_ExtendedShahar <- c(GenesTouse_ExtendedShahar,"SULF1") #fibroblast

##### Run PCA, UMAP

# THIS IS WHERE YOU PICK THE GENESET TO USE
#MBC_ExSeq_QC <- RunPCA(MBC_ExSeq_QC, features = GenesTouse_ExtendedShahar) 
MBC_ExSeq_QC <- RunPCA(MBC_ExSeq_QC) 

# display information about the 'dimensionality' in the PCA space
###




###  
  print(MBC_ExSeq_QC[["pca"]], dims = 1:5, nfeatures = 5)
  DimHeatmap(MBC_ExSeq_QC, dims = 1:15, cells = 500, balanced = TRUE)
  ElbowPlot(MBC_ExSeq_QC)
  
  ###
  MBC_ExSeq_QC<- JackStraw(object = MBC_ExSeq_QC, dims = 15)
  MBC_ExSeq_QC <- ScoreJackStraw(MBC_ExSeq_QC, dims = 1:15)
  JackStrawPlot(object = MBC_ExSeq_QC, dims = 1:15) + theme(legend.position="bottom")
  ####
  score_df <- data.frame()
  for (i in 3:7){

    print(i)
  # Cluster
    MBC_ExSeq_QC <- FindNeighbors(MBC_ExSeq_QC, dims = 1:i)
    ###
    for (j in seq(0.3,1, by = 0.1)) {
      print(j)
    ###
      MBC_ExSeq_QC <- FindClusters(MBC_ExSeq_QC, resolution = j) # was 1 for nCount_RNA > 100 
      
      # Dimensionality reduction 
      MBC_ExSeq_QC <- RunUMAP(MBC_ExSeq_QC, dims = 1:i)
      the_plot <- DimPlot(MBC_ExSeq_QC, reduction = "umap")
      the_plot
      
      # detect marker genes for each cluster
      MBC_ExSeq_QC.markers <- FindAllMarkers(MBC_ExSeq_QC, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
      GenesInCluster <- MBC_ExSeq_QC.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
      
      GenesInCluster_10 <- MBC_ExSeq_QC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
      
      ###
      

      #this vectors will hold the cluster's number ot the cell type
      tumer_clusters <- c()
      Bcell_clusters <- c()
      Tcell_clusters <- c()
      macro_clusters <- c()
      eryth_clusters <- c()
      dendr_clusters <- c()
      cd24_clusters <- c()
      fibroblast_clusters <- c()
      
      #the loop goes over the clusters and looking for the cell type
      for(numCluster in unique(GenesInCluster_10$cluster)){
        
        #optimal score of this iteration
        geneInClusterScore <- c(i,j)
        
        tumer_cells <- GenesInCluster_10 %>% filter(cluster == numCluster & (gene == "KRT19" | gene == "KRT18"))
        if (nrow(tumer_cells) > 0) {
          tumer_clusters <- c(tumer_clusters, numCluster)
        }
        cat("tumer_clusters", tumer_clusters, "\n")
        
        Bcell_cells <- GenesInCluster_10 %>% filter(cluster == numCluster & (gene == "IGHG1" | gene == "IGHG4" | gene == "IGKC"))
        if (nrow(Bcell_cells) > 0) {
          Bcell_clusters <- c(Bcell_clusters, numCluster)
        }
        cat("Bcell_clusters",Bcell_clusters, "\n")
        
        Tcell_cells <- GenesInCluster_10 %>% filter(cluster == numCluster & (gene == "CD3G" | gene == "CD3D" | gene == "CD3E"))
        if (nrow(Tcell_cells) > 0) {
          Tcell_clusters <- c(Tcell_clusters, numCluster)
        }
        cat("Tcell_clusters", Tcell_clusters, "\n")
        
        macro_cells <- GenesInCluster_10 %>% filter(cluster == numCluster & (gene == "HIF1A" | gene == "LGMN"))
        if (nrow(macro_cells) > 0) {
          macro_clusters <- c(macro_clusters, numCluster)
        }
        cat("macro_clusters", macro_clusters, "\n")
        
        
        eryth_cells <- GenesInCluster_10 %>% filter(cluster == numCluster & (gene == "CENPF"))
        if (nrow(eryth_cells) > 0) {
          eryth_clusters <- c(eryth_clusters, numCluster)
        }
        cat("eryth_clusters", eryth_clusters, "\n")
        
        dendr_cells <- GenesInCluster_10 %>% filter(cluster == numCluster & (gene == "TCF4"))
        if (nrow(dendr_cells) > 0) {
          dendr_clusters <- c(dendr_clusters, numCluster)
        }
        cat("dendr_clusters", dendr_clusters, "\n")
        
        
        cd24_cells <- GenesInCluster_10 %>% filter(cluster == numCluster & (gene == "CD24"))
        if (nrow(cd24_cells) > 0) {
          cd24_clusters <- c(cd24_clusters, numCluster)
        }
        cat("cd24_clusters", cd24_clusters, "\n")
        
        fibroblast_cells <- GenesInCluster_10 %>% filter(cluster == numCluster & (gene == "COL3A1"))
        if (nrow(fibroblast_cells) > 0) {
          fibroblast_clusters <- c(fibroblast_clusters, numCluster)
        }
        cat("fibroblast_clusters", fibroblast_clusters, "\n")
        
        #if (is.null(tumer_clusters) | is.null(Bcell_clusters) | is.null(Tcell_clusters) | is.null(macro_clusters)) {
        if (is.null(tumer_clusters) | is.null(Bcell_clusters) | is.null(macro_clusters)) {
          next
        } else geneInClusterScore <- c(geneInClusterScore, 1)
        
        if (!is.null(eryth_clusters)  && !is.null(dendr_clusters)) {
          geneInClusterScore <- c(geneInClusterScore, 3)
        } else if (!is.null(eryth_clusters) | !is.null(dendr_clusters)) {
          geneInClusterScore <- c(geneInClusterScore, 2)
        }
        if (!is.null(cd24_clusters)) {
          geneInClusterScore <- c(geneInClusterScore, 4)
        }
        if (is.null(Tcell_clusters)){
          next
        } else geneInClusterScore <-c(geneInClusterScore, 5)
        
        if (is.null(fibroblast_clusters)){
          next
        } else geneInClusterScore <-c(geneInClusterScore, 6)
        
        cat("geneInClusterScore", geneInClusterScore, "\n")
        
        #browser()
      }
      if(length(geneInClusterScore) < 7){
        geneInClusterScore <- c(geneInClusterScore, rep(NA, 7- length(geneInClusterScore)))
      }
      score_df <-rbind(score_df, geneInClusterScore)
      #browser()
      colnames(score_df) <- c("dim", "res", "tumor, Bcell, macro", "eryth&dendr","cd24", "Tcell", "fibro")
    }
  }


###########################
# P-value preffernce

## in this part i want to check if there are clusters that are common to 2 genes or more and if so to check which p-value is better, to do that im going to make a data fraim out of all the vectors
#max_length <- max(length(tumer_clusters), length(Bcell_clusters), length(Tcell_clusters), length(macro_clusters), length(eryth_clusters), length(dendr_clusters), length(cd24_clusters))
#clusters_df <- data.frame(
#  Tumer_Clusters = c(tumer_clusters, rep(NA, max_length - length(tumer_clusters))),
#  Bcell_Clusters = c(Bcell_clusters, rep(NA, max_length - length(Bcell_clusters))),
#  Tcell_Clusters = c(Tcell_clusters, rep(NA, max_length - length(Tcell_clusters))),
#  Macro_Clusters = c(macro_clusters, rep(NA, max_length - length(macro_clusters))),
#  Eryth_Clusters = c(eryth_clusters, rep(NA, max_length - length(eryth_clusters))),
#  Dendr_Clusters = c(dendr_clusters, rep(NA, max_length - length(dendr_clusters))),
#  CD24_Clusters = c(cd24_clusters, rep(NA, max_length - length(cd24_clusters)))
#)

#duplicates <- find_duplicates(clusters_df)
#print(duplicates)

###


#finds duplicates gene in 2 or more clusters

#find_duplicates <- function(df) {
#  df_long <- stack(df)
#  df_long <- df_long[!is.na(df_long$values), ]
  
#  table_values <- table(df_long$values)
#  duplicate_values <- names(table_values[table_values > 1])
  
#  result <- sapply(duplicate_values, function(value) {
#    columns <- unique(df_long$ind[df_long$values == as.numeric(value)])
#    data.frame(index = as.numeric(value), columns = paste(columns, collapse = ", "))
#  })
  
#  result_df <- do.call(rbind, result)
#  rownames(result_df) <- NULL
  
#  return(result_df)
#}

###
write.csv(GenesInCluster_10,"GenesInCluster_10_res0.4_5dim.csv")
###

