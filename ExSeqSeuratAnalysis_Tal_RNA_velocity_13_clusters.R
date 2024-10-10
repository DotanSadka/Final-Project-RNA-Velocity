rm(list = ls())
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyverse)

# the following is according to:
# https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html


counts_raw <- read.table("conbined_counts_0_271Genes.csv", header=TRUE,row.names = 1, sep = ",") 
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
print(MBC_ExSeq_QC[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(MBC_ExSeq_QC, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(MBC_ExSeq_QC)

###
MBC_ExSeq_QC<- JackStraw(object = MBC_ExSeq_QC, dims = 15)
MBC_ExSeq_QC <- ScoreJackStraw(MBC_ExSeq_QC, dims = 1:15)
JackStrawPlot(object = MBC_ExSeq_QC, dims = 1:15) + theme(legend.position="bottom")
####

# Cluster
MBC_ExSeq_QC <- FindNeighbors(MBC_ExSeq_QC, dims = 1:4)
MBC_ExSeq_QC <- FindClusters(MBC_ExSeq_QC, resolution = 0.7) # was 1 for nCount_RNA > 100 

# Dimensionality reduction 
MBC_ExSeq_QC <- RunUMAP(MBC_ExSeq_QC, dims = 1:4)
the_plot <- DimPlot(MBC_ExSeq_QC, reduction = "umap")
the_plot

# detect marker genes for each cluster
MBC_ExSeq_QC.markers <- FindAllMarkers(MBC_ExSeq_QC, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
GenesInCluster <- MBC_ExSeq_QC.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

GenesInCluster_10 <- MBC_ExSeq_QC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

###
write.csv(GenesInCluster_10,"GenesInCluster_10_res0.7_4dim.csv")
###

FeaturePlot(MBC_ExSeq_QC, features = c("IGHG1","IGHG4","B-cell_IGKC")) # B-cells
FeaturePlot(MBC_ExSeq_QC, features = c("CD8A","CD3D", "CD4")) # T-cells
FeaturePlot(MBC_ExSeq_QC, features = c("HLA-DRA")) # Macrophage
FeaturePlot(MBC_ExSeq_QC, features = c("HSPG2")) # fibroblast
FeaturePlot(MBC_ExSeq_QC, features = c("PGR","EGFR","ALDH1A3","CD44")) # tumor


## 1311 run only: FindVariableFeatures, ScaleData,RunPCA, FindNeighbors, FindClusters, RunUMAP, FindAllMarkers...
## 13 clusters (0.6 resolution)
#new.cluster.ids <- c("Fibroblast_SULF1", "Un", "B-cell_IGHG4", "Tumor_EPCAM", "Tumor_CD24", "Tumor_PGR", "Tumor_CD44", "T-cell_CD4", "T-cell_CD3E", "T-cell_CD8A", "T-cell_CD3D", "Tumor_EGFR",  "Tumor_ALDH1A3") #13 0.6

new.cluster.ids <- c("UnKnown1", "Macrophage_HIF1A" , "Unknown2", "Unknown3", "Tumor_KRT19", "T-cell_CD8A", "Unknown4", "Macrophage_LGMN", "Tumor_CD24", "B-cell_IGHG4", "T-cell_CD8A", "Fibroblast_COL3A1")

names(new.cluster.ids) <- levels(MBC_ExSeq_QC)
MBC_ExSeq_QC <- RenameIdents(MBC_ExSeq_QC, new.cluster.ids)

## plot UMAP - 13 clusters 
DimPlot(MBC_ExSeq_QC, reduction = "umap", label = TRUE, pt.size = 1) 


## MBC_ExSeq_QC - cell types
labels_c = MBC_ExSeq_QC@active.ident
labels_c_add <- sub("_.*", "",labels_c)
labels_c_add <- data.frame(labels_c,labels_c_add)


## PCA visualization
cell_pc<-MBC_ExSeq_QC@reductions$pca@cell.embeddings
cell_pc_df<- data.frame(pc1=cell_pc[,'PC_1'], pc2=cell_pc[,'PC_2'], pc3=cell_pc[,'PC_3'])
cell_pc_df <- merge(cell_pc_df, labels_c_add ,by=0)
library(plotly)
plot_ly(x=cell_pc_df$pc1, y=cell_pc_df$pc2, z=cell_pc_df$pc3, type="scatter3d", mode="markers", color=cell_pc_df$labels_c)

## 
library(dplyr)
cell_pc_df_without_Un <- cell_pc_df %>% filter(labels_c_add=='T-cell')
#colors <- c("red", "green", "blue", "purple")
#colors <- colors[factor(cell_pc_df_without_Un$labels_c)]
#plot_ly(x=cell_pc_df_without_Un$pc1, y=cell_pc_df_without_Un$pc2, z=cell_pc_df_without_Un$pc3, type="scatter3d", mode="markers", color=colors)
plot_ly(x=cell_pc_df_without_Un$pc1, y=cell_pc_df_without_Un$pc2, z=cell_pc_df_without_Un$pc3, type="scatter3d", mode="markers", color=factor(cell_pc_df_without_Un$labels_c), marker=list(line= list(color="black", width=2)))

## add column that identifies only T-cells
cell_pc_df[,'check_T']<- 'other'
cell_pc_df[cell_pc_df['labels_c']=='T-cell_CD8A', 'check_T'] <- 'T-cell_CD8A'
cell_pc_df[cell_pc_df['labels_c']=='T-cell_CD3E', 'check_T'] <- 'T-cell_CD3E'
cell_pc_df[cell_pc_df['labels_c']=='T-cell_CD3D', 'check_T'] <- 'T-cell_CD3D'
cell_pc_df[cell_pc_df['labels_c']=='T-cell_CD3G', 'check_T'] <- 'T-cell_CD3G'
cell_pc_df[cell_pc_df['labels_c']=='T-cell_CD4', 'check_T'] <- 'T-cell_CD4'
plot_ly(x=cell_pc_df$pc1, y=cell_pc_df$pc2, z=cell_pc_df$pc3, type="scatter3d", mode="markers", color=cell_pc_df$check_T)

## save to csv file
write.csv(cell_pc_df,"cell_pc_df.csv", row.names = FALSE)

unique(cell_pc_df$labels_c)

#####
#sorting "cell_pc_df" file by "lable_c" column (cell types) to check hoe many neighbors to take

file_path <- "cell_pc_df.csv"
df <- read.csv(file_path)

df_sorted <- df[order(df$labels_c), ]

labels_count <- table(df_sorted$labels_c)

print(labels_count)

####


#############################
#for (t in seq(0.2, 5.0, by = 0.2)) {
#for (t in seq(1)) {
  ## Plot S(t) from RNA velocity analysis on the same 3D plot of PCA 
  pc_original <- read.table("cell_pc_df.csv", header=TRUE, sep = ",", row.names=1)
  #input_filename <- sprintf("Hyperparameters/t/t=%.1f_Gene_expression_filtered.csv", t)
  input_filename <- sprintf("t=1.0_Gene_expression_filtered_1.csv")
  st_origin <- read.table(input_filename , header=TRUE, sep = ",", row.names=1, , skip=1)
  st_origin <- as.data.frame(t(st_origin))
  
  ## normalization for RNA velocity = spliced count for a given gene per cell / total spliced counts for all genes in that given cell * mean total spliced counts per cell
  ## t1_Gene_expression_origin already divided by total spliced counts for all genes in that given cell, but we have to multiply the results by the mean total spliced counts per cell
  st <- st_origin*(mean(colSums(st_origin))) 
  
  
  ## convert the gene expression data per cell into Seurat object
  MBC_ExSeq_velocity <- CreateSeuratObject(counts = st, project = "BreastCancer", assay = "RNA", min.cells = 1, min.features = 1) # note that only cells with >1 features (expressed genes) are included -- this doesn't matter in practice
  MBC_ExSeq_velocity <- AddMetaData(object = MBC_ExSeq_velocity, metadata = colSums(st), col.name = "counts_per_cell")
  
  counts_dim_velocity <- dim(st)
  
  numGenes_velocity <- counts_dim_velocity[1]
  numCells_velocity <- counts_dim_velocity[2]
  
  
  ###### Choose 100 highest variable genes and plot top 10 on variance vs expression plot
  MBC_ExSeq_velocity_QC <- FindVariableFeatures(MBC_ExSeq_velocity, selection.method = "vst", nfeatures = 100) 
  
  ## scale the expression of each gene so that the mean expression is 0 and the variance is 1 across all cells
  all.genes <- rownames(MBC_ExSeq_velocity_QC) 
  MBC_ExSeq_velocity_QC <- ScaleData(MBC_ExSeq_velocity_QC, features = all.genes)
  
  ##### Pick Gene set to use for PCA
  # run PCA analysis 
  # NOTE (1): this is more than just PCA - this step constructs K-nearest neighbor graph based on the euclidean distance in PCA space, then uses Louvain algorithm to find clusters in the PCA space
  # NOTE (2): this is supervised as we choose the important genes for the PCA using prior knowledge and not according to variance
  # the list of genes that we use for the PCA is according to this paper: 
  # "A Structured Tumor-Immune Microenvironment in Triple Negative Breast Cancer Revealed by Multiplexed Ion Beam Imaging"
  #MBC <- RunPCA(MBC, features = VariableFeatures(object = MBC))
  
  ## Geneset 4: extending Shahar's genelist by adding: "IGHG1", "IGHG4", "IGKC", "IGHM"
  #GenesTouse_ExtendedShahar <- c("CD3G" ,"CD68","FOXP3", "CD4" ,"CD8A" ,"CD3D" ,"CD3E","HLA-DRA") #non-tumor marker genes 
  #GenesTouse_ExtendedShahar <- c(GenesTouse_ExtendedShahar,"EGFR","GRB7","ERBB2","PGR","CD44","CD24","ALDH1A3","EPCAM","KRT19","KRT18","CDH1") #tumor marker genes
  #GenesTouse_ExtendedShahar <- c(GenesTouse_ExtendedShahar,"IGHG1", "IGHG4", "IGKC", "IGHM") #b-cells
  #GenesTouse_ExtendedShahar <- c(GenesTouse_ExtendedShahar,"SULF1") #fibroblast
  
  ##### Run PCA, UMAP
  
  # THIS IS WHERE YOU PICK THE GENESET TO USE
  MBC_ExSeq_velocity_QC <- RunPCA(MBC_ExSeq_velocity_QC) 
  
  # display information about the 'dimensionality' in the PCA space 
  #print(MBC_ExSeq_velocity_QC[["pca"]], dims = 1:5, nfeatures = 5)
  #DimHeatmap(MBC_ExSeq_velocity_QC, dims = 1:15, cells = 500, balanced = TRUE)
  #ElbowPlot(MBC_ExSeq_velocity_QC)
  
  # get all PCs values per cell
  s_all_PCs <-MBC_ExSeq_velocity_QC@reductions$pca@cell.embeddings
  write.csv(s_all_PCs,"s_all_PCs.csv")
  
  ## PCA to df
  cell_pc_df_st<- data.frame(pc1=s_all_PCs[,'PC_1'], pc2=s_all_PCs[,'PC_2'], pc3=s_all_PCs[,'PC_3'])
  ## save to csv file
  #output_filename <- sprintf("Hyperparameters/t/st_filter/st_t=%.1f_filter.csv", t)
  output_filename <- sprintf("st_t=1.0_filter.csv")
  write.csv(cell_pc_df_st, output_filename)
#}




################## Do not run the sections below - not working

## PCA visualization
st_t3_filter <- read.table("st_t3_filter.csv", header=TRUE, sep = ",", row.names=1)

library(plotly)

# Create the scatter plot
fig <- plot_ly() %>%
  add_markers(
    data = pc_original,
    x = ~pc1,
    y = ~pc2,
    z = ~pc3,
    color = I("blue"),
    size = 3,
    name = "Data Points"
  )
fig
# Add arrows as lines to the plot
for (i in 1:nrow(pc_original)) {
  fig <- fig %>%
    add_lines(
      x = c(pc_original$pc1[i], st_t3_filter$pc1[i]),
      y = c(pc_original$pc2[i], st_t3_filter$pc2[i]),
      z = c(pc_original$pc3[i], st_t3_filter$pc3[i]),
      color = I("black"),
      name = "",
      showlegend = FALSE
    )# %>%
  #add_markers(
  #  x = st_t3_filter$pc1[i],
  #  y = st_t3_filter$pc2[i],
  #  z = st_t3_filter$pc3[i],
  #  color = I("black"),
  #  size = I(3),
  #  name = "",
  #  showlegend = FALSE
  #)
}

# Set layout options
fig <- fig %>%
  layout(
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    )
  )

# Show the plot
fig




##################

plot_ly(x=pc_original$pc1, y=pc_original$pc2, z=pc_original$pc3, type="scatter3d", mode="markers", color=pc_original$labels_c)
pc_original[,'check_T_CD8A']<- 'other'
pc_original[pc_original['labels_c']=='T-cell_CD8A', 'check_T_CD8A'] <- 'T-cell_CD8A'
pc_original_8A <- pc_original %>% filter(labels_c=='T-cell_CD8A')

plot_ly(x=pc_original_8A$pc1, y=pc_original_8A$pc2, z=pc_original_8A$pc3, type="scatter3d", mode="markers", color=pc_original_8A$labels_c)


##################################

# Load the required libraries
library(Seurat)
library(dplyr)

#
pc_original <- read.table("cell_pc_df.csv", header=TRUE, sep = ",", row.names=1)
colnames(pc_original)[1] <- "cell_id"
st_t1_nofilter <- read.table("st_t1_nofilter.csv", header=TRUE, sep = ",")
colnames(st_t1_nofilter)[1] <- "cell_id"
s_all_PCs <- read.table("s_all_PCs.csv", header=TRUE, sep = ",", row.names=1)


# Get the index order based on 'cell_id' column
index_order <- match(pc_original$cell_id, st_t1_nofilter$cell_id)

# Reorder the st_t1_nofilter dataframe
st_t1_nofilter_reordered <- st_t1_nofilter[index_order, ]

# Join the two dataframes based on 'cell_id' column
merged_df <- merge(st_t1_nofilter_reordered, pc_original[, c("cell_id", "labels_c")], by="cell_id")

#filter T-cell_CD8A
filtered_df <- subset(merged_df, labels_c == 'T-cell_CD8A')
filtered_row_names <- rownames(filtered_df)

# Subset rows in s_all_PCs

common_rows <- rownames(s_all_PCs) %in% rownames(filtered_df)
s_all_PCs_subset <- s_all_PCs[common_rows, ]
                      
# Create a Seurat object with your PCA data
count_matrix <- matrix(0, nrow = nrow(filtered_df), ncol = 1)  # Replace 1 with the appropriate number of genes
seurat_obj <- CreateSeuratObject(counts = filtered_df[, c("pc1", "pc2", "pc3")])

seurat_obj[["pca"]] <- pc_original

# Assign the 'labels_c' column from pc_original to the Seurat object
seurat_obj[["labels_c"]] <- pc_original[["labels_c"]]

# Compute the k-param nearest neighbors
seurat_obj <- FindNeighbors(
  object = seurat_obj,
  dims = 1:3,  # Specify the dimensions to use for neighbor computation
  k.param = 10  # Specify the number of nearest neighbors to consider
)

# Apply Louvain algorithm for clustering
seurat_obj <- FindClusters(
  object = seurat_obj,
  resolution = 0.5  # Specify the resolution parameter for modularity optimization
)

# Get the cluster labels
cluster_labels <- seurat_obj$seurat_clusters



