i <- 3
j <- 0.3
#this vectors will hold the cluster's number ot the cell type
tumer_clusters <- c()
Bcell_clusters <- c()
Tcell_clusters <- c()
macro_clusters <- c()
eryth_clusters <- c()
dendr_clusters <- c()
cd24_clusters <- c()

#the loop goes over the clusters and looking for the cell type
for(numCluster in unique(GenesInCluster_10$cluster)){
  
    #optimal score of this iteration
    geneInClusterScore <- c(i,j)
    
    tumer_cells <- GenesInCluster_10 %>% filter(cluster == numCluster & (gene == "KRT19" | gene == "KRT18")) %>%
      group_by(cluster) %>% summarise(min_p_val_adj = min(p_val_adj))
    if (nrow(tumer_cells) > 0) {
      p_val_adj_tumer <- tumer_cells$min_p_val_adj 
      tumer_clusters <- c(tumer_clusters, list(numCluster, p_val_adj_tumer))
    }
    str(tumer_clusters)
    
    
    Bcell_cells <- GenesInCluster_10 %>% filter(cluster == numCluster & (gene == "IGHG1" | gene == "IGHG4" | gene == "IGKC")) %>%
      group_by(cluster) %>% summarise(min_p_val_adj = min(p_val_adj))
    if (nrow(Bcell_cells) > 0) {
      p_val_adj_Bcell <- Bcell_cells$min_p_val_adj 
      Bcell_clusters <- c(Bcell_clusters, list(numCluster, p_val_adj_Bcell))
    }
    str(Bcell_clusters)
    
    Tcell_cells <- GenesInCluster_10 %>% filter(cluster == numCluster & (gene == "CD8A")) %>%
      group_by(cluster) %>% summarise(min_p_val_adj = min(p_val_adj))
    if (nrow(Tcell_cells) > 0) {
      p_val_adj_Tcell <- Tcell_cells$min_p_val_adj 
      Tcell_clusters <- c(Tcell_clusters, list(numCluster, p_val_adj_Tcell))
    }
    str(Tcell_clusters)
    
    macro_cells <- GenesInCluster_10 %>% filter(cluster == numCluster & (gene == "HIF1A" | gene == "LGMN"))%>%
      group_by(cluster) %>% summarise(min_p_val_adj = min(p_val_adj))
    if (nrow(macro_cells) > 0) {
      p_val_adj_macro <- macro_cells$min_p_val_adj 
      macro_clusters <- c(macro_clusters, list(numCluster, p_val_adj_macro))
    }
    str(macro_clusters)
    
    
    eryth_cells <- GenesInCluster_10 %>% filter(cluster == numCluster & (gene == "CENPF")) %>%
      group_by(cluster) %>% summarise(min_p_val_adj = min(p_val_adj))   
    if (nrow(eryth_cells) > 0) {
      p_val_adj_eryth <- eryth_cells$min_p_val_adj
      eryth_clusters <- c(eryth_clusters, list(numCluster,p_val_adj_eryth))
    }
    str(eryth_clusters)
    
    dendr_cells <- GenesInCluster_10 %>% filter(cluster == numCluster & (gene == "TCF4"))%>%
      group_by(cluster) %>% summarise(min_p_val_adj = min(p_val_adj))  
    if (nrow(dendr_cells) > 0) {
      p_val_adj_dendr <-dendr_cells$min_p_val_adj
      dendr_clusters <- c(dendr_clusters, list(numCluster, p_val_adj_dendr))
    }
    str(dendr_clusters)
    
    
    cd24_cells <- GenesInCluster_10 %>% filter(cluster == numCluster & (gene == "CD24"))%>%
      group_by(cluster) %>% summarise(min_p_val_adj = min(p_val_adj))  
    if (nrow(cd24_cells) > 0) {
      p_val_adj_cd24 <- cd24_cells$min_p_val_adj
      cd24_clusters <- c(cd24_clusters, list(numCluster, p_val_adj_cd24))
    }
    str(cd24_clusters)
    
    
    if (is.null(tumer_clusters) | is.null(Bcell_clusters) | is.null(Tcell_clusters) | is.null(macro_clusters)) {
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
    cat("geneInClusterScore", geneInClusterScore, "\n")

    browser()
  }
