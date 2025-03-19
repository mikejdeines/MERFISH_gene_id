get_mutual_info_data <- function(seuratobject, de_genes){
  require(Seurat)
  #' Saves DE gene expression data from a Seurat object as a .csv file.
  #' This file can be used as input to python code to calculate mutual information.
  #' @param seuratobject an SCTransformed, integrated Seurat object.
  #' @param de_genes a data.frame of DE genes, generated using the seurat_de_genes function.
  #' @returns a path to the .csv file containing the SCT-corrected counts for all unique DE genes in the input Seurat object.
  unique_de_genes <- unique(de_genes$Genes)
  SCT <- seuratobject[["SCT"]]$counts
  SCT <- SCT[unique_de_genes,]
  out.path <- paste0(seuratobject, "_SCT.csv")
  write.csv(SCT, file = out.path)
  return(out.path)
}
de_combination_genes <- function(de_genes){
  #' Finds the genes with the highest DE score for all pairs of clusters.
  #' If multiple genes are DE between the same pair of clusters, the gene with the highest DE score is retained.
  #' @param de_genes a data.frame of DE genes, generated using the seurat_de_genes function.
  #' @returns a vector of gene names.
  require(tidyr)
  require(dplyr)
  cluster_combos <- unique(de_genes$combination)
  combo_genes <- data.frame()
  for (i in 1:length(cluster_combos)){
    genes <- de_genes[de_genes$combination == cluster_combos[i],]
    genes <- genes[1:10,]
    combo_genes <- rbind(combo_genes, genes)
  }
  combo_genes <- tidyr::drop_na(combo_genes)
  gene_list <- c()
  cluster_list <- unique(de_genes$Cluster_1)
  for (i in 1:length(cluster_list)){
    df <- combo_genes[combo_genes$Cluster_1 == cluster_list[i],]
    df <- df %>% add_count(Genes, name = "Gene_count")
    df <- df %>% group_by(Gene_count, Cluster_2) %>% slice_max(order_by = DE_score, n = 1) %>% ungroup()
    genes <- unique(df$Genes)
    gene_list <- append(gene_list, genes)
  }
  gene_list <- unique(gene_list)
  return(gene_list)
}
assemble_genes <- function(mutualinfopath, manualpath, gene_list, mutual_info_percentile = 0.95){
  #' Assembles a gene list from mutual information output, manual gene selections, and DE genes.
  #' @param mutualinfopath path to output from the mutual information python function.
  #' @param manualpath path to a .xlsx file containing the names of manually-selected genes.
  #' @param gene_list the vector of gene names output from de_combination_genes.
  #' @param mutual_info_percentile percentile cutoff for mutual information genes. Default 0.95.
  #' @returns a vector of gene names.
  require(readxl)
  mutinfo <- read.csv(mutualinfopath)
  mutual_info_cutoff <- quantile(mutinfo$mutual_information, mutual_info_percentile)
  mutual_info_genes <- mutinfo$symbol[mutinfo$mutual_information > mutual_info_cutoff]
  manual_genes <- read_excel(manualpath)
  merfish_genes <- append(gene_list, mutual_info_genes)
  merfish_genes <- append(merfish_genes, manual_genes$Gene)
  merfish_genes <- unique(merfish_genes)
}
seurat_de_genes <- function(seuratobject, log_threshold = 2, min_pct = 0.5){
  #' Finds DE genes between all pairs of clusters in a Seurat object.
  #' Calculates a DE score equal to -log10(padj). Max score = 300.
  #' Use this function for datasets with < 100 clusters.
  #' @param seuratobject a Seurat object with cluster identities in the Idents slot.
  #' @param log_threshold the minimum log2FC for a gene to be considered DE.
  #' @param min_pct the minimum fraction of gene expression in one cluster for a gene to be considered DE.
  #' @returns a data.frame of differentially-expressed genes between all pairs of clusters.
  require(gtools)
  require(Seurat)
  require(dplyr)
  clusters <- levels(Idents(seuratobject))
  combinations <- permutations(length(clusters), 2, clusters)
  output <- data.frame()
  for (i in 1:length(combinations[,1])){
    de_genes <- Seurat::FindMarkers(seuratobject, ident.1 = combinations[i,1], ident.2 = combinations[i,2], only.pos = TRUE, recorrect_umi = FALSE)
    de_genes$Cluster_1 <- combinations[i,1]
    de_genes$Cluster_2 <- combinations[i,2]
    de_genes$Genes <- rownames(de_genes)
    de_genes <- subset(de_genes, p_val_adj < 1e-5 & avg_log2FC >= log_threshold & pct.1 >= min_pct)
    de_genes$p_val_adj[de_genes$p_val_adj == 0] <- 1e-300
    de_genes$DE_score <- -log10(de_genes$p_val_adj)
    de_genes$combination <- paste(de_genes$Cluster_1, de_genes$Cluster_2, sep = "-")
    output <- rbind(output, de_genes)
  }
  return(output)
}
identify_permutations <- function(seuratobject){
  #' Finds all permutations of 2 clusters in a Seurat object.
  #' The cluster identities should be in the Idents slot of the Seurat object.
  #' @param seuratobject a Seurat object with cluster identities in the Idents slot
  #' @returns an n-by-2 integer of permutations between all pairs of clusters
  require(gtools)
  clusters <- levels(Idents(seuratobject))
  combinations <- permutations(length(clusters), 2, clusters)
  return(combinations)
}
find_optimal_chunk_number <- function(combinations){
  #' Identifies the largest number of chunks the permutations can be divided into.
  #' @param combinations an n-by-2 integer of permutations between all pairs of clusters
  #' @returns the number of chunks the permutations should be divided into
  n <- length(combinations[,1])
  max_value <- ceiling(sqrt(n))
  for (i in 1:max_value){
    if (n %% i == 0){
      optimal_chunk_number <- i
    }
  }
  return(optimal_chunk_number)
}
chunk_permutations <- function(combinations, optimal_chunk_number){
  #' Divides the permutations into a list of chunks.
  #' @param combinations an n-by-2 integer of permutations between all pairs of clusters
  #' @param optimal_chunk_number number of chunks to produce
  #' @returns a list of permutations divided into optimal_chunk_number of chunks
  perm_list <- list()
  chunk_size <- length(combinations[,1])/optimal_chunk_number
  for (i in 1:optimal_chunk_number){
    min <- (i-1)*chunk_size + 1
    max <- i*chunk_size
    perm_chunk <- perm[min : max,]
    perm_list[[paste0("chunk_", i)]] <- perm_chunk
  }
  return(perm_list)
}
identify.de.genes.wilcoxon <- function(combinations, seuratobject, log_threshold, min_pct){
  #' Identifies DE genes using a Wilcoxon rank-sum test.
  #' @param combinations a Matrix of identity combinations
  #' @param seuratobject a Seurat object
  #' @param log_threshold the minimum log2FC value
  #' @param min_pct minimum expression value
  #' @returns a data.frame of the 20 genes with the highest DE score for each combination of identities
  output <- data.frame()
  for (i in 1:length(combinations[,1])){
    de_genes <- Seurat::FindMarkers(seuratobject, ident.1 = combinations[i,1], ident.2 = combinations[i,2], only.pos = TRUE)
    de_genes$Cluster_1 <- combinations[i,1]
    de_genes$Cluster_2 <- combinations[i,2]
    de_genes$Genes <- rownames(de_genes)
    de_genes <- subset(de_genes, p_val_adj < 1e-5 & avg_log2FC >= log_threshold & pct.1 >= min_pct)
    de_genes$p_val_adj[de_genes$p_val_adj == 0] <- 1e-300
    de_genes$DE_score <- -log10(de_genes$p_val_adj)
    de_genes$combination <- paste(de_genes$Cluster_1, de_genes$Cluster_2, sep = "-")
    de_genes <- de_genes[1:20,]
    output <- rbind(output, de_genes)
  }
  return(output)
}
identify.de.genes.deseq2 <- function(combinations, seuratobject, log_threshold, min_pct){
  #' Identifies DE genes using DESeq2.
  #' @param combinations a Matrix of identity combinations
  #' @param seuratobject a Pseudobulked Seurat object
  #' @param log_threshold the minimum log2FC value
  #' @param min_pct minimum expression value
  #' @returns a data.frame of the 20 genes with the highest DE score for each combination of identities
  output <- data.frame()
  for (i in 1:length(combinations[,1])){
    de_genes <- Seurat::FindMarkers(seuratobject, ident.1 = combinations[i,1], ident.2 = combinations[i,2], only.pos = TRUE, test.use = "DESeq2")
    de_genes$Cluster_1 <- combinations[i,1]
    de_genes$Cluster_2 <- combinations[i,2]
    de_genes$Genes <- rownames(de_genes)
    de_genes <- subset(de_genes, p_val_adj < 1e-5 & avg_log2FC >= log_threshold & pct.1 >= min_pct)
    de_genes$p_val_adj[de_genes$p_val_adj == 0] <- 1e-300
    de_genes$DE_score <- -log10(de_genes$p_val_adj)
    de_genes$combination <- paste(de_genes$Cluster_1, de_genes$Cluster_2, sep = "-")
    de_genes <- de_genes[1:20,]
    output <- rbind(output, de_genes)
  }
  return(output)
}
