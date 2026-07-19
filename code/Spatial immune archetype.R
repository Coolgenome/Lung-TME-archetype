### Spatial immune archetypes ###
### load common required packages
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(viridis)
library(ape)
library(ggtree)
library(phangorn)
library(e1071)

### Analyzing spatial distribution patterns of TILs in Visium ###
# load annotation data #
# feature and sample annotations #
D2_Visium_feature_annotation <- readRDS("../data/D2_Visium_feature_annotation.rds")
D2_Visium_sample_annotation <- readRDS("../data/D2_Visium_sample_annotation.rds")

annotation_colors <- list(lesion_type = c("AAH"="#e6194b","AIS"="#3cb44b","MIA"="#ffe119","LUAD"="#4363d8"),
                          patient_ID = c("P1"="#FF0000","P2"="#00FF00","P3"="#0000FF","P4"="#FFFF00",
                                         "P5"="#00FFFF","P6"="#FF00FF","P7"="#FFA500","P8"="#800080",
                                         "P9"="#FFC0CB","P10"="#A52A2A","P11"="#32CD32","P12"="#808000",
                                         "P13"="#800000","P14"="#000080","P15"="#008080","P16"="#40E0D0",
                                         "P17"="#EE82EE","P18"="#FFD700","P19"="#FA8072","P20"="#FF7F50",
                                         "P21"="#4B0082","P22"="#7FFF00","P23"="#DC143C","P24"="#F0E68C",
                                         "P25"="#E6E6FA"),
                          TIL_category = c("Limited_TIL_infiltration"="#33a02c", "T_cell_dominant_infiltration"="#ff7f00", 
                                           "B_cell_dominant_infiltration"="#e31a1c", "TIL_Exclusion"="#6a3d9a"),
                          Mac_category = c("Mac_infiltration" = "#33a02c", "Mac_Exclusion" = "#6a3d9a"),
                          Fib_category = c("Fib_infiltration" = "#33a02c", "Fib_Exclusion" = "#6a3d9a"),
                          region = c("Normal_distal"="#4575b4","Normal_peri"="#74add1","Interface_normal"="#abd9e9",
                                     "Interface_tumor"="#ffffbf","Tumor_region"="#d73027")
                        )

# set up the data frame to record clustering information for each major TME compartment
sample_subcluter_df <- data.frame(sample_ID = rownames(D2_Visium_sample_annotation), 
                                  TIL_node = NA, 
                                  Macrophage_node = NA, 
                                  Fib_node = NA)

## TIL spatial patterns ##
# Z score of bivariant Moran's I against five concentric tissue layers, capped at +/- 0.9
D2_Visium_TIL_moran <- readRDS("../data/D2_Visium_TIL_moran.rds")

# plot spatial pattern heatmap
pdf("Visium_TIL_moran_cluster.pdf",width=15,height=5)
pheatmap::pheatmap(t(D2_Visium_TIL_moran),
                   cluster_rows = T,
                   cluster_cols = T,
                   show_rownames=T,
                   show_colnames=T,
                   clustering_distance_rows = "euclidean", 
                   clustering_distance_cols = "euclidean", 
                   clustering_method = "ward.D2",
                   color = c(colorRampPalette(c('#4575b4','white'))(100),
                               colorRampPalette(c('white','#d73027'))(100)),
                   annotation_col = D2_Visium_sample_annotation,
                   annotation_row = D2_Visium_feature_annotation,
                   annotation_colors = annotation_colors,
                   main = "Visium TIL Moran's I",
                   cellwidth = 12, cellheight = 12, border_color = "white",
                   fontsize = 10, fontsize_row = 8, fontsize_col = 8, angle_col = 45
                )
dev.off()
# note that in the final figures presented, the hierarchical tree for feature clustering 
# (row clustering shown on the left side of the heatmap) is flipped to improve the asthetics
# and readability, without changing the results or the samples clusters.

# inspect the sample hierarchical clustering result
TIL_moran_sample_dist <- dist(D2_Visium_TIL_moran)
TIL_moran_sample_hc <- hclust(TIL_moran_sample_dist, method="ward.D2") 
TIL_moran_sample_hc <- as.phylo(TIL_moran_sample_hc)

png("TIL_moran_sample_phylo.png", width = 2000, height = 2000, res = 300)
    par(mar = c(0.5,0.5,0.5,0.5))
    plot(TIL_moran_sample_hc, show.tip.label = TRUE, cex=0.4)
    nodelabels(text = 1:TIL_moran_sample_hc$Nnode + + Ntip(TIL_moran_sample_hc), 
                node= 1:TIL_moran_sample_hc$Nnode + Ntip(TIL_moran_sample_hc), 
                cex=0.4)
dev.off()

# Based on the determined best cluster numbers (see STAR Methods in the manuscript), cut the hierarchical tree
# at the corresponding level, and obtian the corresponding internal node, here we give the results
TIL_internal_node <- c(60, 62, 66, 67, 70, 71, 64, 68, 69)
# record sample clusters
for(i in 1:length(TIL_internal_node)){
    desc <- Descendants(TIL_moran_sample_hc, TIL_internal_node[i], type = "tips")[[1]]
    sample_ID <- TIL_moran_sample_hc$tip.label[desc]
    sample_subcluter_df$TIL_node[which(sample_subcluter_df$sample_ID %in% sample_ID)] <- TIL_internal_node[i]
}

# Based on the biological meaning of each cluster, categorize them as four major classes
sample_subcluter_df$TIL_category <- NA
sample_subcluter_df$TIL_category[sample_subcluter_df$TIL_node %in% c(66,67,70)] <- "B_cell_dominant_infiltration"
sample_subcluter_df$TIL_category[sample_subcluter_df$TIL_node %in% c(71,64)] <- "T_cell_dominant_infiltration"
sample_subcluter_df$TIL_category[sample_subcluter_df$TIL_node %in% c(60)] <- "Limited_TIL_infiltration"
sample_subcluter_df$TIL_category[sample_subcluter_df$TIL_node %in% c(62,68,69)] <- "TIL_Exclusion"


## Macrophage spatial patterns ##
# Z score of bivariant Moran's I against five concentric tissue layers, capped at +/- 0.9
D2_Visium_Macrophage_moran <- readRDS("../data/D2_Visium_Macrophage_moran.rds")

# plot spatial pattern heatmap
pdf("Visium_Macrophage_moran_cluster.pdf",width=15,height=5)
pheatmap::pheatmap(t(D2_Visium_Macrophage_moran),
                   cluster_rows = T,
                   cluster_cols = T,
                   show_rownames=T,
                   show_colnames=T,
                   clustering_distance_rows = "euclidean", 
                   clustering_distance_cols = "euclidean", 
                   clustering_method = "ward.D",
                   color = c(colorRampPalette(c('#4575b4','white'))(100),
                               colorRampPalette(c('white','#d73027'))(100)),
                   annotation_col = D2_Visium_sample_annotation,
                   annotation_row = D2_Visium_feature_annotation,
                   annotation_colors = annotation_colors,
                   main = "Visium Macrophage Moran's I",
                   cellwidth = 12, cellheight = 12, border_color = "white",
                   fontsize = 10, fontsize_row = 8, fontsize_col = 8, angle_col = 45
                )
dev.off()
# note that in the final figures presented, the hierarchical tree for feature clustering 
# (row clustering shown on the left side of the heatmap) is flipped to improve the asthetics
# and readability, without changing the results or the samples clusters.

# inspect the sample hierarchical clustering result
Macrophage_moran_sample_dist <- dist(Macrophage_moran)
Macrophage_moran_sample_hc <- hclust(Macrophage_moran_sample_dist, method="ward.D")
Macrophage_moran_sample_hc <- as.phylo(Macrophage_moran_sample_hc)

png("Macrophage_moran_sample_phylo.png", width = 2000, height = 2000, res = 300)
    par(mar = c(0.5,0.5,0.5,0.5))
    plot(Macrophage_moran_sample_hc, show.tip.label = TRUE, cex=0.4)
    nodelabels(text = 1:Macrophage_moran_sample_hc$Nnode + + Ntip(Macrophage_moran_sample_hc), 
                node= 1:Macrophage_moran_sample_hc$Nnode + Ntip(Macrophage_moran_sample_hc), 
                cex=0.4)
dev.off()

# Based on the determined best cluster numbers (see STAR Methods in the manuscript), cut the hierarchical tree
# at the corresponding level, and obtian the corresponding internal node, here we give the results
Macrophage_internal_node <- c(56,58,60,61)
# record sample clusters
for(i in 1:length(Macrophage_internal_node)){
    desc <- Descendants(Macrophage_moran_sample_hc, Macrophage_internal_node[i], type = "tips")[[1]]
    sample_ID <- Macrophage_moran_sample_hc$tip.label[desc]
    sample_subcluter_df$Macrophage_node[which(sample_subcluter_df$sample_ID %in% sample_ID)] <- Macrophage_internal_node[i]
}

# Based on the biological meaning of each cluster, categorize them as four major classes
sample_subcluter_df$Mac_category <- NA
sample_subcluter_df$Mac_category[sample_subcluter_df$Macrophage_node %in% c(58,61)] <- "Mac_infiltration"
sample_subcluter_df$Mac_category[sample_subcluter_df$Macrophage_node %in% c(56,60)] <- "Mac_Exclusion"


## Fibroblast spatial patterns ##
# Z score of bivariant Moran's I against five concentric tissue layers, capped at +/- 0.9
D2_Visium_Fibroblast_moran <- readRDS("../data/D2_Visium_Fibroblast_moran.rds")

# plot spatial pattern heatmap
pdf("Visium_Fibroblast_moran_cluster.pdf",width=15,height=5)
pheatmap::pheatmap(t(D2_Visium_Fibroblast_moran),
                   cluster_rows = T,
                   cluster_cols = T,
                   show_rownames=T,
                   show_colnames=T,
                   clustering_distance_rows = "euclidean", 
                   clustering_distance_cols = "euclidean", 
                   clustering_method = "ward.D",
                   color = c(colorRampPalette(c('#4575b4','white'))(100),
                               colorRampPalette(c('white','#d73027'))(100)),
                   annotation_col = D2_Visium_sample_annotation,
                   annotation_row = D2_Visium_feature_annotation,
                   annotation_colors = annotation_colors,
                   main = "Visium Fibroblast Moran's I",
                   cellwidth = 12, cellheight = 12, border_color = "white",
                   fontsize = 10, fontsize_row = 8, fontsize_col = 8, angle_col = 45
                )
dev.off()
# note that in the final figures presented, the hierarchical tree for feature clustering 
# (row clustering shown on the left side of the heatmap) is flipped to improve the asthetics
# and readability, without changing the results or the samples clusters.

# inspect the sample hierarchical clustering result
Fib_moran_sample_dist <- dist(D2_Visium_Fibroblast_moran)
Fib_moran_sample_hc <- hclust(Fib_moran_sample_dist, method="ward.D")
Fib_moran_sample_hc <- as.phylo(Fib_moran_sample_hc)

png("Fib_moran_sample_phylo.png", width = 2000, height = 2000, res = 300)
    par(mar = c(0.5,0.5,0.5,0.5))
    plot(Fib_moran_sample_hc, show.tip.label = TRUE, cex=0.4)
    nodelabels(text = 1:Fib_moran_sample_hc$Nnode + + Ntip(Fib_moran_sample_hc), 
                node= 1:Fib_moran_sample_hc$Nnode + Ntip(Fib_moran_sample_hc), 
                cex=0.4)
dev.off()

# Based on the determined best cluster numbers (see STAR Methods in the manuscript), cut the hierarchical tree
# at the corresponding level, and obtian the corresponding internal node, here we give the results
Fib_internal_node <- c(59,58,60,61)
# record sample clusters
for(i in 1:length(Fib_internal_node)){
    desc <- Descendants(Fib_moran_sample_hc, Fib_internal_node[i], type = "tips")[[1]]
    sample_ID <- Fib_moran_sample_hc$tip.label[desc]
    sample_subcluter_df$Fib_node[which(sample_subcluter_df$sample_ID %in% sample_ID)] <- Fib_internal_node[i]
}

# Based on the biological meaning of each cluster, categorize them as four major classes
sample_subcluter_df$Fib_category <- NA
sample_subcluter_df$Fib_category[sample_subcluter_df$Fib_node %in% c(58,60)] <- "Fib_infiltration"
sample_subcluter_df$Fib_category[sample_subcluter_df$Fib_node %in% c(59,61)] <- "Fib_Exclusion"

# based on above analysis,update sample annotations accordingly
D2_Visium_sample_annotation$TIL_category <- sample_subcluter_df$TIL_category
D2_Visium_sample_annotation$Mac_category <- sample_subcluter_df$Mac_category
D2_Visium_sample_annotation$Fib_category <- sample_subcluter_df$Fib_category



### Defining spatial immune archetypes (IMPACT) ###
# construct co-association matrix, ie, binary distance matrix for each major TME cell type category
D_TIL <- outer(sample_subcluter_df$TIL_category, sample_subcluter_df$TIL_category, FUN = "!=") * 1 
D_Mac <- outer(sample_subcluter_df$Mac_category, sample_subcluter_df$Mac_category, FUN = "!=") * 1
D_Fib <- outer(sample_subcluter_df$Fib_category, sample_subcluter_df$Fib_category, FUN = "!=") * 1

# assign the weight for each major TME cell type based on biological importance
w1 = 5; w2 = 4; w3 = 3
# constructed weighted co-association matrix and perform sample clutering accordingly
coassoc <- (w1*(D_TIL == 0) + w2*(D_Mac == 0) + w3*(D_Fib == 0)) / (w1 + w2 + w3)
colnames(coassoc) <- sample_subcluter_df$sample_ID; rownames(coassoc) <- sample_subcluter_df$sample_ID
dist_mat <- 1 - coassoc
hc <- hclust(as.dist(dist_mat),method="ward.D")

# here we flip a few braches in the raw hierarchical tree to improve the asthetics, without changing the inherit hierarchy among samples or the clustering results
sample_coassoc_hier_clust <- ggtree(hc, layout="rectangular") +
                                geom_tree(linewidth = 0.5) +
                                geom_tiplab(angle = 90, hjust = 1.2, vjust = 0.5, size = 3) + 
                                geom_text(aes(label = node), hjust = 0.35, size = 2) +
                                coord_flip(clip = "off") +
                                scale_x_reverse() +
                                xlim_tree(3) +
                                theme_tree() +
                                theme(plot.margin = margin(t=5, r=5, b=5, l=5))
sample_coassoc_hier_clust <- sample_coassoc_hier_clust %>% flip(70,71) %>% flip(66,67)

pdf("sample_coassoc_hier_clust.pdf", height = 3.5, width = 15)
sample_coassoc_hier_clust
dev.off()
# this gives the final hierrarchical tree presented on the top part of Figure 2B

# extract the sample order, for the following heatmap to match the hierarchical tree
sample_order <- rev(get_taxa_name(sample_coassoc_hier_clust))
# select some features to present as heatmap, can be customized based on the features to display
spatial_moran_combine <- cbind(D2_Visium_TIL_moran, D2_Visium_Macrophage_moran, D2_Visium_Fibroblast_moran)[, c(1:3,5,13:16,20:21)]

pdf("combined_heatmap.pdf",width=15,height=5)
pheatmap::pheatmap(t(spatial_moran_combine)[,sample_order],
                   cluster_rows = T, cluster_cols = F,
                   show_rownames=T, show_colnames=T,
                   color = colorRampPalette(viridis(100, option="A"), bias = 0.8)(100),
                   annotation_col = D2_Visium_sample_annotation,
                   annotation_row = D2_Visium_feature_annotation,
                   annotation_colors = annotation_colors,
                   main = "Visium sample Moran's I",
                   cellwidth = 12, cellheight = 12, border_color = "white",
                   fontsize = 10, fontsize_row = 8, fontsize_col = 8, angle_col = 45)
dev.off()
# this gives the heatmap presented on the bottom part of Figure 2B