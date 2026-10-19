### Spatial features of different IMPACT archetypes ###

# load common packages #
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(viridis)

### Analysis of paired Xenium 5K data ###
# load data
D3_Xenium5K_pseudo_spot_moran <- readRDS("../data/D3_Xenium5K_pseudo_spot_moran.rds")
D3_Xenium5K_sample_anno <- readRDS("../data/D3_Xenium5K_sample_anno.rds")

## Moran's I heatmap for Xenium 5K data ##
pdf("Xenium5K_moran_heatmap.pdf",width=8,height=5)
pheatmap::pheatmap(D3_Xenium5K_pseudo_spot_moran,
                    cluster_cols=T, cluster_rows=F,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean",
                    clustering_method = "ward.D",
                    color = colorRampPalette(viridis(100, option="A"), bias = 1.2)(100),
                    annotation_row = D3_Xenium5K_sample_anno,
                    annotation_colors = list(
                        IMPACT = c("IMPACT1a"="#ff7f00","IMPACT1b"="#b2df8a",
                                   "IMPACT2"="#a6cee3","IMPACT3"="#6a3d9a")
                    ),
                    main = "Xenium5K Moran's I",
                    cellwidth = 16, cellheight = 16, border_color = "white",
                    fontsize = 12, fontsize_row = 10, fontsize_col = 10, angle_col = 45
                )
dev.off()

## calculate correlation with paried Visium ST data using key features ##
# 10 key features shown in Figure 2B are used here for correlation analysis
D2_Visium_key_feature_moran <- readRDS("../data/D2_Visium_key_feature_moran.rds") # scaled Moran's I

# scale features from Xenium data
D3_Xenium5K_pseudo_spot_moran_scale <- scale(D3_Xenium5K_pseudo_spot_moran)

# calculate pairwise correlations
cor_matrix <- cor(t(D3_Xenium5K_pseudo_spot_moran_scale), 
                  t(D2_Visium_key_feature_moran), 
                  method = "pearson", 
                  use = "complete.obs")

# plot correlation matrix
Xenium_sample_order <- c("Xenium_P11_AAH","Xenium_P17_AIS","Xenium_P11_ADC","Xenium_P17_ADC","Xenium_P24_ADC","Xenium_P4_ADC","Xenium_P21_ADC","Xenium_P21_AIS")
Visium_sample_order <- c("Visium_P11_AAH","Visium_P17_AIS","Visium_P11_ADC","Visium_P17_ADC","Visium_P24_ADC","Visium_P4_ADC","Visium_P21_ADC","Visium_P21_AIS1")
col_min <- abs(min(cor_matrix))*100
col_max <- abs(max(cor_matrix))*100

pdf("Feature_correlations.pdf",width=8,height=8)
pheatmap::pheatmap(cor_matrix[Xenium_sample_order,Visium_sample_order],
                   cluster_cols=F, cluster_rows=F,
                   color = c(colorRampPalette(c('#4575b4','white'))(col_min),
                             colorRampPalette(c('white','#d73027'))(col_max)),
                   main = "Feature correlations",
                   cellwidth = 16, cellheight = 16, border_color = "white", 
                   fontsize = 12, fontsize_row = 10, fontsize_col = 10, angle_col = 45
                  )
dev.off()



### Analysis of independent Xenium 5K data ###
# load data (z score of bivariant Moran's I, capped at +/- 0.9)
V1_WTS_Xenium5K_pesudo_spot_moran <- readRDS("../data/V1_WTS_Xenium5K_pesudo_spot_moran.rds")

V1_WTS_Xenium5K_pesudo_spot_moran_mat <- t(V1_WTS_Xenium5K_pesudo_spot_moran[,c(1:5)])
V1_WTS_Xenium5K_anno <- dplyr::select(V1_WTS_Xenium5K_pesudo_spot_moran, c("IMPACT_pred"))

# Moran's I heatmap
pdf("WTS_Xenium5K_moran_heatmap.pdf",width=8,height=5)
pheatmap::pheatmap(V1_WTS_Xenium5K_pesudo_spot_moran_mat,
                    cluster_cols=F, cluster_rows=F,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean",
                    clustering_method = "ward.D",
                    color = colorRampPalette(viridis(100, option="A"), bias = 1.2)(100),
                    annotation_col = V1_WTS_Xenium5K_anno,
                    annotation_colors = list(
                        IMPACT_pred = c("IMPACT1b-like"="#b2df8a",
                                        "IMPACT2-like"="#a6cee3",
                                        "IMPACT3-like"="#6a3d9a")
                    ),
                   main = "WTS_Xenium5K_moran_heatmap",
                   cellwidth = 16, cellheight = 16, border_color = "white", 
                   fontsize = 12, fontsize_row = 10, fontsize_col = 10, angle_col = 45
                  )
dev.off()

# correlation between validation (IMPACT-like) and discovery (IMPACT)
D2_Visium_IMPACT_features <- readRDS("../data/D2_Visium_IMPACT_features.rds")
D2_Visium_IMPACT_features <- D2_Visium_IMPACT_features[-1, ] %>% 
        remove_rownames() %>% 
        column_to_rownames(var = "IMPACT")

V1_WTS_Xenium5K_IMPACT_like_features <- V1_WTS_Xenium5K_pesudo_spot_moran %>%
        group_by(IMPACT_pred) %>%
        summarise(across(where(is.numeric), mean, na.rm=T)) %>%
        as.data.frame() %>%
        column_to_rownames(var = "IMPACT_pred")

cor_mat <- cor(t(V1_WTS_Xenium5K_IMPACT_like_features), t(D2_Visium_IMPACT_features))

# plot correlation matrix
pdf("Feature_correlation_V1_D2.pdf",width=4,height=3)
pheatmap::pheatmap(cor_mat,
                   cluster_cols=F, cluster_rows=F,
                   color = c(colorRampPalette(c('#08306b','#2171b5'))(500),
                             colorRampPalette(c('#2171b5','white'))(482),
                             colorRampPalette(c('white','#ffeda0'))(500),
                             colorRampPalette(c('#ffeda0','#ef3b2c'))(300),
                             colorRampPalette(c('#ef3b2c','#cb181d'))(132)),
                   main = "Feature_correlation_V1_D2",
                   cellwidth = 25, cellheight = 25, border_color = "white", 
                   fontsize = 10, fontsize_row = 8, fontsize_col = 8, angle_col = 45
                  )
dev.off()
