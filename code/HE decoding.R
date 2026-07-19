### H&E-based decoding of spatial tumor architectures ###
### load common required packages
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(viridis)
library(ComplexHeatmap)
library(circlize)

### In-house Lung WSI cohort ###
V3_Lung_WSI_patch_compos <- readRDS("../data/V3_Lung_WSI_patch_compos.rds")

# calculate abundance of TME components in the tumor region
patch_type <- unique(V3_Lung_WSI_patch_compos$patch_type)
lesion_type <- unique(V3_Lung_WSI_patch_compos$lesion_type)
patch_res <- matrix(NA, nrow=4, ncol=7)
rownames(patch_res) <- lesion_type
colnames(patch_res) <- patch_type

for(i in 1:length(patch_type)){
    patch = patch_type[i]
    for(j in 1:length(lesion_type)){
        lesion = lesion_type[j]
        df = filter(V3_Lung_WSI_patch_compos, patch_type == patch & lesion_type == lesion)
        patch_count = sum(df$patch_in_tumor)
        patch_res[j,i] <- patch_count
    }
}

# now we focus on key TME components
patch_res_TME <- patch_res[,c(2:4)]
res <- chisq.test(patch_res_TME)
ratio <- res$observed/res$expected

lesion_order <- c("AAH","AIS","MIA","ADC")
ratio <- ratio[lesion_order,]

col_min <- min(ratio)
col_max <- max(ratio)
col_fun = colorRamp2(c(col_min, 1, col_max), c("#2166ac", "white", "#b2182b"))

pdf("Lung_WSI_Roe.pdf",width=5,height=3)
Heatmap(t(ratio),
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        name = "Ro/e", 
        row_names_side = "right",
        row_names_gp = gpar(fontsize = 14),
        column_names_side = "bottom",
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 14),
        column_title = "Enrichment in the tumor region",
        column_title_gp = gpar(fontsize = 18, fontface = "bold"),
        rect_gp = gpar(col = "white", lwd = 2)
       )
dev.off()

# Compare each TME component individually
V3_Lung_WSI_patch_compos_wide <- V3_Lung_WSI_patch_compos %>% spread(patch_type, patch_in_tumor)
V3_Lung_WSI_patch_compos_wide$tumor_region_TME_total <- rowSums(V3_Lung_WSI_patch_compos_wide[,c(4:6,8)])
V3_Lung_WSI_patch_compos_wide$Fib_TME_ratio <- V3_Lung_WSI_patch_compos_wide$`Fibrous tissue` / V3_Lung_WSI_patch_compos_wide$tumor_region_TME_total
V3_Lung_WSI_patch_compos_wide$Mac_TME_ratio <- V3_Lung_WSI_patch_compos_wide$Macrophage / V3_Lung_WSI_patch_compos_wide$tumor_region_TME_total
V3_Lung_WSI_patch_compos_wide$Lymph_TME_ratio <- V3_Lung_WSI_patch_compos_wide$Lymphocyte / V3_Lung_WSI_patch_compos_wide$tumor_region_TME_total
V3_Lung_WSI_patch_compos_wide$lesion_type <- factor(V3_Lung_WSI_patch_compos_wide$lesion_type, levels=c("AAH","AIS","MIA","ADC"))

# here we plot lymphoctye-enriched patches as an example
pdf("Lymphocyte_TME_prop.pdf", width = 4, height = 4)
ggplot(V3_Lung_WSI_patch_compos_wide, aes(x = lesion_type, y = Lymph_TME_ratio, color = lesion_type)) +
    geom_boxplot(width=0.5, outlier.shape = NA) +
    geom_jitter(width=0.1, alpha = 0.8, size=1) +
    theme_classic() +
    xlab("Lesion Type") +
    ylab("Prop in TME") +
    ggtitle("Lymphocyte-enriched patch") +
    theme(legend.position = "bottom")+
    scale_color_manual(values = c("AAH"="#e6194b","AIS"="#3cb44b","MIA"="#ffe119","ADC"="#4363d8"))
dev.off()



### H&E images in the TCGA cohort ###
# We applied OpenTME analysis results to understand the TME compositions in each signature-baed IMPACT archetype.
# OpenTME is a foundation model based framework that generate single-cell-level spatial maps on H&E images
# Load data
V4_TCGA_HE_cell_num <- readRDS("../data/V4_TCGA_HE_cell_num.rds")

# calculate cell enrichment in the tumor region, we we use the carcinoma region and stroma region segmented by OpenTME together as the tumor region
cell_type <- unique(V4_TCGA_HE_cell_num$cell_type)
IMPACT_subtype <- unique(V4_TCGA_HE_cell_num$IMPACT_assignment)
cell_res <- matrix(NA, nrow=4, ncol=8)
rownames(cell_res) <- IMPACT_subtype
colnames(cell_res) <- cell_type

for(i in 1:length(cell_type)){
    cell_id = cell_type[i]
    for(j in 1:length(IMPACT_subtype)){
        IMPACT_id = IMPACT_subtype[j]
        df = filter(V4_TCGA_HE_cell_num, IMPACT_assignment == IMPACT_id & cell_type == cell_id)
        cell_count = sum(df$Tumor_region)
        cell_res[j,i] <- cell_count
    }
}

res <- chisq.test(cell_res)
ratio <- res$observed/res$expected
col_min <- min(ratio)
col_max <- max(ratio)

col_fun = colorRamp2(c(col_min, 1, col_max), c("#2166ac", "white", "#b2182b"))
  
pdf("TCGA_HE_Roe.pdf",width=10,height=5)
Heatmap(ratio,
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        name = "Ro/e",
        row_names_side = "right",
        row_names_gp = gpar(fontsize = 14),
        column_names_side = "bottom",
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 14),
        column_title = "Enrichment of cell types in tumor & stroma",
        column_title_gp = gpar(fontsize = 18, fontface = "bold"),
        rect_gp = gpar(col = "white", lwd = 2)
       )
dev.off()

# Next, we construct the single-cell spatial map as pseudo-bulk and calculate spatial enrichment of different cell types (bivariant Moran's I)
V4_TCGA_HE_moran <- readRDS("../data/V4_TCGA_HE_moran.rds")

# calculate group median
V4_TCGA_HE_moran_median <- V4_TCGA_HE_moran %>% 
                            group_by(IMPACT_assignment) %>% 
                            summarise(across(where(is.numeric), 
                                      median, na.rm = TRUE)) %>% 
                            as.data.frame() %>% 
                            column_to_rownames(var="IMPACT_assignment")

# scale data for visualization
V4_TCGA_HE_moran_median_scale <- scale(V4_TCGA_HE_moran_median)
V4_TCGA_HE_moran_median_scale[V4_TCGA_HE_moran_median_scale < -0.9] <- -0.9
col_min <- abs(min(V4_TCGA_HE_moran_median_scale))*100
col_max <- abs(max(V4_TCGA_HE_moran_median_scale))*100

pdf("TCGA_moran_heatmap.pdf",width=10,height=5)
pheatmap::pheatmap(t(V4_TCGA_HE_moran_median_scale),
                   cluster_rows = F, cluster_cols = F,
                   show_rownames=T, show_colnames=T,
                   clustering_distance_rows = "euclidean",
                   clustering_method = "ward.D",
                   color = c(colorRampPalette(c('#4575b4','white'))(col_min),
                             colorRampPalette(c('white','#d73027'))(col_max)),
                   main = "TCGA HE Moran's I",
                   cellwidth = 30, cellheight = 30, border_color = "white",
                   fontsize = 10, fontsize_row = 12, fontsize_col = 12, angle_col = 45
                  )
dev.off()