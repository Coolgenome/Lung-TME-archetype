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
# The H&E images in the TCGA cohort were analyzed by LungTME-Net, here we provide the analysis results
# Load data
V4_TCGA_HE_pred <- readRDS("../data/V4_TCGA_HE_pred.rds")

# plot the TR density of different patches
pdf("TCGA_IMPACT_patch_prop_lymph.pdf",width=5,height=5)
ggplot(filter(V4_TCGA_HE_pred, cell_type == "Lymphocytes"),
       aes(x = IMPACT_assignment, y = prop_tumor, color = IMPACT_assignment)) +
    geom_boxplot(width=0.5, outlier.shape = NA) +
    geom_jitter(width=0.1, alpha = 0.8, size=1) +
    theme_classic() +
    xlab("IMPACT") +
    ylab("Prop") +
    ggtitle("Lymph-enriched") +
    theme(legend.position = "bottom") +
    scale_color_manual(values = c("IMPACT1a"="#ff7f00","IMPACT1b"="#b2df8a",
                                  "IMPACT2"="#a6cee3","IMPACT3"="#6a3d9a"))
dev.off()

# Calculate the logFC of TME cells in the TR and plot
V4_TCGA_HE_pred$tumor_to_nontumor <- (V4_TCGA_HE_pred$prop_tumor + 1e-05) / (V4_TCGA_HE_pred$prop_nontumor + 1e-05)
V4_TCGA_HE_pred$log_ratio <- log2(V4_TCGA_HE_pred$tumor_to_nontumor)

V4_TCGA_HE_pred_wide <- V4_TCGA_HE_pred %>% 
                            select(slide_id, IMPACT_assignment, cell_type, log_ratio) %>% 
                            spread(cell_type, log_ratio)

V4_TCGA_HE_pred_wide_IMPACT_mean <- V4_TCGA_HE_pred_wide %>% 
                                        group_by(IMPACT_assignment) %>% 
                                        summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
                                        as.data.frame() %>%
                                        column_to_rownames(var="IMPACT_assignment")

V4_TCGA_HE_pred_wide_IMPACT_mean_scale <- scale(V4_TCGA_HE_pred_wide_IMPACT_mean[,1:3])

pdf("TCGA_IMPACT_TME_enrichment.pdf",width=6,height=5)
pheatmap::pheatmap(t(V4_TCGA_HE_pred_wide_IMPACT_mean_scale),
                    cluster_rows = F, cluster_cols = F,
                    show_rownames=T, show_colnames=T,
                    clustering_distance_rows = "euclidean",
                    clustering_method = "ward.D",
                    color =c(colorRampPalette(c('#1B1026','#4B145F'))(50),
                                colorRampPalette(c('#4B145F','#8C1D40'))(30),
                                colorRampPalette(c('#8C1D40','#D94F04'))(30),
                                colorRampPalette(c('#D94F04','#F9A602'))(50),
                                colorRampPalette(c('#F9A602','#FFF176'))(50),
                                colorRampPalette(c('#FFF176','#ffffcc'))(100)
                            ),
                    main = "TCGA TME enrichment in TR",
                    cellwidth = 30, cellheight = 30, border_color = "white",
                    fontsize = 12, fontsize_row = 12, fontsize_col = 12, angle_col = 45
                   )
dev.off()
