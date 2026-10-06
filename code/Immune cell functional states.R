### Immune cell functional states of different IMPACT archetypes ###

# load common packages #
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(viridis)
library(ComplexHeatmap)
library(circlize) 
library(ggbeeswarm)


### Analyzing cell states across different IMPACT archetypes using snRNA-seq data ###
# load data
D1_snRNAseq_cell_metadata <- readRDS("../data/D1_snRNAseq_cell_metadata.rds")


# distribution of major cell types across different lesions
cell_res_major <- table(D1_snRNAseq_cell_metadata$IMPACT, D1_snRNAseq_cell_metadata$major_cell_type)
chisq_res_major <- chisq.test(cell_res_major)
Roe_major <-chisq_res_major$observed / chisq_res_major$expected

col_min <- min(Roe_major)
col_max <- max(Roe_major)
col_fun_major = colorRamp2(c(col_min, 1, col_max), c("#2166ac", "white", "#b2182b"))

pdf("snRNA_seq_Roe_major_IMPACT.pdf",width=6,height=5)
Heatmap(Roe_major,column_title = "snRNAseq major cell type by IMPACT Ro/e", 
        name = "Ro/e", col = col_fun_major,
        cluster_rows = TRUE,  cluster_columns = TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
            val = Roe_major[i, j]
            label = "-"
            if (val >= 2) {
            label = "+++"
            } else if (val >= 1.5) {
            label = "++"
            } else if (val >= 1) {
            label = "+"
            }
            grid.text(label, x, y, gp = gpar(fontsize = 10, col = "black"))
        })
dev.off()

# distribution of immune cell states across different lesions
D1_snRNAseq_cell_metadata_immune <- filter(D1_snRNAseq_cell_metadata, D1_snRNAseq_cell_metadata$major_cell_type %in% c("TNK","BPC","Myeloid"))
cell_res_immune <- table(D1_snRNAseq_cell_metadata_immune$IMPACT, D1_snRNAseq_cell_metadata_immune$cell_state)
chisq_res_immune <- chisq.test(cell_res_immune)
Roe_immune <-chisq_res_immune$observed / chisq_res_immune$expected

col_min <- min(Roe_immune)
col_max <- max(Roe_immune)
col_fun_immune = colorRamp2(c(col_min, 1, col_max), c("#276419", "white", "#8e0152"))

pdf("snRNA_seq_Roe_immune_IMPACT.pdf",width=15,height=5)
Heatmap(Roe_immune,column_title = "snRNAseq immune cell state by IMPACT Ro/e", 
        name = "Ro/e", col = col_fun_immune,
        cluster_rows = FALSE,  cluster_columns = TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
            val = Roe_immune[i, j]
            label = "-"
            if (val >= 2) {
            label = "+++"
            } else if (val >= 1.5) {
            label = "++"
            } else if (val >= 1) {
            label = "+"
            }
            grid.text(label, x, y, gp = gpar(fontsize = 10, col = "black"))
        })
dev.off()

# distribution of stromal cell states across different lesions
D1_snRNAseq_cell_metadata_stromal <- filter(D1_snRNAseq_cell_metadata, D1_snRNAseq_cell_metadata$major_cell_type %in% c("Stromal"))
cell_res_stromal <- table(D1_snRNAseq_cell_metadata_stromal$IMPACT, D1_snRNAseq_cell_metadata_stromal$cell_state)
chisq_res_stromal <- chisq.test(cell_res_stromal)
Roe_stromal <-chisq_res_stromal$observed / chisq_res_stromal$expected

col_min <- min(Roe_stromal)
col_max <- max(Roe_stromal)
col_fun_stromal = colorRamp2(c(col_min, 1, col_max), c("#276419", "white", "#8e0152"))

pdf("snRNA_seq_Roe_stromal_IMPACT.pdf",width=10,height=5)
Heatmap(Roe_stromal,column_title = "snRNAseq stromal cell state by IMPACT Ro/e", 
        name = "Ro/e", col = col_fun_stromal,
        cluster_rows = FALSE,  cluster_columns = TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
            val = Roe_stromal[i, j]
            label = "-"
            if (val >= 2) {
            label = "+++"
            } else if (val >= 1.5) {
            label = "++"
            } else if (val >= 1) {
            label = "+"
            }
            grid.text(label, x, y, gp = gpar(fontsize = 10, col = "black"))
        })
dev.off()



### Validate IMPACT signatures using pseudo-bulk profiles from snRNA-seq data ###
# load data
D1_pseudobulk_IMPACT_LUAD <- readRDS("../data/D1_pseudobulk_IMPACT_LUAD.rds")

# calculate average signature expression
LUAD_IMPACT_sig_z_score <- D1_pseudobulk_IMPACT_LUAD %>%  
                            group_by(Spatial_based_IMPACT) %>% 
                            summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>% 
                            as.data.frame() %>% 
                            column_to_rownames(var="Spatial_based_IMPACT") %>% 
                            scale()

pdf("IMPACT_sig_exp.pdf",width=5,height=4)
pheatmap::pheatmap(LUAD_IMPACT_sig_z_score,
         color = viridis::viridis(100, option="A"),
         cluster_rows = F,cluster_cols = F,
         main = "IMPACT sig exp",
         cellwidth = 25, cellheight = 25, border_color = "white", 
         fontsize = 10, fontsize_row = 8, fontsize_col = 8, angle_col = 45)
dev.off()



### Application of IMPACT signature in TCGA LUAD cohort ###
V4_TCGA_Sig_based_IMPACT <- readRDS("../data/V4_TCGA_Sig_based_IMPACT.rds")
  
# summarize sample compositions  
IMPACT_count <- as.data.frame(table(V4_TCGA_Sig_based_IMPACT$Sig_based_IMPACT))

pdf("TCGA_IMPACT_composition.pdf",width=3,height=3)
ggplot(IMPACT_count, aes("",Freq,fill=Var1))+
    geom_bar(stat = "identity",position = "fill")+
    coord_polar(theta = "y") + 
    theme_void() + 
    scale_fill_manual(values=c("#ff7f00","#b2df8a","#a6cee3","#6a3d9a"))
dev.off()

# TME-related signature score expression across different signature-based IMPACT archetypes
pdf("TLS_score_TCGA_IMPACT.pdf",width=3,height=3)
    ggplot(V4_TCGA_Sig_based_IMPACT, aes(x=Sig_based_IMPACT, y=TLS,color=Sig_based_IMPACT))+
    geom_boxplot(width=0.6,outlier.alpha=0)+
    geom_point(size=1, alpha=0.5,position="jitter")+
    scale_color_manual(values=c("#ff7f00","#b2df8a","#a6cee3","#6a3d9a"))+
    theme_classic()+
    theme(axis.text.x=element_text(angle = 45, hjust = 1))+
    theme(axis.ticks.length = unit(0.1,'cm'), legend.position = "none")
dev.off()

pdf("Macro_score_TCGA_IMPACT.pdf",width=3,height=3)
    ggplot(V4_TCGA_Sig_based_IMPACT, aes(x=Sig_based_IMPACT, y=Macro,color=Sig_based_IMPACT))+
    geom_boxplot(width=0.6,outlier.alpha=0)+
    geom_point(size=1, alpha=0.5,position="jitter")+
    scale_color_manual(values=c("#ff7f00","#b2df8a","#a6cee3","#6a3d9a"))+
    theme_classic()+
    theme(axis.text.x=element_text(angle = 45, hjust = 1))+
    theme(axis.ticks.length = unit(0.1,'cm'), legend.position = "none")
dev.off()

# correlation between IMPACT signatures and TME-related signatures
my_theme <- theme_bw()+
theme(aspect.ratio=1, legend.position="none",
    axis.text.x=element_text(angle = 0, hjust = 0.5),
    axis.ticks.length = unit(0.1,'cm'),
    panel.spacing = unit(0,"cm"))

blue_range <- function(n){
    colorRampPalette(c("#00204d","#00326f","#31446b","#4e576c","#666970","#958f78","#b0a473","#e7d159","#ffea46"))(n)
}

pdf("TCGA_IMPACT_TLS_cor.pdf",width=5,height=5)
ggplot(V4_TCGA_Sig_based_IMPACT, aes(x=IMPACT2, y=TLS))+
    coord_cartesian(expand=FALSE)+
    geom_density2d_filled(bins=10)+
    scale_fill_manual(values=blue_range(10))+
    geom_point(color="lightgrey", size=1, alpha=0.3)+
    geom_smooth(method = lm, se=FALSE, color="#ffffff", linetype="dashed")+
    my_theme
dev.off()

cor.test(V4_TCGA_Sig_based_IMPACT$IMPACT2, V4_TCGA_Sig_based_IMPACT$TLS, method = "spearman")

# similar analyses can be done for other IMPACT signatures and TME-related signatures



### TLS analysis ###
D2_Visium_TLS_summary <- readRDS("../data/D2_Visium_TLS_summary.rds")

# TLS area across different lesion types (TLS spot proportions among all spots)
D2_Visium_TLS_summary$lesion_type <- factor(D2_Visium_TLS_summary$lesion_type, levels = c("AAH", "AIS", "MIA", "LUAD"))
stats_df <- D2_Visium_TLS_summary %>% 
                group_by(lesion_type) %>% 
                summarise(median = median(TLS_spot_prop), 
                          q1 = quantile(TLS_spot_prop, 0.25),
                          q3 = quantile(TLS_spot_prop, 0.75))
            
pdf("TLS spot prop.pdf",width=3,height=3)
ggplot(D2_Visium_TLS_summary)+
    geom_quasirandom(width = 0.3, size = 1.8, alpha = 0.7, aes(x=lesion_type, y=TLS_spot_prop, color=lesion_type)) +
    geom_point(data = stats_df, aes(y = median, x = lesion_type), shape = 20, size = 3, color = "#969696") +
    geom_errorbar(data = stats_df, aes(ymin = q1, ymax = q3, x = lesion_type), width = 0.15, color = "#969696",linewidth = 0.8) +
    scale_color_manual(values=c("AAH"="#e31a1c","AIS"="#33a02c","MIA"="#FEE013","LUAD"="#1f78b4"))+
    ggtitle("TLS spot prop")+
    theme_classic()+
    theme(axis.text.x=element_text(angle = 45, hjust = 1))+
    theme(axis.ticks.length = unit(0.1,'cm'), legend.position = "none")
dev.off()

# We classified TLS regions as E-TLS and mTLS based on the presence of FDC signals (see STAR Methods), subsequently quantify TLS compositions #
D2_Visium_TLS_region_anno <- readRDS("../data/D2_Visium_TLS_region_anno.rds")

# TLS compositions in pre-invasive lesions
D2_Visium_TLS_region_anno_pre_invasive <- filter(D2_Visium_TLS_region_anno, D2_Visium_TLS_region_anno$lesion_class == "Pre_invasive")
Pre_invasive_TLS_composition <- as.data.frame(table(D2_Visium_TLS_region_anno_pre_invasive$Maturation_status))
    
pdf("Pre_invasive_TLS_composition.pdf",width=3,height=3)
ggplot(Pre_invasive_TLS_composition, aes("",Freq,fill=Var1))+
    geom_bar(stat = "identity",position = "fill")+
    coord_polar(theta = "y") + 
    theme_void() + 
    scale_fill_manual(values=c("E-TLS"="#80b1d3","mTLS"="#fdae61"))
dev.off()

# TLS compositions in invasive lesions
D2_Visium_TLS_region_anno_invasive <- filter(D2_Visium_TLS_region_anno, D2_Visium_TLS_region_anno$lesion_class == "Invasive")
Invasive_TLS_composition <- as.data.frame(table(D2_Visium_TLS_region_anno_invasive$Maturation_status))

pdf("Invasive_TLS_composition.pdf",width=3,height=3)
ggplot(Invasive_TLS_composition, aes("",Freq,fill=Var1))+
    geom_bar(stat = "identity",position = "fill")+
    coord_polar(theta = "y") + 
    theme_void() + 
    scale_fill_manual(values=c("E-TLS"="#80b1d3","mTLS"="#fdae61"))
dev.off()   

# TLS count and compositions in different IMPACT archetypes
pdf("IMPACT_TLS_count.pdf",width=4,height=3)
ggplot(filter(D2_Visium_TLS_summary, lesion_class == "Invasive"),
        aes(x=IMPACT, y=TLS_region_count, color=IMPACT))+
    geom_boxplot(width=0.5, outlier.alpha=0)+
    geom_jitter(width=0.1, size=1, alpha=0.7)+
    scale_color_manual(values=c("#ff7f00","#b2df8a","#a6cee3","#6a3d9a","#e31a1c","#253494"))+
    theme_classic()+
    theme(axis.ticks.length = unit(0.1,'cm'),legend.position = "right")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
# similar for plotting pre-invasive lesions or all lesions

TLS_compos_IMPACT <- as.data.frame(table(D2_Visium_TLS_region_anno$IMPACT, D2_Visium_TLS_region_anno$Maturation_status))
colnames(TLS_compos_IMPACT) <- c("IMPACT","Maturation_status","Num")

pdf("prop.pdf",width=7,height=3)
ggplot(TLS_compos_IMPACT, aes(Num, IMPACT, fill=Maturation_status))+
    geom_bar(stat = "identity",position = "fill")+
    scale_fill_manual(values=c("E-TLS"="#80b1d3","mTLS"="#fdae61"))+
    theme_classic()+
    theme(axis.ticks.length = unit(0.1,'cm'),legend.position = "right")+
    theme(axis.text.x=element_text(angle = 45, hjust = 1))+
    guides(fill=guide_legend(title = NULL))
dev.off()

# Spatial distribution of TLS spots 
D2_Visium_TLS_spot_distribution <- readRDS("../data/D2_Visium_TLS_spot_distribution.rds")

pdf("TLS_spot_distribution_Pre_invasive.pdf",width=3,height=3)
ggplot(filter(D2_Visium_TLS_spot_distribution, lesion_class=="Pre_invasive"), 
        aes(Maturation_satus, Proportion, fill=Tissue_layer))+
    geom_bar(stat = "identity",position = "stack")+
    scale_fill_manual(values=c("#d73027","#fee090","#abd9e9","#74add1","#4575b4"))+
    ggtitle("TLS spot distribution Pre-invasive")+
    labs(x="Maturation Status", y="Proportion among TLS spots")+
    theme_classic()+
    theme(axis.ticks.length = unit(0.1,'cm'),legend.position = "right",
            axis.text.x = element_text(angle = 45, hjust = 1))+
    guides(fill=guide_legend(title = NULL))
dev.off()

pdf("TLS_spot_distribution_Invasive.pdf",width=3,height=3)
ggplot(filter(D2_Visium_TLS_spot_distribution, lesion_class=="Invasive"), 
        aes(Maturation_satus, Proportion, fill=Tissue_layer))+
    geom_bar(stat = "identity",position = "stack")+
    scale_fill_manual(values=c("#d73027","#fee090","#abd9e9","#74add1","#4575b4"))+
    ggtitle("TLS spot distribution Invasive")+
    labs(x="Maturation Status", y="Proportion among TLS spots")+
    theme_classic()+
    theme(axis.ticks.length = unit(0.1,'cm'),legend.position = "right",
            axis.text.x = element_text(angle = 45, hjust = 1))+
    guides(fill=guide_legend(title = NULL))
dev.off()

# TLS spot distribution in different IMPACT archetypes
D2_Visium_TLS_spatial_moran <- readRDS("../data/D2_Visium_TLS_spatial_moran.rds")
D2_Visium_TLS_spatial_moran$IMPACT_tissue_layer <- paste0(D2_Visium_TLS_spatial_moran$IMPACT, "_", D2_Visium_TLS_spatial_moran$Tissue_layer)
IMPACT_tissue_layer_anno <- dplyr::select(D2_Visium_TLS_spatial_moran, c("IMPACT_tissue_layer", "IMPACT", "Tissue_layer")) %>% unique()
tissue_layer_order <- c("Normal_distal","Normal_peri","Interface_normal","Interface_tumor","Tumor_region")

Overall_TLS_median <- D2_Visium_TLS_spatial_moran %>% 
                        filter(TLS_states == "TLS" & lesion_class %in% c("Pre_invasive", "Invasive")) %>% 
                        # here "TLS" represents all TLS spots, E-TLS represents spots within E-TLS regions, mTLS represents spots within mTLS regions
                        group_by(IMPACT_tissue_layer) %>% 
                        summarise(across(where(is.numeric), median, na.rm = TRUE)) %>%
                        as.data.frame() %>%
                        left_join(IMPACT_tissue_layer_anno, by="IMPACT_tissue_layer") %>%
                        reshape2::dcast(formula = IMPACT ~ Tissue_layer, value.var="Moran") %>%
                        column_to_rownames(var="IMPACT")

Overall_TLS_median <- Overall_TLS_median[,tissue_layer_order]
col_low <- abs(min(Overall_TLS_median, na.rm = TRUE))*100
col_high <- abs(max(Overall_TLS_median, na.rm = TRUE))*100

pdf("Moran_all_TLS_all_samples.pdf",width=7,height=5)
pheatmap::pheatmap(Overall_TLS_median,
                    cluster_cols = F, cluster_rows = F, 
                    show_colnames = T, show_rownames = T,
                    color = c(colorRampPalette(c('#66BB6A','white'))(col_low),
                                colorRampPalette(c('white','#7E57C2'))(col_high)),
                    main = "Median Moran's I for all TLS spots across all samples",
                    cellwidth = 25, cellheight = 25, border_color = "white", 
                    fontsize = 10, fontsize_row = 8, fontsize_col = 8, angle_col = 45
                    )
dev.off()
# similar for plotting the distribution of E-TLS/mTLS spots within pre-invasive/invasive lesions