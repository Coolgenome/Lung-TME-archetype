### TME cellular landscape ###
### load common required packages
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)

### Analyzing TME cellular landscape of PML and LUAD from Xenium 5K data ###

## analyze TME cellular compositions
# load data
D3_Xenium5K_all_cell_merge_cluster <- readRDS("../data/D3_Xenium5K_all_cell_merge_cluster.rds")
# due to size limit, the D3_Xenium5K_all_cell_merge_cluster.rds file was splitted into two parts in the data folder.
D3_Xenium5K_cell_merge_cluster_TME <- filter(D3_Xenium5K_all_cell_merge_cluster, 
                                             Major_cell_type %in% c("T","NK","B","Plasma","Macro","Neutro","Mast","cDC","pDC","Langerhans","Myeloid_other",
                                                                    "Fibro","Chondrocyte","Vas_Endo","Lym_Endo","Pericyte","SMC"))

# calculate cell proportion
cell.prop <- as.data.frame(prop.table(table(D3_Xenium5K_cell_merge_cluster_TME$sample_ID, D3_Xenium5K_cell_merge_cluster_TME$cell_state),margin=1))
colnames(cell.prop)<-c("sample_ID","cell_state","Proportion")
sample_anno <- dplyr::select(D3_Xenium5K_cell_merge_cluster_TME, c("sample_ID","lesion_type")) %>% unique()
cell.prop <- left_join(cell.prop, sample_anno, by="sample_ID")

# plot
pdf("Xenium5K_CD4_Treg_FOXP3_cell_prop.pdf",width=3,height=2)
  ggplot(filter(cell.prop, cell.prop$cell_state == "CD4_Treg_FOXP3"),
         aes(lesion_type,Proportion,color=lesion_type))+
    geom_boxplot(outlier.alpha=0,aes(color=lesion_type))+
    geom_point(size=2)+
    ggtitle("CD4_Treg_FOXP3")+
    labs(x = "Lesion type", y = "Prop. in TME")+
    scale_color_manual(values=c("#e51a4c","#faaf40","#4f65af"))+
    theme_classic()+
    theme(axis.ticks.length = unit(0.1,'cm'),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 10))+
    guides(fill=guide_legend(title = NULL))
dev.off()

# similar plots can be done for other cell states

## analyze TME cellular distributions across five concentric tissue layers
# load cell loc and tissue layer data
D3_Xenium5K_cell_coord_tissue_layer <- readRDS("../data/D3_Xenium5K_cell_coord_tissue_layer.rds")
# due to size limit, the D3_Xenium5K_cell_coord_tissue_layer.rds file was splitted into four parts in the data folder.

D3_Xenium5K_cell_merge_cluster_TME <- D3_Xenium5K_cell_merge_cluster_TME %>% 
                                        rownames_to_column(var="cell_ID") %>% 
                                        left_join(D3_Xenium5K_cell_coord_tissue_layer, by="cell_ID")

# note: for P19, 2 tissue sections for each lesion were profiled, and only 1 for each had H&E and thus pathology-based tissue segmentations accordingly

D3_Xenium5K_cell_merge_cluster_TME$sample_layer <- paste0(D3_Xenium5K_cell_merge_cluster_TME$sample_ID, ":", D3_Xenium5K_cell_merge_cluster_TME$tissue_layer)
cell_porportion_in_layer_TME <- as.data.frame(prop.table(table(D3_Xenium5K_cell_merge_cluster_TME$sample_layer, D3_Xenium5K_cell_merge_cluster_TME$cell_state), margin=1))
colnames(cell_porportion_in_layer_TME) <- c("sample_layer","cell_state","Proportion")
cell_porportion_in_layer_TME <- separate(cell_porportion_in_layer_TME, col="sample_layer", into=c("sample_ID","tissue_layer"), sep=":")
cell_porportion_in_layer_TME$tissue_layer <- factor(cell_porportion_in_layer_TME$tissue_layer, levels=c("Normal_distal","Normal_peri","Interface_normal","Interface_tumor","Tumor_region"))
cell_porportion_in_layer_TME <- na.omit(cell_porportion_in_layer_TME)

# plot PML and LUAD separately
library(splines)
tissue_layer_order <- c("Normal_distal", "Normal_peri", "Interface_normal", "Interface_tumor", "Tumor_region")

PML <- c("S0040805_AAH","S0041306_AAH","S0042000_AAH","S0041486_AIS","S0041490_AIS","S0041515_AIS")
LUAD <- c("S0040805_ADC","S0041306_ADC","S0042000_ADC","S0041486_ADC","S0041490_ADC","S0041515_ADC")

sample_color <- c("S0040805_AAH"="#fb8072", "S0040805_ADC"="#e31a1c", "S0041306_AAH"="#a6cee3", "S0041306_ADC"="#0570b0", "S0041486_AIS"="#b2df8a", "S0041486_ADC"="#33a02c", 
                "S0041490_AIS"="#fdbf6f", "S0041490_ADC"="#ff7f00", "S0041515_AIS"="#cab2d6", "S0041515_ADC"="#54278f", "S0042000_AAH"="#d9d9d9", "S0042000_ADC"="#969696")

PML_df <- filter(cell_porportion_in_layer_TME, sample_ID %in% PML & cell_state == "CD4_Treg_FOXP3")
PML_df$tissue_layer_num = as.numeric(PML_df$tissue_layer)

pdf("Xenium5K_CD4_Treg_FOXP3_PML_trend.pdf",width=7,height=5)
ggplot(PML_df, 
       aes(x = tissue_layer_num, y = Proportion, group= sample_ID)) +
       geom_point(size=1, aes(color = sample_ID)) +
       stat_smooth(method = "lm", formula = y ~ ns(x, 3), se = FALSE, aes(color = sample_ID))+
       scale_x_continuous(breaks = 1:5, labels = tissue_layer_order) +
       ggtitle("CD4_Treg_FOXP3 in PML")+
       labs(x = "Tissue Layer", y = "Prop. in TME") +
       scale_color_manual(values=sample_color)+
       theme_classic()+
       theme(axis.ticks.length = unit(0.1,'cm'),
             axis.text.x = element_text(angle = 45, hjust = 1),
             plot.title = element_text(size = 10))
dev.off()


LUAD_df <- filter(cell_porportion_in_layer_TME, sample_ID %in% LUAD & cell_state == "CD4_Treg_FOXP3")
LUAD_df$tissue_layer_num = as.numeric(LUAD_df$tissue_layer)

pdf("Xenium5K_CD4_Treg_FOXP3_LUAD_trend.pdf",width=7,height=5)
ggplot(LUAD_df, 
       aes(x = tissue_layer_num, y = Proportion, group= sample_ID)) +
       geom_point(size=1, aes(color = sample_ID)) +
       stat_smooth(method = "lm", formula = y ~ ns(x, 3), se = FALSE, aes(color = sample_ID))+
       scale_x_continuous(breaks = 1:5, labels = tissue_layer_order) +
       ggtitle("CD4_Treg_FOXP3 in LUAD")+
       labs(x = "Tissue Layer", y = "Prop. in TME") +
       scale_color_manual(values=sample_color)+
       theme_classic()+
       theme(axis.ticks.length = unit(0.1,'cm'),
             axis.text.x = element_text(angle = 45, hjust = 1),
             plot.title = element_text(size = 10))
dev.off()

# similar plots can be done for other cell states

# For rate of change of cell states across different distances relative to the tumor boundary, fit spline per sample and get derivatives as an estimate
x_grid <- seq(1, 5, length.out = 100)
PML_deriv_df <- PML_df %>%
    group_by(sample_ID) %>%
    nest() %>%
    mutate(
        model = map(data, ~ lm(Proportion ~ ns(tissue_layer_num, df = 3), data = .x)),
        pred = map(model, ~ predict(.x, newdata = data.frame(tissue_layer_num = x_grid))),
        derivative = map(pred, ~ c(NA, diff(.x) / diff(x_grid)))  # numerical first derivative
    ) %>%
    select(sample_ID, derivative) %>%
    unnest(derivative) %>%
    mutate(tissue_layer_num = rep(x_grid, times = n()/length(x_grid)))

pdf("CD4_Treg_FOXP3 rate of change in PML.pdf",width=7,height=5)
ggplot(PML_deriv_df, aes(x = tissue_layer_num, y = derivative, color = sample_ID)) +
        geom_line(linewidth=1.5) +
        labs(title = "CD4_Treg_FOXP3 rate of change in PML", y = "d(Proportion)/dx", x = "Tissue Layer") +
        scale_x_continuous(breaks = 1:5, labels = tissue_layer_order) +
        scale_color_manual(values=sample_color)+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


LUAD_deriv_df <- LUAD_df %>%
    group_by(sample_ID) %>%
    nest() %>%
    mutate(
    model = map(data, ~ lm(Proportion ~ ns(tissue_layer_num, df = 3), data = .x)),
    pred = map(model, ~ predict(.x, newdata = data.frame(tissue_layer_num = x_grid))),
    derivative = map(pred, ~ c(NA, diff(.x) / diff(x_grid)))  # numerical first derivative
    ) %>%
    select(sample_ID, derivative) %>%
    unnest(derivative) %>%
    mutate(tissue_layer_num = rep(x_grid, times = n()/length(x_grid)))

pdf("CD4_Treg_FOXP3 rate of change in LUAD.pdf",width=7,height=5)
ggplot(LUAD_deriv_df, aes(x = tissue_layer_num, y = derivative, color = sample_ID)) +
        geom_line(linewidth=1.5) +
        labs(title = "CD4_Treg_FOXP3 rate of change in LUAD", y = "d(Proportion)/dx", x = "Tissue Layer") +
        scale_x_continuous(breaks = 1:5, labels = tissue_layer_order) +
        scale_color_manual(values=sample_color)+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# similar plots can be done for other cell states



### Analyzing TME cellular landscape of PML and LUAD from Xenium TMA data ###
## TME cellular compositions of different lesion types
# load data
V2_PCA_Xenium_TME_metadata <- readRDS("../data/V2_PCA_Xenium_TME_metadata.rds")

# calculate enrichment
cell_state_table <- table(V2_PCA_Xenium_TME_metadata$Pathology,V2_PCA_Xenium_TME_metadata$cell_state)
chisq_res <- chisq.test(cell_state_table)
Roe_ratio <- chisq_res$observed/chisq_res$expected
lesion_type_order <- c("No_lesion","AAH","AIS","MIA","LUAD_lepidic","LUAD_acinar","LUAD_papillary","LUAD_solid")
Roe_ratio <- Roe_ratio[lesion_type_order,]

# plot
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, 1, 2.5), c("#2166ac", "white", "#b2182b"))

pdf("Xenium_TMA_TME_Roe.pdf",width=15,height=5)
Heatmap(Roe_ratio,  column_title = "Xenium TMA TME cell state Ro/e", name = "Ro/e",
        cluster_columns = TRUE, cluster_rows = FALSE, col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
            val = Roe_ratio[i, j]
            label = "-"
            if (val >= 5) {
            label = "+++"
            } else if (val >= 3) {
            label = "++"
            } else if (val >= 1) {
            label = "+"
            }
            grid.text(label, x, y, gp = gpar(fontsize = 10, col = "black"))
        }
       )
dev.off()

## tumor cell neighborhood in different lesions
# load data
V2_tumor_epi_neighborhood_TME <- readRDS("../data/V2_tumor_epi_neighborhood_TME.rds")

Avg_epi_neighbor_count_TME <- V2_tumor_epi_neighborhood_TME %>% group_by(Pathology) %>% summarise(across(where(is.numeric), mean, na.rm = TRUE))
Avg_epi_neighbor_count_TME <- as.data.frame(Avg_epi_neighbor_count_TME) %>% column_to_rownames(var="Pathology")
Avg_epi_neighbor_count_TME_scale <- scale(Avg_epi_neighbor_count_TME)
#Avg_epi_neighbor_count_TME_scale_clean <- Avg_epi_neighbor_count_TME_scale %>% as.data.frame() %>% select(where(~ !all(is.na(.x))))
Avg_epi_neighbor_count_TME_scale[Avg_epi_neighbor_count_TME_scale > 1.5] <- 1.5
Avg_epi_neighbor_count_TME_scale[Avg_epi_neighbor_count_TME_scale < -1.5] <- -1.5
Avg_epi_neighbor_count_TME_scale <- Avg_epi_neighbor_count_TME_scale[lesion_type_order[-1],]

# plot
library(pheatmap)
pdf("Xenium TMA avg epi neighbor count TME scale.pdf",width=10,height=3.5)
pheatmap(Avg_epi_neighbor_count_TME_scale,
         cluster_cols = TRUE, cluster_rows = FALSE,
         color = c(colorRampPalette(c("#1b7837",'white'))(150), 
                   colorRampPalette(c('white',"#762a83"))(150)))
dev.off()



### Analyzing TME cellular landscape of PML and LUAD from snRNAseq data ###
# load data
D1_snRNAseq_cell_metadata <- readRDS("../data/D1_snRNAseq_cell_metadata.rds")

library(ComplexHeatmap)
library(circlize) 
lesion_order <- c("AAH","AIS","MIA","LUAD")

# distribution of major cell types across different lesions
cell_res_major <- table(D1_snRNAseq_cell_metadata$lesion_type, D1_snRNAseq_cell_metadata$major_cell_type)
chisq_res_major <- chisq.test(cell_res_major)
Roe_major <-chisq_res_major$observed / chisq_res_major$expected
Roe_major <- Roe_major[lesion_order,]

col_min <- min(Roe_major)
col_max <- max(Roe_major)
col_fun_major = colorRamp2(c(col_min, 1, col_max), c("#2166ac", "white", "#b2182b"))

pdf("snRNA_seq_Roe_major.pdf",width=6,height=5)
Heatmap(Roe_major,column_title = "snRNAseq major cell type Ro/e", 
        name = "Ro/e", col = col_fun_major,
        cluster_rows = FALSE,  cluster_columns = TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
            val = Roe_major[i, j]
            label = "-"
            if (val >= 3) {
            label = "+++"
            } else if (val >= 2) {
            label = "++"
            } else if (val >= 1) {
            label = "+"
            }
            grid.text(label, x, y, gp = gpar(fontsize = 10, col = "black"))
        })
dev.off()

# distribution of immune cell states across different lesions
D1_snRNAseq_cell_metadata_immune <- filter(D1_snRNAseq_cell_metadata, D1_snRNAseq_cell_metadata$major_cell_type %in% c("TNK","BPC","Myeloid"))
cell_res_immune <- table(D1_snRNAseq_cell_metadata_immune$lesion_type, D1_snRNAseq_cell_metadata_immune$cell_state)
chisq_res_immune <- chisq.test(cell_res_immune)
Roe_immune <-chisq_res_immune$observed / chisq_res_immune$expected
Roe_immune <- Roe_immune[lesion_order,]

col_min <- min(Roe_immune)
col_max <- max(Roe_immune)
col_fun_immune = colorRamp2(c(col_min, 1, col_max), c("#276419", "white", "#8e0152"))

pdf("snRNA_seq_Roe_immune.pdf",width=15,height=5)
Heatmap(Roe_immune,column_title = "snRNAseq immune cell state Ro/e", 
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
cell_res_stromal <- table(D1_snRNAseq_cell_metadata_stromal$lesion_type, D1_snRNAseq_cell_metadata_stromal$cell_state)
chisq_res_stromal <- chisq.test(cell_res_stromal)
Roe_stromal <-chisq_res_stromal$observed / chisq_res_stromal$expected
Roe_stromal <- Roe_stromal[lesion_order,]

col_min <- min(Roe_stromal)
col_max <- max(Roe_stromal)
col_fun_stromal = colorRamp2(c(col_min, 1, col_max), c("#276419", "white", "#8e0152"))

pdf("snRNA_seq_Roe_stromal.pdf",width=10,height=5)
Heatmap(Roe_stromal,column_title = "snRNAseq stromal cell state Ro/e", 
        name = "Ro/e", col = col_fun_stromal,
        cluster_rows = FALSE,  cluster_columns = TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
            val = Roe_stromal[i, j]
            label = "-"
            if (val >= 3) {
            label = "+++"
            } else if (val >= 2) {
            label = "++"
            } else if (val >= 1) {
            label = "+"
            }
            grid.text(label, x, y, gp = gpar(fontsize = 10, col = "black"))
        })
dev.off()