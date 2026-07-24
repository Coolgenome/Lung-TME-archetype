### Molecular prolfiles and immune-related cell communications across different lesion types and IMPACT archetypes ###
### load common required packages
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(viridis)

### NMF for identifying robust metaprograms in the Visium ST discovery cohort ###
#load data
D2_Visium_NMF_program_jaccard <- readRDS("../data/D2_Visium_NMF_program_jaccard.rds")
D2_Visium_NMF_program_anno <- readRDS("../data/D2_Visium_NMF_program_anno.rds")

custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

anno_colors = list("Patient_ID"=c("P1"="#FF0000","P2"="#00FF00","P3"="#0000FF","P4"="#FFFF00",
                                  "P5"="#00FFFF","P6"="#FF00FF","P7"="#FFA500","P8"="#800080",
                                  "P9"="#FFC0CB","P10"="#A52A2A","P11"="#32CD32","P12"="#808000",
                                  "P13"="#800000","P14"="#000080","P15"="#008080","P16"="#40E0D0",
                                  "P17"="#EE82EE","P18"="#FFD700","P19"="#FA8072","P20"="#FF7F50",
                                  "P21"="#4B0082","P22"="#7FFF00","P23"="#DC143C","P24"="#F0E68C",
                                  "P25"="#E6E6FA"),
                   "Lesion_type"=c("AAH"="#e6194b","AIS"="#3cb44b","MIA"="#ffe119","LUAD"="#4363d8"))

pdf("D2_Visium_NMF_program_jaccard.pdf",width=6,height=5)
pheatmap::pheatmap(D2_Visium_NMF_program_jaccard, 
                   cluster_rows = F, cluster_cols = F,
                   show_rownames = F, show_colnames = F,
                   color=custom_magma,
                   annotation_col = D2_Visium_NMF_program_anno,
                   annotation_colors = anno_colors
                  )
dev.off()



### Compositions of MPs within the tumor region ###
# Each spot was assigned to the MP with the highest score, and the composition of MPs within the tumor region was calculated for each sample.
# load data
D2_Visium_MP_assignment <- readRDS("../data/D2_Visium_MP_assignment.rds")
D2_Visium_sample_IMPACT <- readRDS("../data/D2_Visium_sample_IMPACT.rds") 

# calculate the proportion of MPs within the tumor region for each sample
D2_Visium_MP_assignment_TR <- filter(D2_Visium_MP_assignment, Tissue_layer=="Tumor_region")
MP_prop_TR <- as.data.frame(prop.table(table(D2_Visium_MP_assignment_TR$sample_ID, D2_Visium_MP_assignment_TR$MP_assignment), margin=1))
colnames(MP_prop_TR) <- c("sample_ID","MP","Prop")
MP_prop_TR <- left_join(MP_prop_TR, D2_Visium_sample_IMPACT, by="sample_ID")
MP_prop_TR_wide <- reshape2::dcast(MP_prop_TR, formula = sample_ID+IMPACT ~ MP, value.var="Prop")

# plot the proportion of MPs across different IMPACT archetypes
pdf("MP4_prop_in_IMPACT.pdf",width=3.5,height=3)
ggplot(MP_prop_TR_wide, aes(IMPACT,MP4,color=IMPACT))+
    geom_boxplot(width=0.6, outlier.alpha=0)+
    geom_point(size=2,position="jitter")+
    theme_classic()+
    scale_color_manual(values=c("IMPACT1a"="#f47e20","IMPACT1b"="#b3d789","IMPACT2"="#a5cde1","IMPACT3"="#6b3e97"))+
    theme(axis.ticks.length = unit(0.1,'cm'),legend.position = "right")+
    theme(axis.text.x = element_text(angle=45,hjust=1))+
    guides(fill=guide_legend(title = NULL))
dev.off()



### Spatial expression gradient in myeloid-enriched spots (MP5 and MP7) across different lesion types and IMPACT archetypes ###
# We will examine the spatial expression gradient of myeloid-related genes across different lesion types and IMPACT archetypes.
# We specifically look at MP5 and MP7 spots, which have a higher abundance of myeloid populations.
# load data
D2_Visium_myeloid_gene_IMPACT <- readRDS("../data/D2_Visium_myeloid_gene_IMPACT.rds")
D2_Visium_myeloid_gene_lesion <- readRDS("../data/D2_Visium_myeloid_gene_lesion.rds")

# To compare the expression of CCR2 and CX3CR1 across different IMPACT archetypes, we will plot the mean signature score of these 
# two genes across different tissue layers for each IMPACT archetype. Here as an example we plot IMPACT1a, similar for others.
pdf("CCR2_CX3CR1_exp_IMPACT3.pdf",width=7,height=5)
ggplot(filter(D2_Visium_myeloid_gene_IMPACT, IMPACT == "IMPACT3" & feature %in% c("CCR2","CX3CR1")),
            aes(x = Tissue_layer, y = mean, group = feature, color = feature, fill = feature)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = mean-sem, ymax = mean+sem), alpha = 0.2, color = NA) +
    labs(x = "Layer", y = "Signature score (mean ± sem)",
         title = "IMPACT3",
         color = "feature", fill = "feature") +
    scale_color_manual(values=c("#a6cee3","#ff7f00"))+
    scale_fill_manual(values=c("#a6cee3","#ff7f00"))+    
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 12))
dev.off()

# To compare the expression of myeloid-related genes across different IMPACT archetypes, we will plot the mean signature score of these
# genes across different tissue layers. Here as an example we plot IL1B, similar for other genes.
pdf("IL1B_exp_IMPACT.pdf",width=7,height=5)
ggplot(filter(D2_Visium_myeloid_gene_IMPACT, feature == "IL1B"),
            aes(x = Tissue_layer, y = mean, group = IMPACT, color = IMPACT, fill = IMPACT)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = mean-sem, ymax = mean+sem), alpha = 0.2, color = NA) +
    labs(x = "Layer", y = "Signature score (mean ± sem)",
            title = "IL1B",
            color = "IMPACT", fill = "IMPACT") +
    scale_color_manual(values=c("IMPACT1a"="#ff7f00","IMPACT1b"="#b2df8a","IMPACT2"="#a6cee3","IMPACT3"="#6a3d9a"))+
    scale_fill_manual(values=c("IMPACT1a"="#ff7f00","IMPACT1b"="#b2df8a","IMPACT2"="#a6cee3","IMPACT3"="#6a3d9a"))+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 12))
dev.off()

# Similarly, we can compare the expression of myeloid-related genes across different lesion types
pdf("IL1B_exp_lesion.pdf",width=7,height=5)
ggplot(filter(D2_Visium_myeloid_gene_lesion, feature == "IL1B"),
            aes(x = Tissue_layer, y = mean, group = lesion_type, color = lesion_type, fill = lesion_type)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = mean-sem, ymax = mean+sem), alpha = 0.2, color = NA) +
    labs(x = "Layer", y = "Signature score (mean ± sem)",
            title = "IL1B",
            color = "Lesion Type", fill = "Lesion Type") +
    scale_color_manual(values=c("AAH"="#e51a4c","AIS"="#3eb44b","MIA"="#fee013","ADC"="#4f65af"))+
    scale_fill_manual(values=c("AAH"="#e51a4c","AIS"="#3eb44b","MIA"="#fee013","ADC"="#4f65af"))+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 12))
dev.off()



### T cell functionality ###
# We analyze this from two aspects: 1) Ecosystem-derived signals that may influence T cell functionality; 2) T cell functional score expression.
# 1) Ecosystem-derived signals
# We will examine the expression of functional molecules that may influence T cell functionality, including chemotaxis, interferon response, checkpoint, 
# T cell functional modulation, antigen presentation machinery, and stress response.
Chemotaxis <- c("CCL2", "CCL3", "CCL4", "CCL5","CCL19","CCL21","CXCL9", "CXCL10", "CXCL11","CXCL12","CCR5","CCR7","CXCR3","CXCR4","CX3CR1")
Interferon_response <- c("IFIT2","IFIT3","IFI44L","OAS3","STAT1","IFIT1","MX1","IFI44","MX2","HERC5","LY6E","IFI6","IFI27","ISG15","SAMD9L","XAF1","EPSTI1","EIF2AK2","RSAD2")
Checkpoint <- c("PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","CD244", "BATF","LAYN","ENTPD1","TBX21","CD274","PDCD1LG2","LGALS9","CD47","CD24","SIRPA")
T_cell_funtional_modulation <- c("ICOS","TNFRSF9","TNFRSF4","TNFRSF18","IL2","IL2RA","LCK","ZAP70","LAT","ITK","TOX","TOX2","NR4A1","NR4A2","NR4A3","EOMES","CD69","CD80","CD86","CD40","CD40LG")
Antigen_presentation_machinery <- c("CD74","B2M","CIITA","NLRC5","RFX5","RFXANK","RFXAP","TAP1","TAP2","TAPBP","ERAP1","ERAP2","CALR","CANX","PDIA3","CTSB","CTSS","CTSL")
Stress_response <- c("HSPA1A","HSPA1B","HSP90AA1","HSP90B1","BAG3","FOS","JUN","JUNB","FOSB","HSPD1","HSPB1","HSPA6","HSPA8","HSPE1","HSPA5","DNAJA1","DNAJB4","DNAJA4","DNAJB1")

Funtional_molecules <- c(Chemotaxis, Interferon_response, Checkpoint, T_cell_funtional_modulation, Antigen_presentation_machinery, Stress_response)

# load data
# this is the mean expression of the above functional molecules across different tissue layers in each sample
D2_Visium_functional_molecules <- readRDS("../data/D2_Visium_functional_molecules.rds") 

# calculate the mean and sem of the above functional molecules across different tissue layers in each IMPACT archetype
D2_Visium_functional_molecules_group_mean <- D2_Visium_functional_molecules %>% 
                                                group_by(IMPACT, Tissue_layer) %>%
                                                dplyr::summarise(across(all_of(Funtional_molecules),
                                                                 list(mean = ~mean(.x, na.rm = TRUE), sem  = ~sd(.x, na.rm = TRUE) / sqrt(length(.x))), 
                                                                 .names = "{.col}_{.fn}"), .groups = "drop") %>%
                                                as.data.frame()

# extract the mean expression of the above functional molecules in the tumor region for each IMPACT archetype
D2_Visium_functional_molecules_group_mean_TR <- D2_Visium_functional_molecules_group_mean %>%
                                                    filter(Tissue_layer=="Tumor_region") %>%
                                                    select(c("IMPACT", paste0(Funtional_molecules, "_mean"))) %>%
                                                    column_to_rownames(var="IMPACT")
colnames(D2_Visium_functional_molecules_group_mean_TR) <- gsub("_mean","",colnames(D2_Visium_functional_molecules_group_mean_TR))

# scale the data for visualization, here we plot chemotaxis molecules as an example
D2_Visium_functional_molecules_group_mean_TR_scale <- scale(D2_Visium_functional_molecules_group_mean_TR)

pdf("Chemotaxis_TR_IMPACT.pdf",width=5,height=3)
pheatmap::pheatmap(D2_Visium_functional_molecules_group_mean_TR_scale[,Chemotaxis],
         cluster_rows = F, cluster_cols = T,
         show_rownames=T, show_colnames=T,
         clustering_distance_cols = "euclidean", 
         clustering_method = "average",
         color = c(colorRampPalette(c("#66BB6A", "white"))(200),
                   colorRampPalette(c('white',"#7E57C2"))(200)),
         main = "Chemotaxis molecules in Tumor region",
         cellwidth = 15, cellheight = 15, border_color = "white",
         fontsize = 10, fontsize_row = 8, fontsize_col = 8, angle_col = 45 
        )
dev.off()
# similar ways to plot for other functional molecules

# 2) T cell functional score expression
# T cell functional scores include the following: CD8_Tex, Tstr, CD8_Tisg, CD8_Teff, CD4_Treg, and CD4_Tfh.
T_sig <- c("CD8_Tex", "Tstr", "CD8_Tisg", "CD8_Teff", "CD4_Treg", "CD4_Tfh")

# load data
# We calculated T cell functional score expression across MP4 spots, and calculated the mean expression across different tissue layers in each sample.
D2_Visium_T_sig <- readRDS("../data/D2_Visium_T_sig.rds")

# calculate the mean and sem of T cell functional score expression across different tissue layers in each IMPACT archetype
D2_Visium_T_sig_group_mean <- D2_Visium_T_sig %>% 
                                group_by(IMPACT,Tissue_layer) %>%
                                dplyr::summarise(across(all_of(T_sig),
                                                 list(mean = ~mean(.x, na.rm = TRUE), sem  = ~sd(.x, na.rm = TRUE) / sqrt(length(.x))), 
                                                 .names = "{.col}_{.fn}"), .groups = "drop") %>%
                                as.data.frame()

# Extract the mean expression of T cell functional scores in MP4 spots within the tumor region for each IMPACT archetype
D2_Visium_T_sig_group_mean_TR <- D2_Visium_T_sig_group_mean %>%
                                    filter(Tissue_layer=="Tumor_region") %>%
                                    select(c("IMPACT", paste0(T_sig, "_mean"))) %>%
                                    column_to_rownames(var="IMPACT")
colnames(D2_Visium_T_sig_group_mean_TR) <- gsub("_mean","",colnames(D2_Visium_T_sig_group_mean_TR))

# Scale the data for visualization
D2_Visium_T_sig_group_mean_TR_scale <- scale(D2_Visium_T_sig_group_mean_TR)
col_min <- abs(as.integer(min(D2_Visium_T_sig_group_mean_TR_scale)*100))
col_max <- abs(as.integer(max(D2_Visium_T_sig_group_mean_TR_scale)*100))

pdf("T_sig_TR_IMPACT.pdf",width=6,height=3)
pheatmap::pheatmap(D2_Visium_T_sig_group_mean_TR_scale,
         cluster_rows = F, cluster_cols = T,
         show_rownames=T, show_colnames=T,
         clustering_distance_cols = "euclidean", 
         clustering_method = "average",
         color = c(colorRampPalette(c("#4575b4", "white"))(col_min),
                   colorRampPalette(c('white',"#d73027"))(col_max)),
         main = "T cell functional scores in MP4 spots of tumor region",
         cellwidth = 15, cellheight = 15, border_color = "white",
         fontsize = 10, fontsize_row = 8, fontsize_col = 8, angle_col = 45
        )
dev.off()



### Imune-related cell communications across different lesion types and IMPACT archetypes ###
## First we compare the relative contributions of immune-related signlaing pathways in TR across different IMPACT archetypes.
# IMPACT1a vs IMPACT1b
D2_Visium_IMPACT1a_vs_1b_information_flow <- readRDS("../data/D2_Visium_IMPACT1a_vs_1b_information_flow.rds")

pdf("IMPACT1a_vs_1b_signal_plot.pdf",width=5,height=4)
ggplot(D2_Visium_IMPACT1a_vs_1b_information_flow, 
       aes(x=pathway, y=Relative_information_flow, fill=group))+
    geom_bar(stat="identity", position="stack")+
    xlab("Immune-related signaling pathways in TR")+
    ylab("Relative information flow")+
    geom_hline(yintercept = 0.5, linetype="dashed", color = "red")+
    ggtitle("IMPACT1a vs IMPACT1b")+
    scale_fill_manual(values=c("IMPACT1a"="#ff7f00","IMPACT1b"="#b2df8a"))+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# IMPACT2 vs IMPACT3
D2_Visium_IMPACT2_vs_3_information_flow <- readRDS("../data/D2_Visium_IMPACT2_vs_3_information_flow.rds")

pdf("IMPACT2_vs_3_signal_plot.pdf",width=6.5,height=4)
ggplot(D2_Visium_IMPACT2_vs_3_information_flow, 
       aes(x=pathway, y=Relative_information_flow, fill=group))+
    geom_bar(stat="identity", position="stack")+
    xlab("Signaling pathways")+
    ylab("Contribution")+
    geom_hline(yintercept = 0.5, linetype="dashed", color = "red")+
    ggtitle("IMPACT2 vs IMPACT3")+
    scale_fill_manual(values=c("IMPACT2"="#a6cee3","IMPACT3"="#6a3d9a"))+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


## Now we compare the incoming and outgoing signaling molecules in myeloid spots across different IMPACT archetypes
D2_Visium_MP5_gene_1a_1b <- readRDS("../data/D2_Visium_MP5_gene_1a_1b.rds")
genes_highlight_MP5_1a_1b <- c("SLITRK","OPIOID","VEGI","CD30","CD70","FLT3","IL4","CD200","OX40","PDL2","FASLG","CD86","TNF","BAFF","PD-L1","SELL","SIRP","PVR")

pdf("molecule_plot_MP5_gene_1a_1b.pdf",width=7,height=5)
ggplot(subset(D2_Visium_MP5_gene_1a_1b, specificity != "Shared"),
        aes(x=outgoing, y=incoming, color = specificity, shape = specificity.out.in, label = labels))+
    geom_point(size=1.5)+
    xlab("Differential outgoing interaction strength")+
    ylab("Differential incoming interaction strength")+
    xlim(c(-0.065,0.015))+
    ylim(c(-0.025,0.01))+
    ggtitle("MP5 IMPACT-1a vs IMPACT-1b")+
    geom_text_repel(data = subset(D2_Visium_MP5_gene_1a_1b, labels %in% genes_highlight_MP5_1a_1b),
                aes(color = specificity), size = 3, point.padding = 0.3, box.padding = 0.7, 
                min.segment.length = 0, segment.size = 0.3, max.overlaps = 500)+
    scale_color_manual(values=c("IMPACT1a specific"="#cc4c02", "IMPACT1b specific"="#006837", 
                                "IMPACT2 specific"="#0570b0", "IMPACT3 specific"="#6a51a3"))+
    theme_classic()
dev.off()
# similar for other comparisons, e.g., MP5 IMPACT-1a vs IMPACT-2, MP7 IMPACT-1a vs IMPACT-2, etc.


## We also specifically focusing on two representative pathways, i.e., IL1 and SPP1
color_for_plot <- c(colorRampPalette(c('white',"#ffffcc"))(10),
                    colorRampPalette(c('#ffffcc',"#ffeda0"))(10),
                    colorRampPalette(c('#ffeda0',"#fed976"))(10),
                    colorRampPalette(c('#fed976',"#feb24c"))(10),
                    colorRampPalette(c('#feb24c',"#fd8d3c"))(10),
                    colorRampPalette(c('#fd8d3c',"#fc4e2a"))(10),
                    colorRampPalette(c('#fc4e2a',"#e31a1c"))(10),
                    colorRampPalette(c('#e31a1c',"#bd0026"))(10),
                    colorRampPalette(c('#bd0026',"#800026"))(10))

# This is the communication probability of IL1 signaling across different IMPACT archetypes in the tumor region.
D2_Visium_IL1_commun_prob_IMPACT <- readRDS("../data/D2_Visium_IL1_commun_prob_IMPACT.rds")

pdf("IL1_IMPACT.pdf",width=5,height=7)
pheatmap::pheatmap(D2_Visium_IL1_commun_prob_IMPACT,
                    cluster_rows = F,cluster_cols = F,
                    show_row_names = T, show_column_names = T,
                    color = color_for_plot,
                    main = "IL1 signaling",
                    cellwidth = 15, cellheight = 15, border_color = "white",
                    fontsize = 10, fontsize_row = 8, fontsize_col = 8, angle_col = 45
                   )
dev.off()

# This is the communication probability of SPP1 signaling across different IMPACT archetypes in the tumor region.
D2_Visium_SPP1_commun_prob_IMPACT <- readRDS("../data/D2_Visium_SPP1_commun_prob_IMPACT.rds")

pdf("SPP1_IMPACT.pdf",width=5,height=7)
pheatmap::pheatmap(D2_Visium_SPP1_commun_prob_IMPACT,
                    cluster_rows = F,cluster_cols = F,
                    show_row_names = T, show_column_names = T,
                    color = color_for_plot,
                    main = "SPP1 signaling",
                    cellwidth = 15, cellheight = 15, border_color = "white",
                    fontsize = 10, fontsize_row = 8, fontsize_col = 8, angle_col = 45
                   )
dev.off()

# similar for plotting the communication probability of IL1 and SPP1 signaling across different lesion types in the tumor region.