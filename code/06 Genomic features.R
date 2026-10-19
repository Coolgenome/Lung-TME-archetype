### Genomic features of different IMPACT archetypes ###
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
library(grid)

### Genomic features of different IMPACT archetypes in the study cohort ###
# set up metadata for the study cohort
# IMPACT archetypes of each sample
IMPACT1a <- c("P2_AAH","P9_AIS","P20_AAH","P22_AAH","P4_AAH2","P6_AAH","P9_AAH","P18_LUAD","P22_AIS","P11_LUAD")
IMPACT1b <- c("P7_LUAD1","P15_MIA","P19_LUAD","P17_AIS","P5_AIS","P10_MIA","P11_AAH","P16_AIS","P19_AIS","P1_AAH")
IMPACT2  <- c("P3_LUAD","P22_LUAD","P21_AIS2","P9_LUAD","P17_LUAD","P6_LUAD","P12_AIS","P13_LUAD","P10_LUAD","P24_LUAD",
              "P8_AIS","P7_LUAD2","P15_LUAD","P20_LUAD","P13_MIA")
IMPACT3  <- c("P21_AIS1","P23_AIS1","P14_AIS","P25_AAH","P24_AAH","P18_MIA","P14_LUAD","P16_LUAD","P25_LUAD","P4_LUAD",
              "P12_LUAD","P23_AIS2","P21_LUAD","P8_LUAD","P3_AIS","P23_LUAD","P2_LUAD","P1_LUAD")

Impact_group <- c(
  setNames(rep("IMPACT-1a", length(IMPACT1a)), IMPACT1a),
  setNames(rep("IMPACT-1b", length(IMPACT1b)), IMPACT1b),
  setNames(rep("IMPACT-2",  length(IMPACT2)),  IMPACT2),
  setNames(rep("IMPACT-3",  length(IMPACT3)),  IMPACT3)
)

Impact_group <- factor(Impact_group, levels = c("IMPACT-1a", "IMPACT-1b", "IMPACT-2", "IMPACT-3"))

# Lesion type of each sample
lesion_type <- gsub("^P[0-9]+_", "", names(Impact_group))
lesion_type <- gsub("[0-9]+$", "", lesion_type) 
names(lesion_type) <- names(Impact_group)
lesion_type <- factor(lesion_type, levels = c("AAH", "AIS", "MIA", "LUAD"))

# Mutation count
Mutation_count <- readRDS("../code/WES_mutation_count_matrix.rds")

## order samples within each IMPACT archetype by mutation count
total_count <- rowSums(Mutation_count)
sample_order <- unlist(lapply(levels(Impact_group), function(g) {
  s <- names(Impact_group)[Impact_group == g]
  s[order(total_count[s], decreasing = TRUE)]
}))

Impact_group <- Impact_group[sample_order]
lesion_type <- lesion_type[sample_order]
Mutation_count <- Mutation_count[sample_order, , drop = FALSE]

# Read in driver mutation matrix
WES_driver_mutation <- readRDS("../data/WES_driver_mutation.rds")

# set up colors for oncoplot
mut_col <- c(Truncating = "#E8383A", Missense_Mutation = "#2C4A9E", In_Frame_InDel = "#5FBB68")
mut_lab <- c(Truncating = "Truncating", Missense_Mutation = "Missense", In_Frame_InDel = "In frame InDel")
lesion_col <- c(AAH = "#E8437A", AIS = "#4CAF50", MIA = "#F5D530", LUAD = "#4B4FA8")
group_col <- c("IMPACT-1a" = "#F07C1F", "IMPACT-1b" = "#A7D48A", "IMPACT-2" = "#9FCDE3", "IMPACT-3" = "#6A3D9A")

alter_fun <- list(
  background = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = "#D9D9D9", col = NA)),
  Truncating = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = mut_col["Truncating"], col = NA)),
  Missense_Mutation = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = mut_col["Missense_Mutation"], col = NA)),
  In_Frame_InDel = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = mut_col["In_Frame_InDel"], col = NA))
)

# set up top annotation
top_ann <- HeatmapAnnotation(
  "Lesion type" = lesion_type,
  "Mutation count" = anno_barplot(count_mat,
    gp = gpar(fill = mut_col[colnames(count_mat)], col = NA),
    border = FALSE, height = unit(2.5, "cm")
  ),
  col = list("Lesion type" = lesion_col),
  annotation_name_side = "left",
  gap = unit(1, "mm"),
  show_legend = FALSE
)

## set up right annotation
right_ann <- rowAnnotation(
  "No. sample" = anno_oncoprint_barplot(
    mut_types, border = FALSE,
    width = unit(2.5, "cm")
  )
)

# plot
oncoplot <- oncoPrint(
  WES_driver_mutation,
  alter_fun = alter_fun,
  col = mut_col,
  alter_fun_is_vectorized = FALSE,
  top_annotation = top_ann,
  right_annotation = right_ann,
  row_order = drivers,
  column_order = samples,
  column_split = Impact_group,
  column_title = levels(Impact_group),
  column_title_gp = gpar(fill = group_col[levels(Impact_group)], fontsize = 12, fontface = "bold", col = "white"),
  column_gap = unit(2, "mm"),
  show_column_names = FALSE,
  show_pct = TRUE,
  pct_side = "right",
  pct_gp = gpar(fontsize = 10),
  show_heatmap_legend = FALSE,
  row_names_side = "left"
)

lgd_lesion <- Legend(title = "Lesion type", at = names(lesion_col), legend_gp = gpar(fill = lesion_col))
lgd_mut <- Legend(title = "Mutation type", at = names(mut_lab), labels = mut_lab, legend_gp = gpar(fill = mut_col))

pdf("Driver_mutation_oncoplot.pdf",width=15,height=5)
draw(oncoplot, annotation_legend_list = list(lgd_lesion, lgd_mut), 
     annotation_legend_side = "bottom", heatmap_legend_side = "bottom", 
     merge_legend = TRUE)
dev.off()



### Tumor genomic metrics: TMB, sMF, FGA ###
WES_genomic_metric <- readRDS("../data/WES_genomic_metric.rds")

# here we plot mutation count
pdf("Mutation_count.pdf",width=3.5,height=3)
ggplot(WES_genomic_metric, aes(IMPACT, mutation_count, color=IMPACT))+
        geom_boxplot(outlier.alpha=0)+
        geom_point(size=1.5, position="jitter")+
        theme_classic()+
        scale_color_manual(values=c("#ff7f00","#b2df8a","#a6cee3","#6a3d9a","#e31a1c","#253494"))+
        theme(axis.ticks.length = unit(0.1,'cm'),legend.position = "right", axis.text.x=element_text(angle=45,hjust=1))+
        guides(fill=guide_legend(title = NULL))
dev.off()
# similar for sMF and FGA


### Tumor evolutionary patterns ###
Tumor_evolutionary_pattern <- readRDS("../data/Tumor_evolutionary_pattern.rds")

pattern_levels <- c("Pattern1a", "Pattern1b", "Pattern2", "Other")
pattern_col <- c(Pattern1a = "#8FC4E3", Pattern1b = "#5BA85A",
                  Pattern2 = "#E24A33", Other = "#CFCFCF")

group_levels <- c("IMPACT1a", "IMPACT1b", "IMPACT2", "IMPACT3")
group_col <- c(IMPACT1a = "#F07C1F", IMPACT1b = "#5CA344",
                IMPACT2 = "#3E7CB1", IMPACT3 = "#6A3D9A")

make_pie <- function(group_id) {
  d <- Tumor_evolutionary_pattern %>% filter(IMPACT == group_id)
  cellnum <- as.data.frame(table(factor(d$Evolutionary_pattern)))
  colnames(cellnum) <- c("evolution_pattern", "Num")

  ggplot(cellnum, aes("", Num, fill = evolution_pattern)) +
    geom_bar(stat = "identity", position = "fill") +
    coord_polar(theta = "y") +
    ggtitle(group_id) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13, color = group_col[group_id]),
          legend.position = "none") +
    scale_fill_manual(values = pattern_col, drop = FALSE)
}

plots <- lapply(group_levels, make_pie)
combined <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]])
pdf("Tumor_evolutionary patterns.pdf",width=7,height=5)
print(combined)
dev.off()

