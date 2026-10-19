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
V3_inhouse_WSI_prediction <- readRDS("../data/V3_inhouse_WSI_prediction.rds")

### sample distribution ###
df <- as.data.frame(table(V3_inhouse_WSI_prediction$lesion_type))
df$Var1 <- factor(df$Var1, levels=c("AAH","AIS","MIA","ADC"))

pdf("inhouse_chort_sample_compos.pdf", width = 4, height = 4)
ggplot(df, aes(x = "", y = Freq, fill = Var1)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  labs(title = "Sample Distribution") +
  scale_fill_manual(values = c("AAH"="#e6194b","AIS"="#3cb44b","MIA"="#ffe119","ADC"="#4363d8")) +
  theme_void()
dev.off()


### LUAD IMPACT distribution ###
V3_inhouse_WSI_prediction_LUAD <- filter(V3_inhouse_WSI_prediction, lesion_type=="ADC")
df <- as.data.frame(table(V3_inhouse_WSI_prediction_LUAD$IMPACT))

pdf("LUAD_sample_compos.pdf", width = 4, height = 4)
ggplot(df, aes(x = "", y = Freq, fill = Var1)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  labs(title = "Sample Distribution") +
  scale_fill_manual(values = c("IMPACT1a/2"="#d73027","IMPACT1b"="#b2df8a","IMPACT3"="#6a3d9a","not_assessable"="#d9d9d9")) +
  theme_void()
dev.off()


### Ro/e analysis of key TME compartment distributions across lesion types ###
cell_res <- V3_inhouse_WSI_prediction %>% group_by(lesion_type) %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm = TRUE)), .groups = "drop")

cell_res <- cell_res[,c(1:4)] %>% as.data.frame() %>% column_to_rownames(var = "lesion_type")
res <- chisq.test(cell_res)
ratio <- res$observed/res$expected

library(ComplexHeatmap)
library(circlize)
col_min = min(ratio); col_max = max(ratio)
col_fun = colorRamp2(c(col_min, 1, col_max), c("#2166ac", "white", "#b2182b"))

ratio <- ratio[c("AAH","AIS","MIA","ADC"),c("Myeloid_in_Tumor","Lymph_in_Tumor","Fibrous_in_Tumor")]

pdf("Roe.pdf",width=5,height=5)
Heatmap(ratio,  name = "Ro/e", col = col_fun,
        cluster_rows = F, cluster_columns = F)
dev.off()

### comparisons of TME compositions across lesion types ###
V3_inhouse_WSI_prediction$TME_total <- rowSums(V3_inhouse_WSI_prediction[,c(3:6)])
V3_inhouse_WSI_prediction$Mac_TME_ratio <- V3_inhouse_WSI_prediction$Myeloid_in_Tumor / V3_inhouse_WSI_prediction$TME_total
V3_inhouse_WSI_prediction$Lymph_TME_ratio <- V3_inhouse_WSI_prediction$Lymph_in_Tumor / V3_inhouse_WSI_prediction$TME_total
V3_inhouse_WSI_prediction$Fib_TME_ratio <- V3_inhouse_WSI_prediction$Fibrous_in_Tumor / V3_inhouse_WSI_prediction$TME_total
V3_inhouse_WSI_prediction$lesion_type <- factor(V3_inhouse_WSI_prediction$lesion_type, levels = c("AAH","AIS","MIA","ADC"))


pdf("Mac_TME_ratio_lesion_type.pdf", width = 3.5, height = 4)
ggplot(V3_inhouse_WSI_prediction, aes(x = lesion_type, y = Mac_TME_ratio, color = lesion_type)) +
  geom_boxplot(width=0.6, outlier.shape = NA) +
  geom_jitter(width=0.1, alpha = 0.8, size=1.5) +
  theme_classic() +
  xlab("Lesion Type") +
  ylab("Prop in TME") +
  ggtitle("Mac_TME_ratio") +
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("AAH"="#e6194b","AIS"="#3cb44b","MIA"="#ffe119","ADC"="#4363d8"))
dev.off()

library(FSA)
kruskal.test(Mac_TME_ratio ~ lesion_type, data = V3_inhouse_WSI_prediction)
dunnTest(Mac_TME_ratio ~ lesion_type, data = V3_inhouse_WSI_prediction)


### distribution of IMPACT archetypes across lesion types ###
V3_inhouse_WSI_prediction_tumor_assessable <- filter(V3_inhouse_WSI_prediction, lesion_type %in% c("AIS","MIA","ADC") & IMPACT!="not_assessable")
V3_inhouse_WSI_prediction_tumor_assessable$lesion_type <- factor(V3_inhouse_WSI_prediction_tumor_assessable$lesion_type, levels=c("AIS","MIA","ADC"))

lesion_IMPACT <- table(V3_inhouse_WSI_prediction_tumor_assessable$lesion_type, V3_inhouse_WSI_prediction_tumor_assessable$IMPACT)
res <- chisq.test(lesion_IMPACT)
ratio <- res$observed/res$expected

col_min = min(ratio); col_max = max(ratio)
col_fun = colorRamp2(c(col_min, 1, col_max), c("#2166ac", "white", "#b2182b"))

pdf("Roe.pdf",width=5,height=5)
Heatmap(ratio,  name = "Ro/e", col = col_fun,
        cluster_rows = F, cluster_columns = T,
        cell_fun = function(j, i, x, y, width, height, fill) {
          val = ratio[i, j]
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



