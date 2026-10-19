### HE-based spatial archetype and survival analysis in TCGA ###

### load packages ###
library(dplyr)
library(tidyverse)
library(data.table)
library(ggplot2)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(grid)

V4_TCGA_prediction <- readRDS("../data/V4_TCGA_prediction.rds")

### calculate slide compositions
df <- as.data.frame(table(V4_TCGA_prediction$IMPACT_slide))

pdf("slide_IMPACT_compos.pdf", width = 4, height = 4)
ggplot(df, aes(x = "", y = Freq, fill = Var1)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  labs(title = "Slide Distribution") +
  scale_fill_manual(values = c("IMPACT1/2"="#d73027","IMPACT3"="#6a3d9a","Not_accessible"="#d9d9d9")) +
  theme_void()
dev.off()


### compare tumor and TME components across predicted IMPACT archetypes
V4_TCGA_prediction_filter <- filter(V4_TCGA_prediction, IMPACT_slide != "Not_accessible")
V4_TCGA_prediction_filter$TME_total <- V4_TCGA_prediction_filter$Fibrous_prop_Tumor + V4_TCGA_prediction_filter$Lymph_prop_Tumor + 
                                       V4_TCGA_prediction_filter$Myeloid_prop_Tumor + V4_TCGA_prediction_filter$Vessel_prop_Tumor
V4_TCGA_prediction_filter$Lymph_prop_TME <- V4_TCGA_prediction_filter$Lymph_prop_Tumor / V4_TCGA_prediction_filter$TME_total
V4_TCGA_prediction_filter$Myeloid_prop_TME <- V4_TCGA_prediction_filter$Myeloid_prop_Tumor / V4_TCGA_prediction_filter$TME_total
V4_TCGA_prediction_filter$Fibrous_prop_TME <- V4_TCGA_prediction_filter$Fibrous_prop_Tumor / V4_TCGA_prediction_filter$TME_total

pdf("TCGA_pathol_Tumor_prop.pdf",width=4,height=2.5)
ggplot(V4_TCGA_prediction_filter, aes(x = IMPACT_slide, y = Tumor_prop_Tumor, color=IMPACT_slide))+
  geom_boxplot(width=0.6, outliers = F)+
  geom_jitter(width=0.1, size=0.5, alpha=0.3)+
  scale_color_manual(values=c("IMPACT1/2"="#d73027","IMPACT3"="#6a3d9a"))+
  theme_classic()
dev.off()

pdf("TCGA_pathol_Lymph_prop_TME.pdf",width=4,height=2.5)
ggplot(V4_TCGA_prediction_filter, aes(x = IMPACT_slide, y = Lymph_prop_TME, color=IMPACT_slide))+
  geom_boxplot(width=0.6, outliers = F)+
  geom_jitter(width=0.1, size=0.5, alpha=0.3)+
  scale_color_manual(values=c("IMPACT1/2"="#d73027","IMPACT3"="#6a3d9a"))+
  theme_classic()
dev.off()

p1 = wilcox.test(Tumor_prop_Tumor ~ IMPACT_slide, data = V4_TCGA_prediction_filter)$p.value
p2 = wilcox.test(Lymph_prop_TME ~ IMPACT_slide, data = V4_TCGA_prediction_filter)$p.value
p3 = wilcox.test(Myeloid_prop_TME ~ IMPACT_slide, data = V4_TCGA_prediction_filter)$p.value
p4 = wilcox.test(Fibrous_prop_TME ~ IMPACT_slide, data = V4_TCGA_prediction_filter)$p.value
p = c(p1, p2, p3, p4)
p.adj = p.adjust(p, method = "BH")
p.adj


### examine IMPACT signature score
# Validate IMPACT signatures using pseudo-bulk profiles from snRNA-seq data ###
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

# expression of IMPACT signature in TCGA cohort
V4_TCGA_IMPACT_signature <- readRDS("../data/V4_TCGA_IMPACT_signature.rds")
IMPACT_final <- select(V4_TCGA_prediction_filter, c("SAMPLE_ID","IMPACT_patient")) %>% unique()
V4_TCGA_IMPACT_signature <- left_join(V4_TCGA_IMPACT_signature, IMPACT_final, by="SAMPLE_ID")

pdf("TCGA_IMPACT2_sig.pdf",width=4,height=2.5)
ggplot(filter(V4_TCGA_IMPACT_signature, !is.na(V4_TCGA_IMPACT_signature$IMPACT_patient)), 
       aes(x=IMPACT_patient, y=IMPACT2_sig, color=IMPACT_patient))+
  geom_boxplot(width=0.6, outliers = F)+
  geom_jitter(width=0.1, size=0.5, alpha=0.3)+
  scale_color_manual(values = c("IMPACT1/2"="#d73027","IMPACT3"="#6a3d9a"))+
  theme_classic()
dev.off()

pdf("TCGA_IMPACT3_sig.pdf",width=4,height=2.5)
ggplot(filter(V4_TCGA_IMPACT_signature, !is.na(V4_TCGA_IMPACT_signature$IMPACT_patient)), 
       aes(x=IMPACT_patient, y=IMPACT3_sig, color=IMPACT_patient))+
  geom_boxplot(width=0.6, outliers = F)+
  geom_jitter(width=0.1, size=0.5, alpha=0.3)+
  scale_color_manual(values = c("IMPACT1/2"="#d73027","IMPACT3"="#6a3d9a"))+
  theme_classic()
dev.off()

# due to the very representation of IMPACT-1a/1b cases in LUAD samples, we don't compare the IMPACT-1a/1b signatures here

p1 = wilcox.test(IMPACT2 ~ IMPACT_patient, data = filter(V4_TCGA_IMPACT_signature, !is.na(V4_TCGA_IMPACT_signature$IMPACT_patient)))$p.value
p2 = wilcox.test(IMPACT3 ~ IMPACT_patient, data = filter(V4_TCGA_IMPACT_signature, !is.na(V4_TCGA_IMPACT_signature$IMPACT_patient)))$p.value
p = c(p1, p2)
p.adj = p.adjust(p, method = "BH")
p.adj


### Genomic features of TCGA cohort
V4_TCGA_driver_mutation <- readRDS("../data/V4_TCGA_driver_mutation.rds")
IMPACT12_sample <- intersect(IMPACT_final$SAMPLE_ID[which(IMPACT_final$IMPACT_patient == "IMPACT1/2")], rownames(V4_TCGA_driver_mutation))
IMPACT3_sample <- intersect(IMPACT_final$SAMPLE_ID[which(IMPACT_final$IMPACT_patient == "IMPACT3")], rownames(V4_TCGA_driver_mutation))

# Specify alteration functions
alter_fun = list(background = function(x, y, w, h) grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = "#f0f0f0", col = NA)),
                 Nonsense_Mutation = function(x, y, w, h) grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = "red", col = NA)),
                 Missense_Mutation = function(x, y, w, h) grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = "blue", col = NA)),
                 In_Frame_Indel = function(x, y, w, h) grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = "green", col = NA)),
                 Frame_Shift_Indel = function(x, y, w, h) grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = "#b2df8a", col = NA)),
                 Splice_Region = function(x, y, w, h) grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = "#fdbf6f", col = NA)),
                 Splice_Site = function(x, y, w, h) grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = "#cab2d6", col = NA)),
                 Multi_Hit = function(x, y, w, h) grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = "#fb9a99", col = NA)))

# Set color legend
col = c("Nonsense_Mutation" = "red", "Missense_Mutation" = "blue", "In_Frame_Indel" = "green", 
        "Frame_Shift_Indel"="#b2df8a", "Splice_Region"="#fdbf6f", "Splice_Site"="#cab2d6",
        "Multi_Hit"="#fb9a99")

# prepare data for plot
V4_TCGA_driver_mutation_plot <- t(V4_TCGA_driver_mutation)[,c(IMPACT12_sample,IMPACT3_sample)]

# prepare annotation
IMPACT = data.frame(SAMPLE_ID = colnames(V4_TCGA_driver_mutation_plot),
                    IMPACT_assignment = c(rep("IMPACT12",162), rep("IMPACT3",263)))

top_anno = HeatmapAnnotation(IMPACT=anno_simple(IMPACT$IMPACT_assignment, col = c("IMPACT12" = "#fc4e2a", "IMPACT3" = "#6a3d9a")),
                             annotation_name_side = "right")

pdf("TCGA_oncoplot.pdf",width=15,height=5)
oncoPrint(V4_TCGA_driver_mutation_plot, 
          alter_fun = alter_fun, 
          col = col,
          remove_empty_columns = FALSE, 
          remove_empty_rows = FALSE,
          show_column_names = FALSE, 
          show_row_names = TRUE,
          row_names_side = "left",
          show_pct = TRUE,
          pct_side = "right",
          pct_gp = gpar(fontsize = 10),
          top_annotation = top_anno,
          column_title = "Mutation Landscape")
dev.off()

# Genomic features (TMB information was from TCGA metadata, downloaded from cBioportal)
V4_TCGA_genomic_metric <- readRDS("../data/V4_TCGA_genomic_metric.rds")
V4_TCGA_genomic_metric <- left_join(V4_TCGA_genomic_metric, IMPACT_final, by="SAMPLE_ID")

pdf("TCGA_TMB.pdf",width=3,height=3)
ggplot(filter(V4_TCGA_genomic_metric, !is.na(V4_TCGA_genomic_metric$IMPACT_patient)),
       aes(x=IMPACT_patient, y=TMB_NONSYNONYMOUS, color=IMPACT_patient))+
  geom_boxplot(width=0.6,outlier.alpha=0)+
  geom_jitter(width=0.1, size=1, alpha=0.5)+
  scale_color_manual(values=c("#e31a1c","#6a3d9a","#253494"))+
  theme_classic()+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))+
  theme(axis.ticks.length = unit(0.1,'cm'), legend.position = "none")
dev.off()

wilcox.test(TMB_NONSYNONYMOUS ~ IMPACT_patient, data = filter(V4_TCGA_genomic_metric, !is.na(V4_TCGA_genomic_metric$IMPACT_patient)))

pdf("TCGA_FGA.pdf",width=3,height=3)
ggplot(filter(V4_TCGA_genomic_metric, !is.na(V4_TCGA_genomic_metric$IMPACT_patient)),
       aes(x=IMPACT_patient, y=FGA, color=IMPACT_patient))+
  geom_boxplot(width=0.6,outlier.alpha=0)+
  geom_jitter(width=0.1, size=1, alpha=0.5)+
  scale_color_manual(values=c("#e31a1c","#6a3d9a","#253494"))+
  theme_classic()+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))+
  theme(axis.ticks.length = unit(0.1,'cm'), legend.position = "none")
dev.off()

wilcox.test(FGA ~ IMPACT_patient, data = filter(V4_TCGA_genomic_metric, !is.na(V4_TCGA_genomic_metric$IMPACT_patient)))


### Survival analysis of TCGA cohort
library(survival)
library(survminer)

V4_TCGA_sample_summary <- V4_TCGA_prediction %>% group_by(SAMPLE_ID) %>% 
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>% as.data.frame()

V4_TCGA_sample_summary <- left_join(V4_TCGA_sample_summary, IMPACT_final, by="SAMPLE_ID")

# TCGA metadata, downloaded from cBioportal
TCGA_metadata_sample <- read.table("../luad_tcga_pan_can_atlas_2018/data_clinical_sample.txt",sep="\t",header=T, quote = "")
TCGA_metadata_pt <- read.table("../luad_tcga_pan_can_atlas_2018/data_clinical_patient.txt",sep="\t",header=T)

V4_TCGA_sample_summary <- left_join(V4_TCGA_sample_summary, TCGA_metadata_sample, by="SAMPLE_ID")
V4_TCGA_sample_summary <- left_join(V4_TCGA_sample_summary, TCGA_metadata_pt, by="PATIENT_ID")

# We cap the survival data to 5 years
V4_TCGA_sample_summary$OS_status <- NA
V4_TCGA_sample_summary$OS_status[V4_TCGA_sample_summary$OS_STATUS == "0:LIVING"] <- 0
V4_TCGA_sample_summary$OS_status[V4_TCGA_sample_summary$OS_STATUS == "1:DECEASED"] <- 1

V4_TCGA_sample_summary$OS_mo_cap_5y <- V4_TCGA_sample_summary$OS_MONTHS
V4_TCGA_sample_summary$OS_mo_cap_5y[which(V4_TCGA_sample_summary$OS_mo_cap_5y >= 60)] <- 60
V4_TCGA_sample_summary$OS_status_cap_5y <- V4_TCGA_sample_summary$OS_status
V4_TCGA_sample_summary$OS_status_cap_5y[which(V4_TCGA_sample_summary$OS_MONTHS > 60 & V4_TCGA_sample_summary$OS_status == 1)] <- 0

surv_object <- Surv(time = V4_TCGA_sample_summary$OS_mo_cap_5y, event = V4_TCGA_sample_summary$OS_status_cap_5y)

# dichotomize intratumoral density of TME cells and perform survival analysis
Lymph_median <- median(V4_TCGA_sample_summary$Lymph_prop_Tumor)
V4_TCGA_sample_summary$Lymphocytes_dichotomize <- "Low"
V4_TCGA_sample_summary$Lymphocytes_dichotomize[which(V4_TCGA_sample_summary$Lymph_prop_Tumor >= Lymph_median)] <- "High"
surv_fit_lymph <- survfit(surv_object ~ Lymphocytes_dichotomize, data = V4_TCGA_sample_summary)

pdf("TCGA_OS_lymph.pdf",width=7,height=5)
ggsurvplot(surv_fit_lymph, 
           data = V4_TCGA_sample_summary,
           legend.labs = levels(V4_TCGA_sample_summary$Lymphocytes_dichotomize), 
           palette=c("#d73027","#4393c3"),
           pval = TRUE, 
           conf.int = FALSE,
           surv.median.line = "hv",
           risk.table = TRUE,
           ggtheme = theme_classic())
dev.off()

survdiff(formula = surv_object ~ Lymphocytes_dichotomize, data = V4_TCGA_sample_summary)


Myeloid_median <- median(V4_TCGA_sample_summary$Myeloid_prop_Tumor)
V4_TCGA_sample_summary$Myeloid_dichotomize <- "Low"
V4_TCGA_sample_summary$Myeloid_dichotomize[which(V4_TCGA_sample_summary$Myeloid_prop_Tumor >= Myeloid_median)] <- "High"
surv_fit_myeloid <- survfit(surv_object ~ Myeloid_dichotomize, data = V4_TCGA_sample_summary)

pdf("TCGA_OS_myeloid.pdf",width=7,height=5)
ggsurvplot(surv_fit_myeloid, 
           data = V4_TCGA_sample_summary,
           legend.labs = levels(V4_TCGA_sample_summary$Myeloid_dichotomize), 
           palette=c("#d73027","#4393c3"),
           pval = TRUE, 
           conf.int = FALSE,
           surv.median.line = "hv",
           risk.table = TRUE,
           ggtheme = theme_classic())
dev.off()

survdiff(formula = surv_object ~ Myeloid_dichotomize, data = V4_TCGA_sample_summary)


Fibrous_median <- median(V4_TCGA_sample_summary$Fibrous_prop_Tumor)
V4_TCGA_sample_summary$Fibrous_dichotomize <- "Low"
V4_TCGA_sample_summary$Fibrous_dichotomize[which(V4_TCGA_sample_summary$Fibrous_prop_Tumor >= Fibrous_median)] <- "High"
surv_fit_fibrous <- survfit(surv_object ~ Fibrous_dichotomize, data = V4_TCGA_sample_summary)

pdf("TCGA_OS_fibrous.pdf",width=7,height=5)
ggsurvplot(surv_fit_fibrous, 
           data = V4_TCGA_sample_summary,
           legend.labs = levels(V4_TCGA_sample_summary$Fibrous_dichotomize), 
           palette=c("#d73027","#4393c3"),
           pval = TRUE, 
           conf.int = FALSE,
           surv.median.line = "hv",
           risk.table = TRUE,
           ggtheme = theme_classic())
dev.off()

survdiff(formula = surv_object ~ Fibrous_dichotomize, data = V4_TCGA_sample_summary)


# survival analysis based on predicted IMPACT archetype
surv_fit_IMPACT <- survfit(surv_object ~ IMPACT_patient, data = V4_TCGA_sample_summary)

pdf("TCGA_OS_fibrous.pdf",width=7,height=5)
ggsurvplot(surv_fit_IMPACT, 
           data = V4_TCGA_sample_summary,
           legend.labs = levels(V4_TCGA_sample_summary$IMPACT_patient), 
           palette=c("#d73027","#6a3d9a"),
           pval = TRUE, 
           conf.int = FALSE,
           surv.median.line = "hv",
           risk.table = TRUE,
           ggtheme = theme_classic())
dev.off()

survdiff(formula = surv_object ~ IMPACT_patient, data = V4_TCGA_sample_summary)

# clean up disease stage information
V4_TCGA_sample_summary$Stage <- NA
V4_TCGA_sample_summary$Stage[which(V4_TCGA_sample_summary$AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE I","STAGE IA","STAGE IB"))] <- "Stage I"
V4_TCGA_sample_summary$Stage[which(V4_TCGA_sample_summary$AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE II","STAGE IIA","STAGE IIB"))] <- "Stage II"
V4_TCGA_sample_summary$Stage[which(V4_TCGA_sample_summary$AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE IIIA","STAGE IIIB"))] <- "Stage III"
V4_TCGA_sample_summary$Stage[which(V4_TCGA_sample_summary$AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE IV"))] <- "Stage IV"
V4_TCGA_sample_summary$Stage <- factor(V4_TCGA_sample_summary$Stage, levels=c("Stage I","Stage II","Stage III","Stage IV"))

# model TME compartments together with clinical factors
cox_model_TME <- coxph(Surv(OS_mo_cap_5y, OS_status_cap_5y) ~ Lymph_prop_Tumor + Myeloid_prop_Tumor + Fibrous_prop_Tumor + AGE + SEX + Stage, data = V4_TCGA_sample_summary)
summary(cox_model_TME)

pdf("TCGA_cox_model_TME_clinical.pdf",width=10,height=6)
ggforest(
  model = cox_model_TME,
  data = V4_TCGA_sample_summary,
  main = "Hazard ratio",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 1.0,
  refLabel = "Reference",
  noDigits = 2
)
dev.off()

# model IMPACT together with clinical factors
cox_model_IMPACT <- coxph(Surv(OS_mo_cap_5y, OS_status_cap_5y) ~ IMPACT_patient + AGE + SEX + Stage, data = V4_TCGA_sample_summary)
summary(cox_model_IMPACT)

pdf("TCGA_cox_model_IMPACT_clinical.pdf",width=10,height=6)
ggforest(
  model = cox_model_IMPACT,
  data = V4_TCGA_sample_summary,
  main = "Hazard ratio",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 1.0,
  refLabel = "Reference",
  noDigits = 2
)
dev.off()

# add IMPACT signature information
V4_TCGA_IMPACT_signature <- readRDS("../data/V4_TCGA_IMPACT_signature.rds")
V4_TCGA_sample_summary <- left_join(V4_TCGA_sample_summary, V4_TCGA_IMPACT_signature, by="SAMPLE_ID")

cox_model_IMPACT_sig <- coxph(Surv(OS_mo_cap_5y, OS_status_cap_5y) ~ IMPACT1a_sig + IMPACT1b_sig + IMPACT2_sig + IMPACT3_sig, data = V4_TCGA_sample_summary)
summary(cox_model_IMPACT_sig)

# now based on the coef from above fitting, we derive a composite score by weighted sum of IMPACT signatures
IMPACT_sig_coeff <- cox_model_IMPACT_sig$coefficients
V4_TCGA_sample_summary$IMPACT_sig_composite <- IMPACT_sig_coeff[1] * V4_TCGA_sample_summary$IMPACT1a_sig +
                                               IMPACT_sig_coeff[2] * V4_TCGA_sample_summary$IMPACT1b_sig +
                                               IMPACT_sig_coeff[3] * V4_TCGA_sample_summary$IMPACT2_sig +
                                               IMPACT_sig_coeff[4] * V4_TCGA_sample_summary$IMPACT3_sig
# dichotomize IMPACT signature composite score
IMPACT_sig_composite_median <- median(V4_TCGA_sample_summary$IMPACT_sig_composite, na.rm=TRUE)
V4_TCGA_sample_summary$IMPACT_sig_composite_binary <- NA
V4_TCGA_sample_summary$IMPACT_sig_composite_binary[which(V4_TCGA_sample_summary$IMPACT_sig_composite < IMPACT_sig_composite_median)] <- "Low"
V4_TCGA_sample_summary$IMPACT_sig_composite_binary[which(V4_TCGA_sample_summary$IMPACT_sig_composite >= IMPACT_sig_composite_median)] <- "High"
# derive IMPACT risk profile from combined IMPACT archetype and IMPACT signature composite score
V4_TCGA_sample_summary$IMPACT_risk <- NA
V4_TCGA_sample_summary$IMPACT_risk[!is.na(V4_TCGA_sample_summary$IMPACT_patient) & !is.na(V4_TCGA_sample_summary$IMPACT_sig_composite_binary)] <- "Medium_risk"
V4_TCGA_sample_summary$IMPACT_risk[which(V4_TCGA_sample_summary$IMPACT_patient=="IMPACT1/2" & V4_TCGA_sample_summary$IMPACT_sig_composite_binary == "Low")] <- "Low_risk"
V4_TCGA_sample_summary$IMPACT_risk[which(V4_TCGA_sample_summary$IMPACT_patient=="IMPACT3" & V4_TCGA_sample_summary$IMPACT_sig_composite_binary == "High")] <- "High_risk"
V4_TCGA_sample_summary$IMPACT_risk <- factor(V4_TCGA_sample_summary$IMPACT_risk, levels=c("High_risk","Medium_risk","Low_risk"))

surv_fit_IMPACT_risk <- survfit(surv_object ~ IMPACT_risk, data = V4_TCGA_sample_summary)

pdf("TCGA_survival_IMPACT_risk.pdf",width=7,height=5)
ggsurvplot(surv_fit_IMPACT_risk, 
           data = V4_TCGA_sample_summary,
           legend.labs = levels(V4_TCGA_sample_summary$IMPACT_risk), 
           palette=c("#e31a1c","#fed976","#01665e"),
           pval = TRUE, 
           conf.int = FALSE,
           surv.median.line = "hv",
           risk.table = TRUE,
           ggtheme = theme_classic())
dev.off()

survdiff(formula = surv_object ~ IMPACT_risk, data = V4_TCGA_sample_summary)
pairwise_survdiff(Surv(OS_mo_cap_5y, OS_status_cap_5y) ~ IMPACT_risk, data = V4_TCGA_sample_summary, p.adjust.method = "BH")


cox_model_IMPACT_risk <- coxph(Surv(OS_mo_cap_5y, OS_status_cap_5y) ~ IMPACT_risk + AGE + SEX + Stage, data = V4_TCGA_sample_summary)
summary(cox_model_IMPACT_risk)

pdf("TCGA_cox_model_IMPACT_risk_clinical.pdf",width=10,height=6)
ggforest(
  model = cox_model_IMPACT_risk,
  data = V4_TCGA_sample_summary,
  main = "Hazard ratio",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 1.0,
  refLabel = "Reference",
  noDigits = 2
)
dev.off()



# Now we examine if IMPACT can be leveraged for further patient stratification on top of TIL
# We obtained TIL information from a previous study (Immunity 2018,48(4):812-830.e14), which was based on bulk deconvolution
TCGA_TIL_fraction <- readRDS("../data/V4_TCGA_TIL_fraction.rds")
V4_TCGA_sample_summary <- left_join(V4_TCGA_sample_summary, TCGA_TIL_fraction, by="SAMPLE_ID")
TIL_Regional_Fraction_median <- median(V4_TCGA_sample_summary$TIL.Regional.Fraction, na.rm=TRUE)
V4_TCGA_sample_summary$TIL_Regional_Fraction_dichom <- NA
V4_TCGA_sample_summary$TIL_Regional_Fraction_dichom[which(V4_TCGA_sample_summary$TIL.Regional.Fraction < TIL_Regional_Fraction_median)] <- "Low"
V4_TCGA_sample_summary$TIL_Regional_Fraction_dichom[which(V4_TCGA_sample_summary$TIL.Regional.Fraction >= TIL_Regional_Fraction_median)] <- "High"
# group patients by TIL fraction and IMPACT archetype
V4_TCGA_sample_summary_filter <- filter(V4_TCGA_sample_summary, !is.na(V4_TCGA_sample_summary$TIL_Regional_Fraction_dichom) & !is.na(V4_TCGA_sample_summary$IMPACT_patient))
V4_TCGA_sample_summary_filter$IMPACT_TIL <- paste0(V4_TCGA_sample_summary_filter$IMPACT_patient, "_", V4_TCGA_sample_summary_filter$TIL_Regional_Fraction_dichom)

surv_object2 <- Surv(time = V4_TCGA_sample_summary_filter$OS_mo_cap_5y, event = V4_TCGA_sample_summary_filter$OS_status_cap_5y)
surv_fit_IMPACT_TIL <- survfit(surv_object2 ~ IMPACT_TIL, data = V4_TCGA_sample_summary_filter)

pdf("TCGA_survival_IMPACT_TIL.pdf",width=7,height=5)
ggsurvplot(surv_fit_IMPACT_TIL, 
           data = V4_TCGA_sample_summary_filter,
           legend.labs = levels(V4_TCGA_sample_summary_filter$IMPACT_TIL), 
           palette=c("#e31a1c","#fed976","#01665e","#80cdc1"),
           pval = TRUE, 
           conf.int = FALSE,
           surv.median.line = "hv",
           risk.table = TRUE,
           ggtheme = theme_classic())
dev.off()

survdiff(formula = surv_object2 ~ IMPACT_TIL, data = V4_TCGA_sample_summary_filter)
pairwise_survdiff(Surv(OS_mo_cap_5y, OS_status_cap_5y) ~ IMPACT_TIL, data = V4_TCGA_sample_summary_filter, p.adjust.method = "BH")


V4_TCGA_sample_summary_TIL_high <- filter(V4_TCGA_sample_summary, TIL_Regional_Fraction_dichom == "High")
surv_object3 <- Surv(time = V4_TCGA_sample_summary_TIL_high$OS_mo_cap_5y, event = V4_TCGA_sample_summary_TIL_high$OS_status_cap_5y)
surv_fit_IMPACT_risk_TIL_high <- survfit(surv_object3 ~ IMPACT_risk, data = V4_TCGA_sample_summary_TIL_high)

pdf("TCGA_survival_IMPACT_risk_TIL_high.pdf",width=7,height=5)
ggsurvplot(surv_fit_IMPACT_risk_TIL_high, 
           data = V4_TCGA_sample_summary_TIL_high,
           legend.labs = levels(V4_TCGA_sample_summary_TIL_high$IMPACT_risk), 
           palette=c("#e31a1c","#fed976","#01665e"),
           pval = TRUE, 
           conf.int = FALSE,
           surv.median.line = "hv",
           risk.table = TRUE,
           ggtheme = theme_classic())
dev.off()

survdiff(formula = surv_object3 ~ IMPACT_risk, data = V4_TCGA_sample_summary_TIL_high)
pairwise_survdiff(Surv(OS_mo_cap_5y, OS_status_cap_5y) ~ IMPACT_risk, data = V4_TCGA_sample_summary_TIL_high, p.adjust.method = "BH")


V4_TCGA_sample_summary_TIL_low <- filter(V4_TCGA_sample_summary, TIL_Regional_Fraction_dichom == "Low")
surv_object4 <- Surv(time = V4_TCGA_sample_summary_TIL_low$OS_mo_cap_5y, event = V4_TCGA_sample_summary_TIL_low$OS_status_cap_5y)
surv_fit_IMPACT_risk_TIL_low <- survfit(surv_object4 ~ IMPACT_risk, data = V4_TCGA_sample_summary_TIL_low)

pdf("TCGA_survival_IMPACT_risk_TIL_low.pdf",width=7,height=5)
ggsurvplot(surv_fit_IMPACT_risk_TIL_low, 
           data = V4_TCGA_sample_summary_TIL_low,
           legend.labs = levels(V4_TCGA_sample_summary_TIL_low$IMPACT_risk), 
           palette=c("#e31a1c","#fed976","#01665e"),
           pval = TRUE, 
           conf.int = FALSE,
           surv.median.line = "hv",
           risk.table = TRUE,
           ggtheme = theme_classic())
dev.off()

survdiff(formula = surv_object4 ~ IMPACT_risk, data = V4_TCGA_sample_summary_TIL_low)
pairwise_survdiff(Surv(OS_mo_cap_5y, OS_status_cap_5y) ~ IMPACT_risk, data = V4_TCGA_sample_summary_TIL_low, p.adjust.method = "BH")


