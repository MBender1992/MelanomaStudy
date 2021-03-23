# <<<<<<<<<<<<< HEAD

# load packages
library(tidyverse)
library(ggpubr)
library(rstatix)
library(devtools)
library(data.table)
library(table1)
library(ComplexHeatmap)
library(circlize)
library(RBiomirGS)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

# load data with custom function for melanoma data only for Responders
dat <- load_melanoma_data() # n = 101 patients



#####################################
#                                   #
#         1. patient table          #
#                                   #
#####################################

dat_table1 <- dat
setDT(dat_table1)

# define which factors to display in table
dat_table1$sex <- factor(dat_table1$sex, levels = c("m", "w") , labels = c("Male", "Female"))
dat_table1$miRExpAssess <- factor(dat_table1$miRExpAssess, levels = c(0, 1) , labels = c("no", "yes"))
dat_table1$Responder <- factor(dat_table1$Responder, levels = c("no", "yes",2) , labels = c("no", "yes","P-value"))
dat_table1$adjuvant_IFN <- factor(dat_table1$adjuvant_IFN, levels = c("nein", "ja") , labels = c("no", "yes"))
dat_table1$brainMet <- factor(dat_table1$brainMet, levels = c("nein", "ja") , labels = c("no", "yes"))
dat_table1$subtype <- factor(dat_table1$subtype, levels = c("cutanes Melanom", "Schleimhautmelanom") , labels = c("cutaneous", "mucosal"))
dat_table1$ECOG <- factor(dat_table1$ECOG, levels = c(0,1,2) , labels = c("0", "1", "2"))
dat_table1$Stadium <- factor(dat_table1$Stadium, levels = c("II", "III","IV") , labels = c("II", "III","IV"))
dat_table1$prior_BRAF_therapy <- factor(dat_table1$prior_BRAF_therapy, levels = c(0, 1) , labels = c("no", "yes"))


# define labels for the table
label(dat_table1$Alter)      <- "Age (years)"
label(dat_table1$BRAF)      <- "BRAF-status"
label(dat_table1$Stadium)  <- "AJCC stage" # add Stadium to source table
label(dat_table1$therapy_at_blood_draw) <- "Therapy at blood draw"
label(dat_table1$sex)  <- "Sex"
label(dat_table1$Responder)  <- "Immunotherapy response"
label(dat_table1$ECOG)      <- "ECOG"
label(dat_table1$breslow_thickness_mm)      <- "Breslow thickness (mm)" # change to double
label(dat_table1$subtype) <- "Subtype"
label(dat_table1$localization) <- "Localization"
label(dat_table1$brainMet) <- "Brain metastasis"
label(dat_table1$miRExpAssess) <- "miRNA expression measured"
label(dat_table1$adjuvant_IFN) <- "Received adjuvant IFN treatment"
label(dat_table1$prior_BRAF_therapy) <- "Received prior anti-BRAF therapy"


# define text for footnote
fn <- "Statistical test: Unequal variance t-test (welch's t-test) for numerical data and chi² test for categorical data. Raw p-values are shown."

table1(~ Alter + BRAF +  prior_BRAF_therapy + Stadium + miRExpAssess + adjuvant_IFN + brainMet + sex + ECOG + breslow_thickness_mm + subtype + localization | Responder,
       data=dat_table1, droplevels=F , render=rndr, render.strat=rndr.strat, footnote = fn)



#####################################
#                                   #
#         2. Serum markers          #
#                                   #
#####################################

# change data structure for easier statistical comparison
dat_serum_markers_tidy <- dat %>% 
  filter(!is.na(Responder)) %>%
  select(c(ID, Responder,Baseline, Eosinophile, CRP, LDH, S100)) %>% 
  gather(serum_marker, value,-c(ID, Responder,Baseline)) %>%
  mutate(log_val = ifelse(is.infinite(log2(value)), 0, log2(value))) %>% 
  filter(!is.na(log_val)) 

dat_plot <- dat_serum_markers_tidy %>% mutate(
  serum_marker = str_replace_all(serum_marker, "CRP", "CRP (mg/L)"),
  serum_marker = str_replace_all(serum_marker, "LDH", "LDH (U/L)"),
  serum_marker = str_replace_all(serum_marker, "S100", "S100 (µg/L)"),
  serum_marker = str_replace_all(serum_marker, "Eosinophile", "Eosinophile (%)")
)

# plot 4 markers in separate plots and calculate statistics
plot_serum_markers <- signif_plot_Melanoma(dat_plot, x="Responder", y="log_val", p.adj = "fdr", 
                     plot.type = "boxplot", significance=FALSE, Legend = FALSE, ylab = "log2 serum marker concentration",
                     method ="t.test", p.label="{p.signif}", facet="serum_marker")

png("serum_markers.png", units="in", width=5, height=4, res=1200)
plot_serum_markers$graph
dev.off()





#####################################
#                                   #
#           3. miRNAs               #
#                                   #
#####################################

# tidy miRNA data.....................................................................................................
dat_miRNA_tidy <- dat %>% 
  # only use data where miRNA data was measured and responder status is known
  filter(miRExpAssess == 1 & !is.na(Responder)) %>%
  gather(miRNA, expression, contains("hsa")) %>%
  mutate(miRNA = str_replace_all(.$miRNA, "hsa-","")) %>%
  mutate(log_exp = log2(expression))

# Plot miRNA data
plot_miRNA <- signif_plot_Melanoma(dat_miRNA_tidy, x="Responder", y="log_exp", signif=0.05, p.adj = "fdr", 
                     plot.type = "boxplot", significance=F, Legend = F, var.equal = F,
                     method ="t.test", p.label="p = {round(p,4)}",p.size = 3, facet="miRNA")

png("miRNAs.png", units="in", width=5.5, height=4, res=1200)
plot_miRNA$graph
dev.off()

#
# res.signif <- str_replace_all(c(plot_miRNA$stat_test_results$miRNA,plot_serum_markers$stat_test_results$serum_marker),"-",".")
# saveRDS(res.signif, "significant_features.rds")





#####################################
#                                   #
#           4. BRAF                 #
#                                   #
#####################################

# 
dat_BRAF <- dat %>% 
  filter(miRExpAssess == 1 & !is.na(prior_BRAF_therapy) & !is.na(BRAF) & Stadium == "IV")


# transform data
dat_BRAF_tidy <- dat_BRAF %>%
  mutate(BRAFnew = ifelse(BRAF == "neg", "wt", ifelse(prior_BRAF_therapy == 1, "mut_BRAF_Inhibitor", "mut_no_BRAF_Inhibitor"))) %>%
  mutate(BRAFnew = factor(BRAFnew, levels = c("wt", "mut_no_BRAF_Inhibitor", "mut_BRAF_Inhibitor"))) %>%
  mutate(BRAF = factor(BRAF, levels = c("neg","pos"), labels = c("wt", "mut"))) %>%
  gather(miRNA, expression, contains("hsa")) %>%
  mutate(logexp = log2(expression)) 

# calculate statistics and add ypositions
stat_test <- stat_test_BRAF(dat_BRAF_tidy, var = "BRAFnew", p.adj.anova = "holm")

# plot
plot_BRAF <- dat_BRAF_tidy %>% 
  filter(miRNA %in% unique(stat_test$miRNA))  %>%
  ggplot(aes(BRAFnew, expression)) + 
  stat_boxplot(geom='errorbar', linetype=1, width=0.4)+
  geom_boxplot(outlier.shape = NA, aes(fill = BRAF))+
  geom_jitter(size = 1.2, shape = 1, position = position_jitter(0.1)) + 
  facet_wrap(~miRNA, scales = "free")+
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        legend.key.size = unit(1,"line"),
        panel.grid.minor=element_blank(),
        strip.background=element_blank()) +
  scale_fill_manual(values=c("#3288bd","#d53e4f")) +
  scale_y_continuous(expand = expansion(mult = c(.05, .15))) +
  ylab("miRNA Expression (a.u.)") +
  stat_pvalue_manual(stat_test) 

png("Results/BRAF/BRAF.png", units="in", width=4, height=2, res=1200)
plot_BRAF
dev.off()


# scale data for Heatmap
dat_scaled <- dat_BRAF %>% 
  column_to_rownames("ID") %>%
  select(contains("hsa")) %>%
  scale() %>%
  t()
  

# set ID as rownames so colorbar works properly
dat_colorbar <- dat_BRAF %>% 
  select(c(ID, BRAF, prior_BRAF_therapy, brainMet)) %>% 
  column_to_rownames("ID") %>%
  mutate(BRAF = (ifelse(BRAF == "neg", "wt", "mut")),
         prior_BRAF_therapy = (ifelse(prior_BRAF_therapy == 1, "yes", "no")),
         brainMet = str_replace_all(brainMet, "ja", "yes"),
         brainMet = str_replace_all(brainMet, "nein", "no"))
         
# define colors for colorbar
colorbar <- HeatmapAnnotation(
  df =dat_colorbar,
  col = list(
    BRAF = c("wt" = "#3288bd","mut"="#d53e4f"),
    prior_BRAF_therapy = c("yes" = "#542788","no" = "#b35806"),
    brainMet = c("yes" = "#d73027","no" = "#1a9850")
  ),
  annotation_legend_param = list(
    BRAF = list(nrow=1),
    prior_BRAF_therapy  = list(nrow=1),
    brainMet = list(nrow=1)
  )
)

col_fun = colorRamp2(c(-1.8, 0,0.3, 1.8), c("#4575b4", "white","#ffffbf", "#d73027"))


Ht <- Heatmap(
  dat_scaled,
  col= col_fun,
  top_annotation = colorbar,
  column_title = c("A", "B"),
  border = T,
  row_split = 3,
  column_km = 2,
  column_km_repeats = 100,
  clustering_method_row = "average",
  clustering_method_columns = "average",
  clustering_distance_row = "pearson",
  clustering_distance_column = "euclidean",
  rect_gp = gpar(col = "white",lty = 1, lwd = 1),
  row_names_gp = gpar(fontsize = 10),
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(
    title = "row Z-score",
    at = seq(-2,2,by=1),
    color_bar="continuous",
    title_position ="topcenter",
    legend_direction = "horizontal",
    legend_width = unit(4, "cm")
  ))


png("Results/Heatmap_BRAF.png", units="in", width=9, height=8, res=600)
draw(Ht, annotation_legend_side = "bottom", heatmap_legend_side = "bottom")
dev.off()



# extract clusters out of Heatmap object
Ht_clusters <- extract_clusters(dat_scaled, Ht, sampleName = "ID", sampleClust = "sampleCluster", geneName = "miRNA", geneClust = "miRCluster")

# transform miRNA data
dat_tmp <- dat_BRAF %>% mutate(ID = as.character(ID)) %>% 
  gather(miRNA, expression, contains("hsa")) %>%
  select(ID, miRNA, expression) %>%
  mutate(log_exp = log2(expression))

# combine miRNA data with cluster information
dat_clusters <- left_join(dat_tmp, Ht_clusters$sampleCluster) %>%
  left_join(Ht_clusters$miRCluster)

# create list for each miRCluster
ls_miRCluster <- split(dat_clusters, f = dat_clusters$miRCluster)         

# extract cluster 1 B (upregulated miRNAs in BRAF wildtype)
cl1B <- as.data.frame(lapply(summary_clusters(ls_miRCluster, 1, "B"), FUN=drop_attr))


# collect mRNA targets of the miRNAs upregulated in cluster 1B
rbiomirgs_mrnascan(objTitle = "cl_1B_predicted", mir = cl1B$miRNA, sp = "hsa", 
                   queryType = "predicted",parallelComputing = TRUE,clusterType = "PSOCK")


# calculate GS by logistic regression
rbiomirgs_logistic(objTitle = "cl_1B_predicted_mirna_mrna_iwls_KEGG",mirna_DE = cl1B, 
                   var_mirnaName = "miRNA",var_mirnaFC = "FC",var_mirnaP = "pvalue", mrnalist = cl_1B_predicted_mrna_entrez_list, 
                   mrna_Weight = NULL, gs_file = "Data/Pathway Analysis/c2.cp.kegg.v7.2.entrez.gmt", optim_method = "IWLS", 
                   p.adj = "fdr", parallelComputing = FALSE, clusterType = "PSOCK")



# total number of significantly enriched pathways
sum(cl_1B_predicted_mirna_mrna_iwls_KEGG_GS$adj.p.val < 0.05)

# remove randomly enriched pathways
ctrl_list_KEGG <- readRDS(file = "Data/Pathway Analysis/ctrl_list_KEGG.rds")
n <- 50
res_ctrl_KEGG <- pathway_ctrl_summary(ctrl_list_KEGG, n=n)
bias_KEGG <- names(res_ctrl_KEGG$bias[res_ctrl_KEGG$bias > 0.1])
cl_1B_KEGG_plot <- cl_1B_predicted_mirna_mrna_iwls_KEGG_GS %>% filter(!GS %in% bias_KEGG)


#plot results (volcano plot)
png("Results/BRAF/Pathway/volcano_cl_1B_KEGG.png", units="in", width=5, height=4, res=600)
rbiomirgs_volcano(gsadfm = cl_1B_KEGG_plot,topgsLabel = TRUE,n = 15,gsLabelSize = 2,
                  sigColour = "red",plotWidth = 250,plotHeight = 220,xLabel = "model coefficient")
dev.off()

# plot distribution of enriched gene sets
png("Results/BRAF/Pathway/volcano_bar_dist_cl_1B_KEGG.png", units="in", width=6, height=3, res=600)
rbiomirgs_bar(gsadfm = cl_1B_KEGG_plot,signif_only = F,gs.name = F,
              n = "all",xLabel = "gene set", yLabel = "model coefficient", plotWidth = 250, plotHeight = 220)
dev.off()

# plot top enriched gene sets
png("Results/BRAF/Pathway/volcano_bar_top15cl_1B_KEGG.png", units="in", width=5, height=4, res=600)
rbiomirgs_bar(gsadfm = cl_1B_KEGG_plot,signif_only = 15,gs.name = T,xLabel = "model coefficient",
              yTxtSize = 7, n = 15, plotWidth = 250, plotHeight = 220)
dev.off()
