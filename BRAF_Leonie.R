# <<<<<<<<<<<<< HEAD

# load packages
library(tidyverse)
library(ggpubr)
library(rstatix)
library(devtools)
library(ComplexHeatmap)
library(circlize)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  


#####################################
#                                   #
#         1. Load Data              #
#                                   #
#####################################

dat_miR   <- read_csv("Data/miRNA_Expression_Fireplex_Melanoma_Study.csv")

# change ID column to uniform capital letters for later filtering
names(dat_miR) <- c("miRNA", toupper(names(dat_miR)[-1]))

# define IDs to be dropped for further analyses
duplicates <- c("22B","38B","39B","47B")

# wide miR data (78 patients with miRNA data)
dat_miR_trans <- transpose_dataframe(colnames = c("ID",dat_miR$miRNA), data = dat_miR) %>%
  filter(!ID %in% duplicates) 

# additional data
dat_meta  <- read_xlsx("Data/Metadata_Melanoma_Study.xlsx") %>%
  select(-c(therapy_start, Abnahmedatum)) %>%
  mutate(TRIM_PDL1_Expression = str_replace_all(TRIM_PDL1_Expression,"\\++","+")) %>% 
  mutate(TRIM_PDL1_Expression = ifelse(TRIM_PDL1_Expression == "o", NA,TRIM_PDL1_Expression)) %>%
  mutate(Stadium = toupper(Stadium)) %>%
  mutate(Stadium = str_extract(Stadium, "^[IV]{1,3}")) %>%
  mutate(BRAF = str_replace_all(BRAF, "\\.", "")) %>% 
  mutate(breslow_thickness_mm = parse_number(breslow_thickness_mm))

# combine data with control data added
dat_ctrl <- left_join(dat_miR_trans, dat_meta %>% mutate(ID = as.character(ID)), by = "ID") %>%
  mutate(group = ifelse(str_detect(ID, "K"),"control", "patient")) %>% 
  mutate(BRAF = ifelse(BRAF == "pos", "mut", "wt")) %>% 
  mutate(subgroup = ifelse(BRAF == "mut" & brainMet == "nein", "NBM_mut", 
                           ifelse(BRAF == "mut" & brainMet == "ja", "BM_mut",
                                  ifelse(BRAF == "wt" & brainMet == "nein", "NBM_wt",
                                         ifelse(BRAF == "wt" & brainMet == "ja", "BM_wt",NA))))) %>%
  mutate(subgroup = ifelse(group == "control", "control", subgroup)) %>% 
  mutate(prior_BRAF_therapy = ifelse(str_detect(Vorbehandlung,"Mek|Dabra|Tafinlar|Tefinlar|MEK|BRAF|Vemu|[zZ]ellboraf"), "yes", "no")) %>%
  mutate(subgroup = ifelse(subgroup == "BM_mut" & prior_BRAF_therapy == "yes", "BM_mut_therapy_yes", subgroup),
         subgroup = ifelse(subgroup == "BM_mut" & prior_BRAF_therapy == "no", "BM_mut_therapy_no", subgroup)) %>%
  mutate(subgroup = ifelse(subgroup == "NBM_mut" & prior_BRAF_therapy == "yes", "NBM_mut_therapy_yes", subgroup),
         subgroup = ifelse(subgroup == "NBM_mut" & prior_BRAF_therapy == "no", "NBM_mut_therapy_no", subgroup)) %>%
  filter(!is.na(subgroup))




#####################################
#                                   #
#  2. Compare NBM_mut and control   #
#                                   #
#####################################


# Compare NBM_mut and control
NBM_mut_vs_ctrl <- dat_ctrl %>% filter(subgroup %in% c("NBM_mut_therapy_yes","NBM_mut_therapy_no", "control"))
NBM_mut_vs_ctrl_tidy <- NBM_mut_vs_ctrl %>% 
  gather(miRNA, expression, contains("hsa")) %>%
  mutate(logexp = log2(expression))

# calculate statistics and add ypositions
stat_test <- stat_test_BRAF(NBM_mut_vs_ctrl_tidy, var = "subgroup", p.adj.anova = "holm")

# scale data for Heatmap
dat_scaled <- NBM_mut_vs_ctrl %>%
  select(ID,unique(stat_test$miRNA)) %>%
  column_to_rownames("ID") %>%
  scale() %>%
  t()

# set ID as rownames so colorbar works properly
dat_colorbar <- NBM_mut_vs_ctrl %>% 
  select(c(ID, subgroup)) %>% 
  column_to_rownames("ID")

# define colors for colorbar
colorbar <- HeatmapAnnotation(
  df =dat_colorbar,
  col = list(subgroup = c("NBM_mut_therapy_yes" = "#3288bd","NBM_mut_therapy_no"="#d53e4f", "control"= "#1a9850")),
  annotation_legend_param = list(subgroup = list(nrow=1), labels = c("healthy control", "NBM_mut: no BRAF therapy","NBM_mut: BRAF therapy"))
)


Ht <- Heatmap(
  dat_scaled,
  col= col_fun,
  top_annotation = colorbar,
  column_title = c("1", "2"),
  border = T,
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

png("Results/BRAF/Heatmap_NBM_mut_vs_ctrl.png", units="in", width=10, height=3, res=600)
draw(Ht, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()

# plot significant miRNAs
p1 <- NBM_mut_vs_ctrl_tidy  %>%
  filter(miRNA %in% unique(stat_test$miRNA)) %>%
  ggplot(aes(subgroup, expression)) +
  facet_wrap(~miRNA, scales = "free", nrow = 2) + 
  stat_boxplot(geom='errorbar', linetype=1, width=0.4)+
  geom_boxplot(outlier.shape = NA,aes(fill = group)) +
  geom_jitter(size = 1.2, shape = 1, position = position_jitter(0.1)) + 
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        legend.key.size = unit(1,"line"),
        panel.grid.minor=element_blank(),
        strip.background=element_blank())+
  scale_fill_manual(values=c("#3288bd","#d53e4f")) +
  scale_y_continuous(expand = expansion(mult = c(.05, .15))) +
  stat_pvalue_manual(stat_test) 


png("Results/BRAF/Boxplot_NBM_mut_vs_ctrl.png", units="in", width=8, height=3, res=1200)
p1
dev.off()



#####################################
#                                   #
#  3. Compare BM_wt and control     #
#                                   #
#####################################



# Compare NBM_mut and control
BM_wt_vs_ctrl <- dat_ctrl %>% filter(subgroup %in% c("BM_wt", "control"))
BM_wt_vs_ctrl_tidy <- BM_wt_vs_ctrl %>% 
  gather(miRNA, expression, contains("hsa")) %>%
  mutate(logexp = log2(expression))

# calculate statistics and add ypositions
stat_test <- stat_test_BRAF(BM_wt_vs_ctrl_tidy, var = "subgroup", p.adj.anova = "holm")

# plot significant miRNAs
p2 <- BM_wt_vs_ctrl_tidy  %>%
  filter(miRNA %in% unique(stat_test$miRNA)) %>%
  ggplot(aes(subgroup, expression)) +
  stat_boxplot(geom='errorbar', linetype=1, width=0.4)+
  geom_boxplot(outlier.shape = NA,aes(fill = group)) +
  geom_jitter(size = 1.2, shape = 1, position = position_jitter(0.1)) + 
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        legend.key.size = unit(1,"line"),
        panel.grid.minor=element_blank(),
        strip.background=element_blank())+
  scale_fill_manual(values=c("#3288bd","#d53e4f")) +
  scale_y_continuous(expand = expansion(mult = c(.05, .15))) +
  stat_pvalue_manual(stat_test) 

png("Results/BRAF/Boxplot_BM_wt_vs_ctrl.png", units="in", width=3, height=2, res=1200)
p2
dev.off()




#####################################
#                                   #
#  4. Compare BM_mut and control    #
#                                   #
#####################################


# Compare NBM_mut and control
BM_mut_vs_ctrl <- dat_ctrl %>% filter(subgroup %in% c("BM_mut_therapy_yes", "BM_mut_therapy_no", "control"))
BM_mut_vs_ctrl_tidy <- BM_mut_vs_ctrl %>% 
  gather(miRNA, expression, contains("hsa")) %>%
  mutate(logexp = log2(expression))

# calculate statistics and add ypositions
stat_test <- stat_test_BRAF(BM_mut_vs_ctrl_tidy, var = "subgroup", p.adj.anova = "holm")

# scale data for Heatmap
dat_scaled <- BM_mut_vs_ctrl %>%
  select(ID,unique(stat_test$miRNA)) %>%
  column_to_rownames("ID") %>%
  scale() %>%
  t()

# set ID as rownames so colorbar works properly
dat_colorbar <- BM_mut_vs_ctrl %>% 
  select(c(ID, subgroup)) %>% 
  column_to_rownames("ID")

# define colors for colorbar
colorbar <- HeatmapAnnotation(
  df =dat_colorbar,
  col = list(subgroup = c("BM_mut_therapy_yes" = "#3288bd","BM_mut_therapy_no"="#d53e4f", "control"= "#1a9850")),
  annotation_legend_param = list(subgroup = list(nrow=1), labels = c("BM_mut: no BRAF therapy", "BM_mut: BRAF therapy","healthy control"))
)


Ht <- Heatmap(
  dat_scaled,
  col= col_fun,
  top_annotation = colorbar,
  column_title = c("1", "2"),
  border = T,
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

png("Results/BRAF/Heatmap_BM_mut_vs_ctrl.png", units="in", width=10, height=2, res=600)
draw(Ht, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()


# plot significant miRNAs
p3 <- BM_mut_vs_ctrl_tidy  %>%
  filter(miRNA %in% unique(stat_test$miRNA)) %>%
  ggplot(aes(subgroup, expression)) +
  facet_wrap(~miRNA, scales = "free") +
  stat_boxplot(geom='errorbar', linetype=1, width=0.4)+
  geom_boxplot(outlier.shape = NA,aes(fill = group)) +
  geom_jitter(size = 1.2, shape = 1, position = position_jitter(0.1)) + 
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        legend.key.size = unit(1,"line"),
        panel.grid.minor=element_blank(),
        strip.background=element_blank())+
  scale_fill_manual(values=c("#3288bd","#d53e4f")) +
  scale_y_continuous(expand = expansion(mult = c(.05, .15))) +
  stat_pvalue_manual(stat_test) 

png("Results/BRAF/Boxplot_BM_mut_vs_ctrl.png", units="in", width=4, height=2, res=1200)
p3
dev.off()




#####################################
#                                   #
#  5. Compare BM_mut and NBM mut    #
#                                   #
#####################################


# Compare NBM_mut and control
NBM_wt_vs_NBM_mut <- dat_ctrl %>% filter(subgroup %in% c("NBM_wt", "NBM_mut_therapy_yes", "NBM_mut_therapy_no"))
NBM_wt_vs_NBM_mut_tidy <- NBM_wt_vs_NBM_mut %>% 
  gather(miRNA, expression, contains("hsa")) %>%
  mutate(logexp = log2(expression))

# calculate statistics and add ypositions
# stat_test <- stat_test_BRAF(NBM_wt_vs_NBM_mut_tidy, var = "subgroup", p.adj.anova = "holm")
#
test <- c("hsa-mir-23a-3p","hsa-mir-23b-3p","hsa-mir-145-5p","hsa-mir-150-5p","hsa-mir-185-5p","hsa-mir-205-5p","hsa-mir-214-3p","hsa-mir-342-3p")

# plot significant miRNAs
p4 <- NBM_wt_vs_NBM_mut_tidy  %>%
  filter(miRNA %in% test) %>%
  ggplot(aes(BRAF, expression)) +
  facet_wrap(~miRNA, scales = "free") +
  stat_boxplot(geom='errorbar', linetype=1, width=0.4)+
  geom_boxplot(outlier.shape = NA,aes(fill = BRAF)) +
  geom_jitter(size = 1.2, shape = 1, position = position_jitter(0.1)) + 
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        legend.key.size = unit(1,"line"),
        panel.grid.minor=element_blank(),
        strip.background=element_blank())+
  scale_fill_manual(values=c("#3288bd","#d53e4f")) +
  scale_y_continuous(expand = expansion(mult = c(.05, .15))) +
  stat_pvalue_manual(stat_test) 


# scale data for Heatmap
dat_scaled <- NBM_wt_vs_NBM_mut %>%
  select(ID,test) %>%
  column_to_rownames("ID") %>%
  scale() %>%
  t()

# set ID as rownames so colorbar works properly
dat_colorbar <- NBM_wt_vs_NBM_mut %>% 
  select(c(ID, subgroup)) %>% 
  column_to_rownames("ID")

# define colors for colorbar
colorbar <- HeatmapAnnotation(
  df =dat_colorbar,
  col = list(subgroup = c("NBM_mut_therapy_yes" = "#3288bd","NBM_mut_therapy_no"="#d53e4f", "NBM_wt"= "#1a9850")),
  annotation_legend_param = list(subgroup = list(nrow=1), labels = c("NBM_mut: no BRAF therapy", "NBM_mut: BRAF therapy","NBM_wt"))
)


Ht <- Heatmap(
  dat_scaled,
  col= col_fun,
  top_annotation = colorbar,
  column_title = c("1", "2"),
  border = T,
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

png("Results/BRAF/Heatmap_NBM_mut_vs_NBM_wt.png", units="in", width=10, height=3, res=600)
draw(Ht, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()






#####################################
#                                   #
#  6. Compare BM_wt and NBM_wt      #
#                                   #
#####################################

# Compare NBM_mut and control
BM_wt_vs_NBM_wt <- dat_ctrl %>% filter(subgroup %in% c("NBM_wt", "BM_wt"))
BM_wt_vs_NBM_wt_tidy <- BM_wt_vs_NBM_wt %>% 
  gather(miRNA, expression, contains("hsa")) %>%
  mutate(logexp = log2(expression))

# calculate statistics and add ypositions
stat_test <- stat_test_BRAF(BM_wt_vs_NBM_wt_tidy, var = "subgroup", p.adj.anova = "holm")

p5 <- BM_wt_vs_NBM_wt_tidy %>%
  filter(miRNA %in% unique(stat_test$miRNA)) %>%
  ggplot(aes(subgroup, expression))  +
  stat_boxplot(geom='errorbar', linetype=1, width=0.4)+
  geom_boxplot(outlier.shape = NA,aes(fill = brainMet)) +
  geom_jitter(size = 1.2, shape = 1, position = position_jitter(0.1)) + 
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        legend.key.size = unit(1,"line"),
        panel.grid.minor=element_blank(),
        strip.background=element_blank())+
  scale_fill_manual(values=c("#3288bd","#d53e4f")) +
  scale_y_continuous(expand = expansion(mult = c(.05, .15))) +
  stat_pvalue_manual(stat_test) 

png("Results/BRAF/Boxplot_BM_wt_vs_NBM_wt.png", units="in", width=3, height=2, res=1200)
p5
dev.off()



# #####################################
# #                                   #
# #  7. BM vs NBM                     #
# #                                   #
# #####################################
# 
# # Compare NBM_mut and control
# BM_vs_NBM_tidy <-  dat_ctrl %>% 
#   mutate(brainMet = ifelse(group == "control", "control", brainMet)) %>%
#   gather(miRNA, expression, contains("hsa")) %>%
#   mutate(logexp = log2(expression))
# 
# # calculate statistics and add ypositions
# stat_test <- stat_test_BRAF(BM_vs_NBM_tidy, var = "brainMet", p.adj.anova = "holm")
# 
# 
# p5 <- BM_vs_NBM_tidy  %>% filter(group != "control")
#   filter(miRNA %in% unique(stat_test$miRNA)) %>%
#   ggplot(aes(brainMet, expression))  +
#   facet_wrap(~miRNA, scales = "free") + 
#   stat_boxplot(geom='errorbar', linetype=1, width=0.4)+
#   geom_boxplot(outlier.shape = NA,aes(fill = brainMet)) +
#   geom_jitter(size = 1.2, shape = 1, position = position_jitter(0.1)) + 
#   theme_bw() +
#   theme(axis.title.x=element_blank(),
#         axis.text.x = element_blank(),
#         legend.key.size = unit(1,"line"),
#         panel.grid.minor=element_blank(),
#         strip.background=element_blank())+
#   # scale_fill_manual(values=c("#3288bd","#d53e4f")) +
#   scale_y_continuous(expand = expansion(mult = c(.05, .15))) +
#   stat_pvalue_manual(stat_test) 
