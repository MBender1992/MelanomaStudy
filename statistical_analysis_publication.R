# <<<<<<<<<<<<< HEAD

# load packages
library(tidyverse)
library(ggpubr)
library(rstatix)
library(devtools)
library(data.table)
library(table1)

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
dat_BRAF_tidy$BRAFnew %>% table()
# transform data
dat_BRAF_tidy <- dat %>% 
  filter(miRExpAssess == 1 & !is.na(prior_BRAF_therapy) & !is.na(BRAF) & Stadium == "IV") %>%
  mutate(BRAFnew = ifelse(BRAF == "neg", "wt", ifelse(prior_BRAF_therapy == 1, "mut_BRAF_Inhibitor", "mut_no_BRAF_Inhibitor"))) %>%
  mutate(BRAFnew = factor(BRAFnew, levels = c("wt", "mut_no_BRAF_Inhibitor", "mut_BRAF_Inhibitor"))) %>%
  mutate(BRAF = factor(BRAF, levels = c("neg","pos"), labels = c("wt", "mut"))) %>%
  gather(miRNA, expression, contains("hsa")) %>%
  mutate(logexp = log2(expression)) 

# calculate miRNAs with significant differences 
hccmAnova <- dat_BRAF_tidy %>% 
  group_by(miRNA, Stadium) %>% 
  anova_test(logexp ~BRAFnew, white.adjust = T) %>% 
  filter(p < 0.05) %>% 
  print(n="all")

# calculate GH post-hoc test 
GH_test <- dat_BRAF_tidy %>% 
  group_by(miRNA, Stadium) %>% 
  games_howell_test(logexp ~BRAFnew)%>% 
  filter(p.adj < 0.05) 

# calculate maximum value for each miRNA
maxexp <- dat_BRAF_tidy %>% 
  filter(miRNA %in% GH_test$miRNA) %>%
  group_by(miRNA) %>% 
  summarize(y.position = max(expression)*1.06)

# calculate y.positions
ypos <- left_join(data.frame(miRNA =stat_test$miRNA),  data.frame(maxexp)) %>%
  mutate(y.position = ifelse(duplicated(y.position), y.position*1.11, y.position))

# change y.positions to previously calculated y positions (as stat test was performed on logs so the y.positions are off)
stat_test <- GH_test %>% 
  add_xy_position(x = "BRAFnew") %>%
  mutate(y.position = ypos$y.position)

# plot
plot_BRAF <- dat_BRAF_tidy %>% 
  filter(miRNA %in% hccmAnova$miRNA)  %>%
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

png("Results/BRAF.png", units="in", width=7, height=5, res=1200)
plot_BRAF
dev.off()




