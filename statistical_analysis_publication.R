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


# table(dat$Responder, dat$prior_BRAF_therapy)
# chisq.test(dat$Responder, dat$prior_BRAF_therapy)

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
dat_table1$Responder <- factor(dat_table1$Responder, levels = c("nein", "ja",2) , labels = c("no", "yes","P-value"))
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

table1(~ Alter + BRAF + prior_BRAF_therapy + Stadium + miRExpAssess + adjuvant_IFN + brainMet + sex + ECOG + breslow_thickness_mm + subtype + localization | Responder,
       data=dat_table1, droplevels=F, render=rndr, render.strat=rndr.strat, footnote = fn)



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
  mutate(log_val = ifelse(is.infinite(log2(value)), 0, log2(value)),
         Responder =  factor(Responder, levels = c("nein", "ja") , labels = c("no", "yes"))) %>% 
  filter(!is.na(log_val)) 

# plot 4 markers in separate plots and calculate statistics
plot_serum_markers <- signif_plot_Melanoma(dat_serum_markers_tidy, x="Responder", y="log_val", p.adj = "fdr", 
                     plot.type = "dotplot", significance=FALSE, Legend = FALSE, ylab = "log2 serum marker concentration",
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
  mutate(log_exp = log2(expression),
         Responder =  factor(Responder, levels = c("nein", "ja") , labels = c("no", "yes")))

# Plot miRNA data
plot_miRNA <- signif_plot_Melanoma(dat_miRNA_tidy, x="Responder", y="log_exp", signif=0.05, p.adj = "fdr", 
                     plot.type = "dotplot", significance=F, Legend = F, var.equal = F,
                     method ="t.test", p.label="p = {round(p,4)}",p.size = 3, facet="miRNA")

png("miRNAs.png", units="in", width=5.5, height=4, res=1200)
plot_miRNA$graph
dev.off()

#
# res.signif <- str_replace_all(c(plot_miRNA$stat_test_results$miRNA,plot_serum_markers$stat_test_results$serum_marker),"-",".")
# saveRDS(res.signif, "significant_features.rds")
