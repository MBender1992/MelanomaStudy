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

# Tabelle 1
dat_table1 <- dat
setDT(dat_table1)

dat_table1$sex <- factor(dat_table1$sex, levels = c("m", "w") , labels = c("Male", "Female"))
dat_table1$miRExpAssess <- factor(dat_table1$miRExpAssess, levels = c(0, 1) , labels = c("no", "yes"))


# dat$trt     <- factor(dat$trt, levels=1:2, labels=c("D-penicillamine", "Placebo"))
# dat$sex     <- factor(dat$sex, levels=c("m", "f"), labels=c("Male", "Female"))
# dat$stage   <- factor(dat$stage, levels=1:4, labels=paste("Stage", 1:4))
# dat$edema   <- factor(dat$edema, levels=c(0, 0.5, 1),
#                       labels=c("No edema",
#                                "Untreated or successfully treated",
#                                "Edema despite diuretic therapy"))
# dat$spiders <- as.logical(dat$spiders)
# dat$hepato  <- as.logical(dat$hepato)
# dat$ascites <- as.logical(dat$ascites)


# define labels for the table
label(dat_table1$Alter)      <- "age (years)"
label(dat_table1$BRAF)      <- "BRAF-status"
label(dat_table1$Stadium)  <- "AJCC stage" # add Stadium to source table
label(dat_table1$therapy_at_blood_draw) <- "therapy at blood draw"
label(dat_table1$sex)  <- "sex"
label(dat_table1$Responder)  <- "immunotherapy response"
label(dat_table1$ECOG)      <- "ECOG"
label(dat_table1$breslow_thickness_mm)      <- "breslow thickness (mm)" # change to double
label(dat_table1$subtype) <- "subtype"
label(dat_table1$localization) <- "localization"
label(dat_table1$Hirnmetastase) <- "brain metastasis"
label(dat_table1$miRExpAssess) <- "miRNA expression measured"
label(dat_table1$adjuvant_IFN) <- "received adjuvant IFN treatment"


table1(~ Alter + BRAF + Stadium + miRExpAssess + adjuvant_IFN + Hirnmetastase + sex + ECOG + breslow_thickness_mm + subtype + localization | Responder, data=dat_table1)

 



# tidy miRNA data.....................................................................................................
dat_miRNA_tidy <- dat %>% 
  gather(miRNA, expression, contains("mir")) %>%
  mutate(miRNA = str_replace_all(.$miRNA, "hsa-","")) %>%
  mutate(log_exp = log2(expression)) 
















# Plot miRNA data
plot_miRNA <- signif_plot_Melanoma(dat_miRNA_tidy, x="Responder", y="log_exp", signif=0.05,
                     plot.type = "dotplot", significance=T, Legend = F, var.equal = F,
                     method ="t.test", p.label="p = {round(p,4)}",p.size = 3, facet="miRNA")

png("miRNAs.png", units="in", width=7, height=7, res=1200)
plot_miRNA$graph
dev.off()





# tidy lab parameter data.............................................................................................
dat_lab_pars_tidy <- dat %>%
  filter(!is.na(CRP) & !is.na(LDH)  &!is.na(S100)) %>%
  select(c(ID, Responder,Baseline, Eosinophile, CRP, LDH, S100)) %>% 
  gather(lab_parameter, value,-c(ID, Responder,Baseline)) %>%
  mutate(log_val = ifelse(is.infinite(log2(value)), 0, log2(value))) 


# Plot lab parameter data
plot_lab_pars <- signif_plot_Melanoma(dat_lab_pars_tidy, x="Responder", y="log_val", 
                                   plot.type = "dotplot", significance=F, Legend = F, ylab = "log2 lab parameter concentration",
                                   method ="t.test", p.label="{p.signif}", facet="lab_parameter")

png("lab_pars.png", units="in", width=5, height=4, res=1200)
plot_lab_pars$graph
dev.off()





