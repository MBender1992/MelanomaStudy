# <<<<<<<<<<<<< HEAD

# load packages
library(tidyverse)
library(ggpubr)
library(rstatix)
library(devtools)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

#load data
url_file <- "https://raw.githubusercontent.com/MBender1992/PhD/Marc/Data/200619_chronic_irr_normalized.csv" 
dat <-  load_Fireplex_data_PhD(filename = url(url_file), threshold = 2.5)

load_melanoma_data(characterAsFactor = TRUE) 
# set working directory
setwd("Z:/Aktuell/Eigene Dateien/Eigene Dateien_Marc/R/Projekte/Doktorarbeiten_Melanom_Liquid_Biopsies/Daten")

# load data with custom function for melanoma data only for Responders
dat_combined <- load_melanoma_data(characterAsFactor = TRUE) %>% filter(!is.na(Responder))

test <- dat_combined

# tidy miRNA data.....................................................................................................
dat_miRNA_tidy <- dat_combined %>% 
  gather(miRNA, expression, contains("mir")) %>%
  mutate(miRNA = str_replace_all(.$miRNA, "hsa-","")) %>%
  mutate(log_exp = log2(expression)) 

dat_combined %>% select(contains("mir"))

# tidy lab parameter data.............................................................................................
dat_lab_pars_tidy <- dat_combined %>%
  filter(!is.na(CRP) & !is.na(LDH)  &!is.na(S100)) %>%
  select(c(ID, Responder,Baseline, Eosinophile, CRP, LDH, S100)) %>% 
  gather(lab_parameter, value,-c(ID, Responder,Baseline)) %>%
  mutate(log_val = ifelse(is.infinite(log2(value)), 0, log2(value))) 





# Plot miRNA data
plot_miRNA <- signif_plot_Melanoma(dat_miRNA_tidy, x="Responder", y="log_exp", signif=0.05,
                     plot.type = "dotplot", significance=T, Legend = F, var.equal = F,
                     method ="t.test", p.label="p = {round(p,4)}",p.size = 3, facet="miRNA")

png("miRNAs.png", units="in", width=7, height=7, res=1200)
plot_miRNA$graph
dev.off()



# Plot lab parameter data
plot_lab_pars <- signif_plot_Melanoma(dat_lab_pars_tidy, x="Responder", y="log_val", 
                                   plot.type = "dotplot", significance=F, Legend = F, ylab = "log2 lab parameter concentration",
                                   method ="t.test", p.label="{p.signif}", facet="lab_parameter")

png("lab_pars.png", units="in", width=5, height=4, res=1200)
plot_lab_pars$graph
dev.off()





