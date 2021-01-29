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

# define which factors to display in table
dat_table1$sex <- factor(dat_table1$sex, levels = c("m", "w") , labels = c("Male", "Female"))
dat_table1$miRExpAssess <- factor(dat_table1$miRExpAssess, levels = c(0, 1) , labels = c("no", "yes"))
dat_table1$Responder <- factor(dat_table1$Responder, levels = c("nein", "ja",2) , labels = c("no", "yes","P-value"))
dat_table1$adjuvant_IFN <- factor(dat_table1$adjuvant_IFN, levels = c("nein", "ja") , labels = c("no", "yes"))
dat_table1$Hirnmetastase <- factor(dat_table1$Hirnmetastase, levels = c("nein", "ja") , labels = c("no", "yes"))
dat_table1$subtype <- factor(dat_table1$subtype, levels = c("cutanes Melanom", "Schleimhautmelanom") , labels = c("cutaneous", "mucosal"))
dat_table1$ECOG <- factor(dat_table1$ECOG, levels = c(0,1,2) , labels = c("0", "1", "2"))
dat_table1$Stadium <- factor(dat_table1$Stadium, levels = c("II", "III","IV") , labels = c("II", "III","IV"))

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
label(dat_table1$Hirnmetastase) <- "Brain metastasis"
label(dat_table1$miRExpAssess) <- "miRNA expression measured"
label(dat_table1$adjuvant_IFN) <- "Received adjuvant IFN treatment"

# function to display p-values  
rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- dat_table1[[name]] 
    ind <- !is.na(y)
    y <- y[ind]
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- t.test(y ~ dat_table1$Responder[ind])$p.value
    } else {
      p <- chisq.test(table(y, droplevels(dat_table1$Responder[ind])))$p.value
    }
    s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
    s
  } else {
    render.default(x=x, name=name, ...)
  }
}

rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}

table1(~ Alter + BRAF + Stadium + miRExpAssess + adjuvant_IFN + Hirnmetastase + sex + ECOG + breslow_thickness_mm + subtype + localization | Responder,
       data=dat_table1, droplevels=F, render=rndr, render.strat=rndr.strat)















# tidy miRNA data.....................................................................................................
dat_miRNA_tidy <- dat %>% 
  # only use data where miRNA data was measured 
  filter(miRExpAssess == 1) %>%
  gather(miRNA, expression, contains("hsa")) %>%
  mutate(miRNA = str_replace_all(.$miRNA, "hsa-","")) %>%
  mutate(log_exp = log2(expression))







# Plot miRNA data
plot_miRNA <- signif_plot_Melanoma(dat_miRNA_tidy, x="Responder", y="log_exp", signif=0.05,
                     plot.type = "dotplot", significance=F, Legend = F, var.equal = F,
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





