# <<<<<<<<<<<<< HEAD

# load packages
library(tidyverse)
library(devtools)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

# load data with custom function for melanoma data only for Responders
dat <- load_melanoma_data() %>% 
  filter(!is.na(Responder)) # n = 81

dat2 <- dat %>% dplyr::select(-c(TRIM_PDL1_Expression , miRExpAssess, therapy_at_blood_draw)) %>%
  mutate( across(c(Responder, Stadium,BRAF, Baseline,  ECOG, subtype, localization,
                   sex, Hirnmetastase, adjuvant_IFN, befallen_Organe, nras), as.factor))

xtabs(~ Responder + Stadium, data=dat2) 
xtabs(~ Responder + BRAF, data=dat2)
xtabs(~ Responder + Baseline, data=dat2)
xtabs(~ Responder + ECOG, data=dat2)
xtabs(~ Responder + subtype, data=dat2)
xtabs(~ Responder + localization, data=dat2)
xtabs(~ Responder + sex, data=dat2)
xtabs(~ Responder + Hirnmetastase, data=dat2)
xtabs(~ Responder + adjuvant_IFN, data=dat2)
xtabs(~ Responder + befallen_Organe, data=dat2)
xtabs(~ Responder + nras, data=dat2)


dat_model <- dat2 %>% dplyr::select(c(Responder, breslow_thickness_mm, sex, Hirnmetastase, Alter, adjuvant_IFN, Eosinophile, CRP, LDH, S100, BRAF,`hsa-mir-137`, `hsa-mir-514a-3p` ))

logistic <- glm(Responder ~ log(`hsa-mir-137`) + log(LDH),data=dat_model, family="binomial")
summary(logistic)

# Unterteilung in Training und Test?
# inner loop for feature selection im Trainingsset?
# anschlieﬂend glm model auf ganzes set mit den features, um model coefficients zu erhalten 
# independent test auf Testset?

# oder LOOCV mit Feature selection process f¸r jede LOOCV iteration anders --> anschlieﬂend Mittel der Features und daraus ein Model bilden? 
# vorher Gedanken machen ob logarithmieren oder nicht


  
  
