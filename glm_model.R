# <<<<<<<<<<<<< HEAD

# load packages
library(tidyverse)
library(devtools)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

# load data with custom function for melanoma data only for Responders
dat <- load_melanoma_data() %>% 
  filter(!is.na(Responder)) %>%# n = 81
  mutate(Responder = ifelse(Responder == "ja", 1, 0))





