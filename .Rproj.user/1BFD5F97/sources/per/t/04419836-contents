# <<<<<<<<<<<<< HEAD

# load packages
library(missForest)
library(tidyverse)
library(devtools)
library(caret)
library(doParallel)
library(pROC)
library(pbapply)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

#####################################
## 
## 1.Data loading and preprocessing
##
#####################################

# load data with custom function for melanoma data only for Responders
dat <- load_melanoma_data() %>% 
  filter(!is.na(Responder)) # n = 81

dat_fct <- dat %>%
  filter(miRExpAssess == 1) %>%
  select(-c(TRIM_PDL1_Expression , miRExpAssess, therapy_at_blood_draw)) %>%
  mutate( across(c(Responder, Stadium, Baseline, BRAF, ECOG, subtype, localization,
                   sex, Hirnmetastase, adjuvant_IFN, befallen_Organe, nras), as.factor)) 

xtabs(~ Responder + Stadium, data=dat_fct) 
xtabs(~ Responder + BRAF, data=dat_fct)
xtabs(~ Responder + Baseline, data=dat_fct)
xtabs(~ Responder + ECOG, data=dat_fct) # not enough samples in ECOG2
xtabs(~ Responder + subtype, data=dat_fct) # too many groups with few samples
xtabs(~ Responder + localization, data=dat_fct) # too many groups with few samples
xtabs(~ Responder + sex, data=dat_fct)
xtabs(~ Responder + Hirnmetastase, data=dat_fct)
xtabs(~ Responder + adjuvant_IFN, data=dat_fct)
xtabs(~ Responder + befallen_Organe, data=dat_fct)
xtabs(~ Responder + nras, data=dat_fct) # too few observations

# remove columns that yield high uncertainty
dat_fct$ECOG <- NULL
dat_fct$subtype <- NULL
dat_fct$localization <- NULL
dat_fct$nras <- NULL
dat_fct$Baseline <- NULL
dat_fct$ID <- NULL

#####################################
## 
## 1.a Imputation of missing values
##
#####################################

# detect percentage of NAs in each column
NAs <- sapply(dat_fct, function(df){
  sum(is.na(df) ==TRUE)/length(df);
})

# remove columns with more than 5 % NAs
dat_fct <- dat_fct[, -which(NAs > 0.05)]

# convert factor columns to numerical 
dat_fct$BRAF <- ifelse(dat_fct$BRAF == "pos", 1, 0)
dat_fct$Stadium <- ifelse(dat_fct$Stadium == "II", 2,ifelse(dat_fct$Stadium == "III", 3, 4))

# impute missing values with random forest algorithm
set.seed(25)
dat_imp <- dat_fct %>% 
  select_if(is.numeric) %>%
  as.data.frame() %>%
  missForest() %>%
  .$ximp %>%
  # replace calculated probabilities by the factor
  mutate(BRAF = ifelse(BRAF > 0.5, 1,0),
         Stadium = round(Stadium))

# replace numerical values by factor for encoding later
dat_imp$BRAF <- factor(dat_imp$BRAF, levels = c(0,1), labels = c("neg", "pos"))
dat_imp$Stadium <- factor(dat_imp$Stadium, levels = c(2,3,4), labels = c("II", "III", "IV"))

# replacing NAs with imputed values
dat_fct$BRAF <- dat_imp$BRAF
dat_fct$Stadium <- dat_imp$Stadium
dat_fct$S100 <- dat_imp$S100






#####################################
##
##  2. Modeling process  
##
#####################################

# define test and training set
set.seed(123)
ind.train <- createDataPartition(dat_fct$Responder, p = 0.7, list = FALSE)

train.EDA  <- dat_fct[ind.train, ] # n = 43
test.EDA <- dat_fct[-ind.train, ] # n = 18

#####################################
##
##  2.a EDA on training set (to avoid drawing conclusions including the test set)
##
#####################################

# change data structure for ggplot
dat_miR <- train.EDA %>% 
  select(contains("mir")) %>% 
  gather("miRNA", "expression")

# draw histograms for all miRNAs
miR_hist <- dat_miR %>% 
  ggplot(aes(expression)) +
  geom_histogram() +
  facet_wrap(~miRNA, scales = "free")

# draw qqplots for all miRNAs
miR_qq <- dat_miR %>% 
  ggplot(aes(sample = expression)) +
  geom_qq() +
  geom_qq_line() +
  facet_wrap(~miRNA, scales = "free")



# draw histograms for all miRNAs log-transformed
miR_hist_log <- dat_miR %>% 
  ggplot(aes(log(expression))) +
  geom_histogram() +
  facet_wrap(~miRNA, scales = "free")

# draw qqplots for all miRNAs log-transformed
miR_qq_log <- dat_miR %>% 
  ggplot(aes(sample = log(expression))) +
  geom_qq() +
  geom_qq_line() +
  facet_wrap(~miRNA, scales = "free")

## log-transforming miRNA expression improves approximation to normality and gene expression data is known to be log-normal distributed 
## log-transformed miRNA values were used for ML

# histogram of numerical variables that are not miRNAs
par(mfrow = c(2,4))
par(mar=c(0.5, 4.5, 0.5, 0.5))

# original expression values
hist(train.data$LDH)
hist(train.data$Eosinophile)
hist(train.data$S100)
hist(train.data$CRP)

# log-transformed expression values
hist(log(train.data$LDH))
hist(log(train.data$Eosinophile))
hist(log(train.data$S100))
hist(log(train.data$CRP))

## lab parameters were also used in log-transformed space

#####################################
##
## 2.b log-transformation and conversion of factors to dummy variables
##
#####################################

# transform the whole dataset
tmp <- dat_fct %>% select(where(is.numeric))
fctrs <- dat_fct %>% select(!where(is.numeric))
dat_log <- data.frame(cbind(log(tmp+1), fctrs))


#
y <- dat_log$Responder







#####################################
##
## 2.c Fitting different models
##
#####################################


#####################################
##
## c.1 complete model
##
#####################################


# define parameters for 10 fold cross validation repeated 10 times
k <- 10
rep <- 10
# models.lasso.complete <- lassoEval("complete", dat_log, rep = rep, k = k)
models.lasso.complete <- readRDS("models/models_lasso_complete.rds")

# set names of list elements
models.lasso.complete <- setNames(lapply(models.lasso.complete, setNames, folds), reps)

# extract metrics for inner fold (training) and outer fold (testing) from list and convert to df
df.train.complete <- trainDF(models.lasso.complete)
df.test.complete <- testDF(models.lasso.complete)

# extract important coefficients
extract.coefs.complete <- extractCoefs(models.lasso.complete) %>% do.call(rbind,.) %>% table() 


# calculate percentages
feat.freq.complete <- data.frame(sort(extract.coefs.complete/100)) %>% 
  setNames(c("coef", "freq"))

# plot important features
ggplot(data = feat.freq.complete, aes(coef, freq, fill = ifelse(freq > 0.5, "red", "blue"))) +
  geom_bar(stat = "identity",  color = "black") + 
  coord_flip() +
  xlab("") +
  ylab("fraction of cv-models using this feature (relative feature importance)") +
  theme_bw() +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), expand = c(0,0), labels = scales::percent_format()) +
  geom_hline(yintercept = 0.5, lty = 2, color = "red") + 
  scale_fill_manual(labels = c("< 50 %", "> 50 %"), values = c("gray95", "lightblue")) +
  labs(fill = "frequency")





#####################################
##
## c.2 baseline model
##
#####################################

# model process and evaluation, k and rep define fold and repeats in outer loop 
# models.lasso.baseline <- lassoEval("baseline", dat_log, rep = rep, k = k)
models.lasso.baseline <- readRDS("models/models_lasso_baseline.rds")

# set names of list elements
models.lasso.baseline <- setNames(lapply(models.lasso.baseline, setNames, folds), reps)

# extract metrics for inner fold (training) and outer fold (testing) from list and convert to df
df.train.baseline <- trainDF(models.lasso.baseline)
df.test.baseline <- testDF(models.lasso.baseline)

# extract important coefficients
extract.coefs.baseline <- extractCoefs(models.lasso.baseline) %>% do.call(rbind,.) %>% table() 


# calculate percentages
feat.freq <- data.frame(sort(extract.coefs.baseline/100)) %>% 
  setNames(c("coef", "freq"))

# plot important features
ggplot(data = feat.freq, aes(coef, freq, fill = ifelse(freq > 0.5, "red", "blue"))) +
  geom_bar(stat = "identity",  color = "black") + 
  coord_flip() +
  xlab("") +
  ylab("fraction of cv-models using this feature (relative feature importance)") +
  theme_bw() +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), expand = c(0,0), labels = scales::percent_format()) +
  geom_hline(yintercept = 0.5, lty = 2, color = "red") + 
  scale_fill_manual(labels = c("< 50 %", "> 50 %"), values = c("gray95", "lightblue")) +
  labs(fill = "frequency")









#####################################
##
## c.3 signif
##
#####################################

# 
# models.lasso.signif <- lassoEval("signif", dat_log, rep = rep, k = k)
models.lasso.signif <- readRDS("models/models_lasso_signif.rds")

# set names of list elements
models.lasso.signif <- setNames(lapply(models.lasso.signif, setNames, folds), reps)

# extract metrics for inner fold (training) and outer fold (testing) from list and convert to df
df.train.signif <- trainDF(models.lasso.signif)
df.test.signif <- testDF(models.lasso.signif)

# extract important coefficients
extract.coefs.signif <- extractCoefs(models.lasso.signif) %>% do.call(rbind,.) %>% table() 


# calculate percentages
feat.freq <- data.frame(sort(extract.coefs.signif/100)) %>% 
  setNames(c("coef", "freq"))

# plot important features
ggplot(data = feat.freq, aes(coef, freq)) +
  geom_bar(stat = "identity",  color = "black", fill = "lightblue") + 
  coord_flip() +
  xlab("") +
  ylab("fraction of cv-models using this feature (relative feature importance)") +
  theme_bw() +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), expand = c(0,0), labels = scales::percent_format()) +
  geom_hline(yintercept = 0.5, lty = 2, color = "red") + 
  labs(fill = "frequency")








#####################################
##
## c.4 miRNA
##
#####################################

#
# models.lasso.miRNA <- lassoEval("miRNA", dat_log, rep = 10, k = 10)
models.lasso.miRNA <- readRDS("models/models_lasso_miRNA.rds")

# set names of list elements
models.lasso.miRNA <- setNames(lapply(models.lasso.miRNA, setNames, folds), reps)

# extract metrics for inner fold (training) and outer fold (testing) from list and convert to df
df.train.miRNA <- trainDF(models.lasso.miRNA)
df.test.miRNA <- testDF(models.lasso.miRNA)

# extract important coefficients
extract.coefs.miRNA <- extractCoefs(models.lasso.miRNA) %>% do.call(rbind,.) %>% table() 


# calculate percentages
feat.freq <- data.frame(sort(extract.coefs.miRNA/100)) %>% 
  setNames(c("coef", "freq"))

# plot important features
ggplot(data = feat.freq, aes(coef, freq, fill = ifelse(freq > 0.5, "red", "blue"))) +
  geom_bar(stat = "identity",  color = "black") + 
  coord_flip() +
  xlab("") +
  ylab("fraction of cv-models using this feature (relative feature importance)") +
  theme_bw() +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), expand = c(0,0), labels = scales::percent_format()) +
  geom_hline(yintercept = 0.5, lty = 2, color = "red") + 
  scale_fill_manual(labels = c("< 50 %", "> 50 %"), values = c("gray95", "lightblue")) +
  labs(fill = "frequency")






#####################################
##
## c.5 relaxed LASSO 
##
#####################################

# obtain features for relaxed LASSO analysis (features with importance > 0.5, BRAF added manually within the function)
feat.relaxed <-  feat.freq.complete[feat.freq.complete$freq > 0.5,]
feat.relaxed <- feat.relaxed[feat.relaxed$coef != "BRAFpos",]

# modelling and evaluation
# models.lasso.relaxedLasso <- lassoEval("relaxedLasso", dat_log, rep = rep, k = k)
models.lasso.relaxedLasso <- readRDS("models/models_lasso_relaxedLasso.rds")

# set names of list elements
models.lasso.relaxedLasso <- setNames(lapply(models.lasso.relaxedLasso, setNames, folds), reps)

# extract metrics for inner fold (training) and outer fold (testing) from list and convert to df
df.train.relaxedLasso <- trainDF(models.lasso.relaxedLasso)
df.test.relaxedLasso <- testDF(models.lasso.relaxedLasso)

# extract important coefficients
extract.coefs.relaxedLasso <- extractCoefs(models.lasso.relaxedLasso) %>% do.call(rbind,.) %>% table() 

# calculate percentages
feat.freq <- data.frame(sort(extract.coefs.relaxedLasso/100)) %>% 
  setNames(c("coef", "freq"))

# plot important features
ggplot(data = feat.freq, aes(coef, freq)) +
  geom_bar(stat = "identity",  color = "black", fill = "lightblue") + 
  coord_flip() +
  xlab("") +
  ylab("fraction of cv-models using this feature (relative feature importance)") +
  theme_bw() +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), expand = c(0,0), labels = scales::percent_format()) +
  geom_hline(yintercept = 0.5, lty = 2, color = "red") 



## to do
# ROC Kurve
# Model process strukturieren
# konfidenzintervalle konstruieren
# relaxed miRNA model?
# feature selection
  # nur serum parameter
  # nur miRNAs
  # signifikante features
  # complete model with lasso for feature selection and subsequent relaxed lasso







## inner loop: 10-fold cv repeated 5 times, outer loop: 10-fold cv repeated 10 times
## An Analysis on Better Testing than Training Performances on the Iris Dataset
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7565855/
## Modelling Process angucken
## calibration curve





# cross validation to find the optimal parameters for elastic net regression

# nested cv??
# https://stackoverflow.com/questions/62276400/how-to-do-nested-cross-validation-with-lasso-in-caret-or-tidymodels
# https://www.tidymodels.org/learn/work/nested-resampling/
# https://stats.stackexchange.com/questions/65128/nested-cross-validation-for-model-selection

# internal validation via bootstrap optimism? https://stats.stackexchange.com/questions/103411/internal-validation-via-bootstrap-what-roc-curve-to-present
# https://stats.stackexchange.com/questions/61344/getting-the-bootstrap-validated-auc-in-r
# https://www.rdocumentation.org/packages/rms/versions/6.1-1/topics/lrm mal durcharbeiten
# concordance measure anstatt ROC?
# cost function?

# 1 model mit signif
# 1 model mit LASSO feature selection
# 1 model nur mit miRNas
# 1 model nur mit herkömmlichen Pr?diktoren

# define test and training set
# use glmnet for feature selection
# train glm to get coefficients and quantitative probabilities of response 
# think about a proper cross-validation loop 


# As a side comment, if you want to interpret the result be sure to demonstrate the that set of
# variables selected by lasso is stable. This can be done using Monte Carlo simulation or by bootstrapping 
# your own dataset.

#Bootstrapping is the process of resampling with replacement (all values in the sample have an equal probability of being selected, 
# including multiple times, so a value could have a duplicate).

# glmnet for feature selection and then rf or xgbTree? or just a simple glm to have coefficients? 
# or run another glmnet with the reduced variables

# https://stats.stackexchange.com/questions/25305/is-this-a-correct-procedure-for-feature-selection-using-cross-validation?rq=1
# https://stats.stackexchange.com/questions/60692/not-all-features-selected-by-glmnet-considered-signficant-by-glm-logistic-regre

# run on bootstrap samples not just different data splits of test and training


# Unterteilung in Training und Test?
# inner loop for feature selection im Trainingsset?
# anschlie?end glm model auf ganzes set mit den features, um model coefficients zu erhalten 
# independent test auf Testset?







# simple logistic regression (Not recommended when features have been chosen by LASSO or elastic net regularization)
logistic <- glm(Responder ~ LDH + BRAF + Eosinophile +`hsa-mir-514a-3p`,data=dat3, family="binomial")
summary(logistic)


ll.null <- logistic$null.deviance/-2
ll.proposed <- logistic$deviance/-2

## McFadden's Pseudo R^2 = [ LL(Null) - LL(Proposed) ] / LL(Null)
(ll.null - ll.proposed) / ll.null


## The p-value for the R^2
1 - pchisq(2*(ll.proposed - ll.null), df=(length(logistic$coefficients)-1))

nullmodel <- glm(Responder ~ 1, data = dat3, family = "binomial")
anova(nullmodel,logistic, test = 'LRT')




