# <<<<<<<<<<<<< HEAD

# load packages
library(missForest)
library(tidyverse)
library(devtools)
library(caret)
library(pROC)
library(DescTools)
library(pbapply)
library(cvAUC)
library(doParallel)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

#####################################
## 
## 1.Data loading and preprocessing
##
#####################################

# load data with custom function for melanoma data only for Responders
dat <- load_melanoma_data()

dat_fct <- dat %>%
  filter(miRExpAssess == 1 & !is.na(Responder)) %>%
  select(-c(TRIM_PDL1_Expression , miRExpAssess, therapy_at_blood_draw)) %>%
  mutate( across(c(Responder, Stadium, Baseline, BRAF, ECOG, subtype, localization,
                   sex, brainMet, adjuvant_IFN, organsInvolved, nras, prior_BRAF_therapy), as.factor)) 

xtabs(~ Responder + Stadium, data=dat_fct) 
xtabs(~ Responder + BRAF, data=dat_fct)
xtabs(~ Responder + Baseline, data=dat_fct)
xtabs(~ Responder + ECOG, data=dat_fct) # not enough samples in ECOG2
xtabs(~ Responder + subtype, data=dat_fct) # too many groups with few samples
xtabs(~ Responder + localization, data=dat_fct) # too many groups with few samples
xtabs(~ Responder + sex, data=dat_fct)
xtabs(~ Responder + brainMet, data=dat_fct)
xtabs(~ Responder + adjuvant_IFN, data=dat_fct)
xtabs(~ Responder + organsInvolved, data=dat_fct) # too few observations
xtabs(~ Responder + nras, data=dat_fct) # too few observations
xtabs(~ Responder + prior_BRAF_therapy, data=dat_fct)



# remove columns that yield high uncertainty
dat_fct$ECOG <- NULL
dat_fct$subtype <- NULL
dat_fct$localization <- NULL
dat_fct$nras <- NULL
dat_fct$Baseline <- NULL
dat_fct$ID <- NULL
dat_fct$organsInvolved <- NULL
dat_fct$breslow_thickness_mm <- NULL # included in stage variable
dat_fct$BRAF <- NULL # highly correlated with prior BRAF therapy and therefore rather adds noise to the model


#####################################
## 
## 1.a Impute missing values
##
#####################################

# detect percentage of NAs in each column
NAs <- sapply(dat_fct, function(df){
  sum(is.na(df) ==TRUE)/length(df);
})

# remove columns with more than 20 % NAs
dat_fct <- dat_fct[, -which(NAs > 0.2)]

# convert factor columns to numerical 
dat_fct$Stadium <- ifelse(dat_fct$Stadium == "II", 2,ifelse(dat_fct$Stadium == "III", 3, 4))
dat_fct$sex <- ifelse(dat_fct$sex == "m", 1, 0)
dat_fct$brainMet <- ifelse(dat_fct$brainMet == "ja", 1, 0)
dat_fct$prior_BRAF_therapy <- parse_number(as.character(dat_fct$prior_BRAF_therapy))

# impute missing values with random forest algorithm
set.seed(25)
dat_imp <- dat_fct %>% 
  select_if(is.numeric) %>%
  as.data.frame() %>%
  missForest() %>%
  .$ximp %>%
  # replace calculated probabilities by the factor
  mutate(Stadium = round(Stadium), 
         Alter = round(Alter), 
         brainMet = ifelse(brainMet > 0.5, 1,0), 
         prior_BRAF_therapy = ifelse(prior_BRAF_therapy > 0.5, 1, 0))



# replace numerical values by factor for encoding later
dat_imp$sex <- factor(dat_imp$sex, levels = c(0,1), labels = c("w", "m"))
dat_imp$brainMet <- factor(dat_imp$brainMet, levels = c(0,1), labels = c("no", "yes"))
dat_imp$prior_BRAF_therapy  <- factor(dat_imp$prior_BRAF_therapy, levels = c(0,1), labels = c("no", "yes"))



# replacing NAs with imputed values
dat_fct$Stadium <- dat_imp$Stadium
dat_fct$S100 <- dat_imp$S100
dat_fct$Alter<- dat_imp$Alter
dat_fct$brainMet <- dat_imp$brainMet
dat_fct$prior_BRAF_therapy <- dat_imp$prior_BRAF_therapy
dat_fct$sex <- dat_imp$sex


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
  geom_histogram(color = "black", fill = "grey") +
  facet_wrap(~miRNA, scales = "free") + 
  theme_bw()

# png("miRNA_histogram.png", units="in", width=12, height=8, res=1200)
# miR_hist
# dev.off()


# draw qqplots for all miRNAs
miR_qq <- dat_miR %>% 
  ggplot(aes(sample = expression)) +
  geom_qq() +
  geom_qq_line() +
  facet_wrap(~miRNA, scales = "free")+
  theme_bw()

# png("miRNA_qq.png", units="in", width=12, height=8, res=1200)
# miR_qq
# dev.off()


# draw histograms for all miRNAs log-transformed
miR_hist_log <- dat_miR %>% 
  ggplot(aes(log(expression))) +
  geom_histogram(color = "black", fill = "grey") +
  facet_wrap(~miRNA, scales = "free") + 
  theme_bw()

# png("miRNA_histogram_log.png", units="in", width=12, height=8, res=1200)
# miR_hist_log
# dev.off()


# draw qqplots for all miRNAs log-transformed
miR_qq_log <- dat_miR %>% 
  ggplot(aes(sample = log(expression))) +
  geom_qq() +
  geom_qq_line() +
  facet_wrap(~miRNA, scales = "free") +
  theme_bw()

# png("miRNA_qq_log.png", units="in", width=12, height=8, res=1200)
# miR_qq_log
# dev.off()


## log-transforming miRNA expression improves approximation to normality and gene expression data is known to be log-normal distributed 
## log-transformed miRNA values were used for ML

# histogram of numerical variables that are not miRNAs
par(mfrow = c(2,4))
par(mar=c(0.5, 4.5, 0.5, 0.5))

# original expression values
hist(train.data$LDH, main = "LDH")
hist(train.data$Eosinophile, main = "Eosinophile")
hist(train.data$S100, main = "S100")
hist(train.data$CRP, main = "CRP")

# log-transformed expression values
hist(log(train.data$LDH), main = "log-transformed LDH")
hist(log(train.data$Eosinophile), main = "log-transformed Eosinophile")
hist(log(train.data$S100), main = "log-transformed S100")
hist(log(train.data$CRP), main = "log-transformed CRP")

## lab parameters were also used in log-transformed space

#####################################
##
## 2.b log-transformation 
##
#####################################

# transform the whole dataset
tmp <- dat_fct %>% select(where(is.numeric))
fctrs <- dat_fct %>% select(!where(is.numeric))
dat_log <- data.frame(cbind(log(tmp+1), fctrs)) 



#####################################
##
## 2.c Fitting different models
##
#####################################


#####################################
##
## c.1a complete model
##
#####################################

#
y <- dat_log$Responder

# define parameters for 10 fold cross validation repeated 10 times (outer loop)
k <- 10
rep <- 10

reps <- paste0("Rep", 1:rep)
folds <- paste0("Fold", 1:k)

model.matrix.complete <- model.matrix.subset("complete", data = dat_log)
models.lasso.complete <- mlEval(model.matrix.complete, dat_log, rep = rep, k = k)
# saveRDS(models.lasso.complete, "models/models_lasso_complete.rds")
# models.lasso.complete <- readRDS("models/models_lasso_complete.rds")


# set names of list elements
models.lasso.complete <- setNames(lapply(models.lasso.complete, setNames, folds), reps)

## confidence interval for the cv.train folds in the inner loop 
ci.complete <- rbind.model.ci(models.lasso.complete)

# extract important coefficients
extract.coefs.complete <- extractCoefs(models.lasso.complete) %>% do.call(rbind,.) %>% table() 

# calculate percentages
feat.freq.complete <- data.frame(sort(extract.coefs.complete/100)) %>% 
  setNames(c("coef", "freq"))

# plot important features
# ggplot(data = feat.freq.complete, aes(coef, freq, fill = ifelse(freq > 0.5, "red", "blue"))) +
#   geom_bar(stat = "identity",  color = "black") +
#   coord_flip() +
#   xlab("") +
#   ylab("fraction of cv-models using this feature (relative feature importance)") +
#   theme_bw() +
#   scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), expand = c(0,0), labels = scales::percent_format()) +
#   geom_hline(yintercept = 0.5, lty = 2, color = "red") +
#   scale_fill_manual(labels = c("< 50 %", "> 50 %"), values = c("gray95", "lightblue")) +
#   labs(fill = "frequency")



#####################################
##
## c.1b relaxed LASSO of the complete model
##
#####################################

# obtain features for relaxed LASSO analysis (features with importance > 0.5, BRAF added manually within the function)
feat.relaxed <-  feat.freq.complete[feat.freq.complete$freq > 0.5,]
feat.relaxed <- feat.relaxed[as.character(feat.relaxed$coef) %like any% names(dat_log),]


# modelling and evaluation
model.matrix.relaxedLasso <- model.matrix.subset("relaxedLasso", data = dat_log)
models.lasso.relaxedLasso <- mlEval(model.matrix.relaxedLasso, dat_log, rep = rep, k = k)
# saveRDS(models.lasso.relaxedLasso, "models/models_lasso_relaxedLasso.rds")
# models.lasso.relaxedLasso <- readRDS("models/models_lasso_relaxedLasso.rds")

# set names of list elements
models.lasso.relaxedLasso <- setNames(lapply(models.lasso.relaxedLasso, setNames, folds), reps)

## confidence interval for the cv.train folds in the inner loop 
ci.relaxedLasso <- rbind.model.ci(models.lasso.relaxedLasso)

# extract important coefficients
extract.coefs.relaxedLasso <- extractCoefs(models.lasso.relaxedLasso) %>% do.call(rbind,.) %>% table() 

# calculate percentages
feat.freq <- data.frame(sort(extract.coefs.relaxedLasso/100)) %>% 
  setNames(c("coef", "freq"))

# plot important features
# ggplot(data = feat.freq, aes(coef, freq)) +
#   geom_bar(stat = "identity",  color = "black", fill = "lightblue") + 
#   coord_flip() +
#   xlab("") +
#   ylab("fraction of cv-models using this feature (relative feature importance)") +
#   theme_bw() +
#   scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), expand = c(0,0), labels = scales::percent_format()) +
#   geom_hline(yintercept = 0.5, lty = 2, color = "red") 





#####################################
##
## c.2 baseline model
##
#####################################

# model process and evaluation, k and rep define fold and repeats in outer loop 
model.matrix.baseline <- model.matrix.subset("baseline", data = dat_log)
models.lasso.baseline <- mlEval(model.matrix.baseline, dat_log, rep = rep, k = k)
# saveRDS(models.lasso.baseline, "models/models_lasso_baseline.rds")
# models.lasso.baseline <- readRDS("models/models_lasso_baseline.rds")

# set names of list elements
models.lasso.baseline <- setNames(lapply(models.lasso.baseline, setNames, folds), reps)

## confidence interval for the cv.train folds in the inner loop 
ci.baseline <- rbind.model.ci(models.lasso.baseline)

# extract important coefficients
extract.coefs.baseline <- extractCoefs(models.lasso.baseline) %>% do.call(rbind,.) %>% table() 

# calculate percentages
feat.freq <- data.frame(sort(extract.coefs.baseline/100)) %>% 
  setNames(c("coef", "freq"))

# plot important features
# ggplot(data = feat.freq, aes(coef, freq, fill = ifelse(freq > 0.5, "red", "blue"))) +
#   geom_bar(stat = "identity",  color = "black") + 
#   coord_flip() +
#   xlab("") +
#   ylab("fraction of cv-models using this feature (relative feature importance)") +
#   theme_bw() +
#   scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), expand = c(0,0), labels = scales::percent_format()) +
#   geom_hline(yintercept = 0.5, lty = 2, color = "red") + 
#   scale_fill_manual(labels = c("< 50 %", "> 50 %"), values = c("gray95", "lightblue")) +
#   labs(fill = "frequency")






#####################################
##
## c.3 signif
##
#####################################

# 
model.matrix.signif <- model.matrix.subset("signif", data = dat_log)
models.lasso.signif <- mlEval(modelM, dat_log, rep = rep, k = k)
# saveRDS(models.lasso.signif, "models/models_lasso_signif.rds")
# models.lasso.signif <- readRDS("models/models_lasso_signif.rds")

# set names of list elements
models.lasso.signif <- setNames(lapply(models.lasso.signif, setNames, folds), reps)

## confidence interval for the cv.train folds in the inner loop 
ci.signif <- rbind.model.ci(models.lasso.signif)

# extract important coefficients
extract.coefs.signif <- extractCoefs(models.lasso.signif) %>% do.call(rbind,.) %>% table() 


# calculate percentages
feat.freq <- data.frame(sort(extract.coefs.signif/100)) %>% 
  setNames(c("coef", "freq"))

# plot important features
# ggplot(data = feat.freq, aes(coef, freq)) +
#   geom_bar(stat = "identity",  color = "black", fill = "lightblue") +
#   coord_flip() +
#   xlab("") +
#   ylab("fraction of cv-models using this feature (relative feature importance)") +
#   theme_bw() +
#   scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), expand = c(0,0), labels = scales::percent_format()) +
#   geom_hline(yintercept = 0.5, lty = 2, color = "red") +
#   labs(fill = "frequency")





#####################################
##
## c.4a miRNA
##
#####################################

#
model.matrix.miRNA <- model.matrix.subset("miRNA", data = dat_log)
models.lasso.miRNA <- mlEval(model.matrix.miRNA, dat_log, rep = rep, k = k)
# saveRDS(models.lasso.miRNA, "models/models_lasso_miRNA.rds")
# models.lasso.miRNA <- readRDS("models/models_lasso_miRNA.rds")

# set names of list elements
models.lasso.miRNA <- setNames(lapply(models.lasso.miRNA, setNames, folds), reps)

## confidence interval for the cv.train folds in the inner loop 
ci.miRNA <- rbind.model.ci(models.lasso.miRNA)

# extract important coefficients
extract.coefs.miRNA <- extractCoefs(models.lasso.miRNA) %>% do.call(rbind,.) %>% table() 

# calculate percentages
feat.freq.miRNA <- data.frame(sort(extract.coefs.miRNA/100)) %>% 
  setNames(c("coef", "freq"))

# plot important features
# ggplot(data = feat.freq.miRNA, aes(coef, freq, fill = ifelse(freq > 0.5, "red", "blue"))) +
#   geom_bar(stat = "identity",  color = "black") +
#   coord_flip() +
#   xlab("") +
#   ylab("fraction of cv-models using this feature (relative feature importance)") +
#   theme_bw() +
#   scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), expand = c(0,0), labels = scales::percent_format()) +
#   geom_hline(yintercept = 0.5, lty = 2, color = "red") +
#   scale_fill_manual(labels = c("< 50 %", "> 50 %"), values = c("gray95", "lightblue")) +
#   labs(fill = "frequency")



#####################################
##
## c.4b relaxed LASSO miRNA
##
#####################################

feat.relaxed.miRNA <-  feat.freq.miRNA[feat.freq.miRNA$freq > 0.5,]

model.matrix.relaxed.miRNA <- model.matrix.subset("relaxedLassomiRNA", data = dat_log)
models.lasso.relaxed.miRNA <- mlEval(model.matrix.relaxed.miRNA, dat_log, rep = 10, k = 10)
#saveRDS(models.lasso.relaxed.miRNA, "models/models_lasso_relaxed_miRNA.rds")
# models.lasso.relaxed.miRNA <- readRDS("models/models_lasso_relaxed_miRNA.rds")

# set names of list elements
models.lasso.relaxed.miRNA <- setNames(lapply(models.lasso.relaxed.miRNA, setNames, folds), reps)

## confidence interval for the cv.train folds in the inner loop 
ci.relaxed.miRNA <- rbind.model.ci(models.lasso.relaxed.miRNA)

# extract important coefficients
extract.coefs.relaxed.miRNA <- extractCoefs(models.lasso.relaxed.miRNA) %>% do.call(rbind,.) %>% table() 

# calculate percentages
feat.freq.relaxed.miRNA <- data.frame(sort(extract.coefs.relaxed.miRNA/100)) %>% 
  setNames(c("coef", "freq"))

# plot important features
ggplot(data = feat.freq.relaxed.miRNA, aes(coef, freq)) +
  geom_bar(stat = "identity",  color = "black", fill = "skyblue") + 
  coord_flip() +
  xlab("") +
  ylab("fraction of cv-models using this feature (relative feature importance)") +
  theme_bw() +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), expand = c(0,0), labels = scales::percent_format()) +
  geom_hline(yintercept = 0.5, lty = 2, color = "red") 
  

#####################################
##
## d. Model comparison
##
#####################################

# combine data to compare models
dat_compare <- rbind(complete = ci.complete,
      relaxedLasso = ci.relaxedLasso,
      baseline = ci.baseline,
      signif = ci.signif,
      miRNA = ci.miRNA,
      relaxedmiRNA = ci.relaxed.miRNA) %>% 
  rownames_to_column("tmp") %>%
  separate(tmp,c("model", "results"), extra = "merge") %>%
  mutate(model = factor(model),
         model = reorder(model, cvAUC))

# train inner cv ROC
ggplot(filter(dat_compare, results == "cv.AUC.inner"), aes(x=model, y=cvAUC)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width = 0.3, size = 1) +
  geom_point(size = 4, shape = 18, color = "red") +
  coord_flip() + 
  theme_bw() + 
  scale_y_continuous(breaks = seq(0.5, 0.9, 0.05))+
  ylab("ROC")


# comparison of train.outer and test.outer
ggplot(dat_compare, aes(x=model, y=cvAUC, color= results)) + 
  geom_point(position=position_dodge(.5)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width = 0.2,position=position_dodge(.5)) + 
  coord_flip() + 
  theme_bw() + 
  scale_y_continuous(breaks = seq(0.5, 0.9, 0.05))+
  ylab("ROC")


## to do
# ROC Kurve
# ci for model coefficients to assess stability of the best model? 


#####################################
##
## e. ROC curve
##
#####################################



# calcuate AUC for different folds
ls <- ls_cvAUC(models.lasso.relaxedLasso)
out <- cvAUC(ls$predictions, ls$labels)
res <- ci.cvAUC(ls$predictions, ls$labels)

#Plot CV AUC
plot(out$perf, col="grey82", lty=3, main="10-fold CV AUC (repeated 10 times)")
plot(out$perf, col="blue", avg="vertical",add =T)
abline(0,1, col = "red", lty = 2)
text(0.8, 0.2, paste("AUC: ",round(res$cvAUC,3), " (",round(res$ci[1],3),"; ",round(res$ci[2],3),")", sep =""))


#####################################
##
## f. final model
##
#####################################

feat.final <- names(select(dat_log, c(feat.relaxed$coef, Alter, prior_BRAF_therapy)))
model.formula <- as.formula(paste("Responder~",paste(feat.final, collapse ="+")))

x <- model.matrix(model.formula, dat_log)[,-1]
y <- dat_log$Responder

set.seed(27)
# 
final <- train(x, y, method = "glmnet",preProcess = c("center","scale"), 
      trControl = cctrl1, metric = "ROC", tuneGrid = expand.grid(alpha = 1, lambda = 0.01))

# ROC of final model
final$results[final$results$lambda == final$finalModel$lambdaOpt,]

# coefficients of final model 
coef(final$finalModel, final$finalModel$lambdaOpt)





test2 <- bind_rows(lapply(1:10, function(x){
    tmp <- bind_rows(sapply(models.lasso.relaxedLasso[[x]], '[', 'coefficients'))
  }))


test2 %>% group_by(coefs) %>% summarize(mean = mean(vals), ci = confInt(vals)) %>%
  mutate(lower = mean - ci, upper = mean + ci) %>% select(-ci)
 

unlist.model(models.lasso.complete, "lambda", "train.metrics") %>% table()



# simple logistic regression (Not recommended when features have been chosen by LASSO or elastic net regularization)
logistic <- glm(Responder ~ LDH + prior_BRAF_therapy + Eosinophile +hsa.mir.514a.3p,data=dat_log, family="binomial")
summary(logistic)



ll.null <- logistic$null.deviance/-2
ll.proposed <- logistic$deviance/-2

## McFadden's Pseudo R^2 = [ LL(Null) - LL(Proposed) ] / LL(Null)
(ll.null - ll.proposed) / ll.null


## The p-value for the R^2
1 - pchisq(2*(ll.proposed - ll.null), df=(length(logistic$coefficients)-1))

nullmodel <- glm(Responder ~ 1, data = dat3, family = "binomial")
anova(nullmodel,logistic, test = 'LRT')


