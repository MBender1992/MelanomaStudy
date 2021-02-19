# <<<<<<<<<<<<< HEAD

# load packages
library(missForest)
library(tidyverse)
library(devtools)
library(caret)
library(doParallel)

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

train.data  <- dat_fct[ind.train, ] # n = 43
test.data <- dat_fct[-ind.train, ] # n = 18

#####################################
##
##  2.a EDA on training set
##
#####################################

# change data structure for ggplot
dat_miR <- train.data %>% 
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

# transform training data
tmp.train <- train.data %>% select(where(is.numeric)) 
fctrs <- train.data %>% select(!where(is.numeric))
train <- data.frame(cbind(log(tmp.train  + 1), fctrs))

# transform testing data
tmp.test  <- test.data %>% select(where(is.numeric))
fctrs <- test.data %>% select(!where(is.numeric))
test <- data.frame(cbind(log(tmp.test  + 1), fctrs))


#####################################
##
## 2.b conversion of factors to dummy variables
##
#####################################

# dummy encoding of factors with model.matrix
x.train <- model.matrix(Responder~.,data=train)[,-1] 
y.train <- train$Responder

# dummy encoding of factors with model.matrix for test set
x.test <- model.matrix(Responder~.,data=test)[,-1] 
y.test <- test$Responder




#####################################
##
## 2.c Fitting the model
##
#####################################

# nested cv??
# https://stackoverflow.com/questions/62276400/how-to-do-nested-cross-validation-with-lasso-in-caret-or-tidymodels
# https://www.tidymodels.org/learn/work/nested-resampling/
# https://stats.stackexchange.com/questions/65128/nested-cross-validation-for-model-selection

# activate parallel computing
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)

# define ctrl function
cctrl1 <- trainControl(method="repeatedcv", number=10,repeats = 5, returnResamp="all", 
                       classProbs=TRUE, summaryFunction=twoClassSummary)

# run glmnet model
set.seed(849)
md <- train(x.train, y.train, method = "glmnet",preProcess = c("center","scale"),
            trControl = cctrl1,metric = "ROC",tuneGrid = expand.grid(alpha = seq(0,1,0.1),
                                                                     lambda = seq(0.001,0.2,by = 0.001)))

max(md$results$ROC)
coef(md$finalModel, md$finalModel$lambdaOpt)






# bootstrap optimized auc
library(nlpred)
x <- model.matrix(Responder~.,data=dat_fct)[,-1] 
y <- ifelse(dat_fct$Responder == "ja", 1,0)
boot_auc(y, x, B = 10, learner = "glm_wrapper")



pred <- md %>% predict(x.test, type = "prob")
obs <- y.test
library(pROC)
roc_obj <- roc(obs, pred$ja)
auc(roc_obj)

# cross validation to find the optimal parameters for elastic net regression


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

# oder LOOCV mit Feature selection process f?r jede LOOCV iteration anders --> anschlie?end Mittel der Features und daraus ein Model bilden? 
# vorher Gedanken machen ob logarithmieren oder nicht




library(glmnet)




dat2 <- dat2%>% select(-c(ECOG,Hirnmetastase, ID, nras, Baseline, subtype, localization, Alter, breslow_thickness_mm,adjuvant_IFN, befallen_Organe))
dat3 <- dat2 %>% filter(!is.na(BRAF) & !is.na(S100)& !is.na(Stadium))

ind <- which(names(dat3) == "Responder")
dat3$Responder <- ifelse(dat3$Responder == "ja", 1,0)


# Lasso 
x <- model.matrix(as.formula(paste("Responder ~.")),data=dat3)
x <- x[,-1]

set.seed(30)
test.cv <- cv.glmnet(x, dat3$Responder,type.measure = "auc",family = "binomial", nfolds = 5,alpha = 1)
plot(test)
coef(test, s = "lambda.1se")


lasso.model <- glmnet(x, dat3$Responder,type.measure = "auc", alpha = 1, family = "binomial",
                      lambda = test.cv$lambda.1se)
coef(lasso.model)







library(doParallel)
# cross validation to find the optimal parameters for elastic net regression
a <- seq(0.1, 0.9, 0.05)
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- cv.glmnet(x, dat3$Responder,type.measure = "auc", family = "binomial", nfold = 10,  parallel = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
}
cv3 <- search[search$cvm == min(search$cvm), ]
md3 <- glmnet(x, dat3$Responder, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
coef(md3)







library(caret)

# Using caret to perform CV
cctrl1 <- trainControl(method="cv", number=10, returnResamp="all",
                       classProbs=TRUE, summaryFunction=twoClassSummary)

x <- model.matrix(Responder ~.,data=dat3)
x <- x[,-1]

set.seed(849)
test_class_cv_model <- train(x, dat3$Responder, method = "glmnet", 
                             trControl = cctrl1,metric = "ROC",tuneGrid = expand.grid(alpha = seq(0,1,0.1),
                                                                                      lambda = seq(0.001,0.2,by = 0.001)))
coef(test_class_cv_model$finalModel, test_class_cv_model$finalModel$lambdaOpt)








x <- model.matrix(as.formula(paste("Responder ~ LDH + BRAF + Eosinophile +`hsa-mir-514a-3p`")),data=dat3)
x <- x[,-1]

set.seed(849)
test_class_cv_model <- train(x, dat3$Responder, method = "xgbTree", 
                             trControl = cctrl1,metric = "ROC")









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




