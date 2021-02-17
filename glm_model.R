# <<<<<<<<<<<<< HEAD

# load packages
library(tidyverse)
library(devtools)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

# load data with custom function for melanoma data only for Responders
dat <- load_melanoma_data() %>% 
  filter(!is.na(Responder)) # n = 81

dat2 <- dat %>%
  filter(miRExpAssess == 1) %>%
  select(-c(TRIM_PDL1_Expression , miRExpAssess, therapy_at_blood_draw)) %>%
  mutate( across(c(Responder, Stadium,BRAF, Baseline,  ECOG, subtype, localization,
                   sex, Hirnmetastase, adjuvant_IFN, befallen_Organe, nras), as.factor)) 
  
  

xtabs(~ Responder + Stadium, data=dat2) 
xtabs(~ Responder + BRAF, data=dat2)
xtabs(~ Responder + Baseline, data=dat2)
xtabs(~ Responder + ECOG, data=dat2) # not enough samples in ECOG2
xtabs(~ Responder + subtype, data=dat2)
xtabs(~ Responder + localization, data=dat2)
xtabs(~ Responder + sex, data=dat2)
xtabs(~ Responder + Hirnmetastase, data=dat2)
xtabs(~ Responder + adjuvant_IFN, data=dat2)
xtabs(~ Responder + befallen_Organe, data=dat2)
xtabs(~ Responder + nras, data=dat2) # too few observations





# remove columns like ECOG and nras with too few observations
# impute missing values for columns with less than 5 % NA, omit other columns
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









predict

library(glmnet)

df_test <- dat2 %>% filter(!is.na(S100)) %>%
  select(LDH, S100, CRP, `hsa-mir-132-3p`, `hsa-mir-137`,`hsa-mir-197-3p`,`hsa-mir-214-3p`, `hsa-mir-514a-3p`, Responder)

test <- cv.glmnet(as.matrix(df_test[,1:8]), df_test$Responder,type.measure = "deviance",family = "binomial", nfolds = 5)



dat2 <- dat2%>% select(-c(ECOG,Hirnmetastase, ID, nras, Baseline, subtype, localization, Alter, breslow_thickness_mm,adjuvant_IFN, befallen_Organe))

dat3 <- dat2 %>% filter(!is.na(BRAF) & !is.na(S100)& !is.na(Stadium))

ind <- which(names(dat3) == "Responder")
dat3$Responder <- ifelse(dat3$Responder == "ja", 1,0)

x <- model.matrix(as.formula(paste("Responder ~.")),data=dat3)
x <- x[,-1]

set.seed(30)
test <- cv.glmnet(x, dat3$Responder,type.measure = "auc",family = "binomial", nfolds = 5,alpha = 0.6)
coef(test, s = "lambda.1se")

library(doParallel)
# cross validation to find the optimal parameters
a <- seq(0.1, 0.9, 0.05)
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- cv.glmnet(x, dat3$Responder,type.measure = "auc", family = "binomial", nfold = 10,  parallel = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
}
cv3 <- search[search$cvm == min(search$cvm), ]
md3 <- glmnet(x, dat3$Responder, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
coef(md3)



# Unterteilung in Training und Test?
# inner loop for feature selection im Trainingsset?
# anschlieﬂend glm model auf ganzes set mit den features, um model coefficients zu erhalten 
# independent test auf Testset?

# oder LOOCV mit Feature selection process f¸r jede LOOCV iteration anders --> anschlieﬂend Mittel der Features und daraus ein Model bilden? 
# vorher Gedanken machen ob logarithmieren oder nicht

library(caret)

# Using caret to perform CV
cctrl1 <- trainControl(method="cv", number=10, returnResamp="all",
                       classProbs=TRUE, summaryFunction=twoClassSummary)

x <- model.matrix(as.formula(paste("Responder ~.")),data=dat3)
x <- x[,-1]

set.seed(849)
test_class_cv_model <- train(x, dat3$Responder, method = "glmnet", 
                             trControl = cctrl1,metric = "ROC",tuneGrid = expand.grid(alpha = seq(0,1,by=0.1),
                                                                                      lambda = seq(0.001,0.1,by = 0.001)))

x <- model.matrix(as.formula(paste("Responder ~ LDH + BRAF + Eosinophile +`hsa-mir-514a-3p`")),data=dat3)
x <- x[,-1]

set.seed(849)
test_class_cv_model <- train(x, dat3$Responder, method = "xgbTree", 
                             trControl = cctrl1,metric = "ROC")



logistic <- glm(Responder ~ LDH + BRAF + Eosinophile +`hsa-mir-514a-3p`,data=dat3, family="binomial")
summary(logistic)
library(pscl)

mcFR2 <- pR2(logistic)

# log likelihood of fit
llfit <- mcFR2[1]
#log likelihood of null model
llnull <- mcFR2[2]

# calculating chi≤
#2*(llfit - llnull)

# calculating p
#chi≤ = 26 with df = 4

#--> p < 0.0001

nullmodel <- glm(Responder ~ 1, data = dat3, family = "binomial")
anova(nullmodel,logistic, test = 'LRT')




