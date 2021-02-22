setwd("Z:/Aktuell/Eigene Dateien/Eigene Dateien_Marc/R/base_scripts")
#loading custom functions
source("R_functions_Melanoma.R")
# TeST

library(caret)
library(AppliedPredictiveModeling)
library(car)
library(pROC)
library(doParallel) 
library(gbm)
library(tidyverse)
library(corrplot)
library(PerformanceAnalytics)
library(caretEnsemble)
library(kernlab)
library(xgboost)
library(RANN) 
library(mice)
library(missForest)
library(cvAUC)
library(devtools)
library(cutpointr)

# functions used 
#....................................................................................................................
# this function calculates accuracy, F1, spec, sens,... for the Folds    #
# within cross validation using the resamples with the optimal tuning    #
# parameters

cv_res <- function(modelNames){
  
  # show ROC of models
  model_results <- lapply(names(modelNames), function(x){
    modelNames[[x]]$results %>%
      filter(ROC == max(ROC)) %>%
      select(c(ROC, Sens, Spec)) %>% .[1,]
  }) 
  
  model_results <- do.call(rbind.data.frame, model_results) 
  
  
  ## extract those resample with the optimal tuning parameters from the model list for each model
  model_best_tune <- lapply(names(modelNames), function(x){
    tmp <- modelNames[[x]]$pred
    
    for (i in 1:length(modelNames[[x]]$bestTune)){
      tmp   <- tmp %>%
        filter(!!sym(names(modelNames[[x]]$bestTune)[i]) == unlist(modelNames[[x]]$bestTune[i]))
    }
    
    return(tmp)
  })
  
  
  # calculate accuracy, sensitivity and so on for the cross validated models 
  cv_results  <- sapply(1:length(modelNames), function(y){
    
    resamples <- modelNames[[1]]$pred$Resample %>% unique()
    
    sapply(resamples, function(x){
      tmp_res <- model_best_tune[[y]] %>%
        filter(Resample == x)
      
      conf_mat <- confusionMatrix(tmp_res$pred, tmp_res$obs)
      conf_mat$byClass
    }) %>%
      t() %>%
      data.frame() %>%
      summarize_all(mean, na.rm =T)
  })%>%
    data.frame() %>%
    setNames(names(modelNames)) %>%
    t()
  
  cbind(ROC = model_results$ROC,cv_results)
  
}












#....................................................................................................................
# this function averages results from cv_res function for different Seeds    #
modelSeed_meanSD <- function(y){
  model_metrics <- lapply(1:length(y), function(x){
    cv_res(y[[x]])
  })
  
  metricList <- lapply(names(y[[1]]),function(modelName){
    tmp <- do.call(rbind, model_metrics) %>%
      data.frame() %>%
      rownames_to_column("model")%>%
      filter(str_detect(model, modelName)) %>%
      select(-model) %>%
      mutate_all(as.numeric) 
    
    tmp_mean <- tmp %>%
      summarize_all(mean)
    tmp_sd <- tmp %>%
      summarize_all(sd)
    
    
    data.frame(mean = t(tmp_mean), sd = t(tmp_sd)) %>%
      rownames_to_column("metric") %>%
      mutate(error = qt(0.975, length(y) -1)*sd/sqrt(length(y)),
             lower = mean - error,
             upper = mean + error) %>%
      column_to_rownames("metric")
  })
  
  names(metricList) <- names(y[[1]])
  return(metricList)
}





#....................................................................................................................
# calculate model metrics 
calcTestMetrics <- function(data, models){
  metricList_test <- lapply(models,function(x){
    tmp <- do.call(rbind, data) %>%
      data.frame() %>%
      rownames_to_column("model")%>%
      filter(str_detect(model, x)) %>%
      select(-model) %>%
      mutate_all(as.numeric) 
    
    tmp_mean <- tmp %>%
      summarize_all(mean)
    tmp_sd <- tmp %>%
      summarize_all(sd)
    
    n <- length(names(model_list[[1]]))
    
    data.frame(mean = t(tmp_mean), sd = t(tmp_sd)) %>%
      rownames_to_column("metric") %>%
      mutate(error = qt(0.975,n-1)*sd/sqrt(n),
             lower = mean - error,
             upper = mean + error) %>%
      column_to_rownames("metric")
  })
  
  names(metricList_test) <- names(model_list[[1]])
  return(metricList_test)
}





#....................................................................................................................
# function to combine predicted values with real values from a list
# x: name of list element
# model: name of model
predictList <- function(x,model){
  data.frame(Sample = x,
             Responder = predict.train(model_list[[x]][[model]], newdata = dat_split[[x]]$testing,type="prob"),
             RealClass = dat_split[[x]]$testing$Responder) 
}



#....................................................................................................................
# function to calculate TPR and FPR for different cutoffs averaged over a distinct number of samples defined by "supgroup"
# x: predicted values
# class: class containing real values
mean_roc <- function(data, cutoffs = seq(from = 0, to = 1, by = 0.1)) {
  map_df(cutoffs, function(cp) {
    out <- cutpointr(data = data, x = PredictionValues, class = RealClass,
                     subgroup = Sample, method = oc_manual, cutpoint = cp,
                     pos_class = "neg", direction = ">=")
    data.frame(cutoff = cp, 
               sensitivity = mean(out$sensitivity),
               specificity = mean(out$specificity))
  })
}






#....................................................................................................................
# function to calculate AUC using predicted and observed values
# probs: predicted values (probabilities)
# true_Y: observed values
getROC_AUC <- function(probs, true_Y){
  probsSort <- sort(probs, decreasing = TRUE, index.return = TRUE)
  val <- unlist(probsSort$x)
  idx <- unlist(probsSort$ix)  
  
  roc_y <- true_Y[idx];
  stack_x <- cumsum(roc_y == "pos")/sum(roc_y == "pos")
  stack_y <- cumsum(roc_y == "neg")/sum(roc_y == "neg")    
  
  auc <- sum((stack_x[2:length(roc_y)]-stack_x[1:length(roc_y)-1])*stack_y[2:length(roc_y)])
  return(list(stack_x=stack_x, stack_y=stack_y, auc=auc))
}








#....................................................................................................................
# function to calculate the average AUC of multiple ROC models
# samples: number of different models
mean_AUC <- function(data,samples){
  tmp <- sapply(samples, function(x){
    dat <- data %>% 
      filter(Sample == x)
    getROC_AUC(dat$PredictionValues, dat$RealClass) %>%
      .$auc
  }) 
  data.frame(mean.auc = mean(tmp),
             error = confInt(tmp)) %>%
    mutate(lower = mean.auc- error, 
           upper = mean.auc+error)
  
}









#....................................................................................................................
# function to combine predictions from different models (generated by different seeds) within a list
# out object gives the predictions table, ROC and AUC 
pred_summary <-function(model, list.element = 1:length(n)){
  pred <- map_df(list.element,model = model, predictList) %>% 
    select(-Responder.pos) %>%
    setNames(c("Sample","PredictionValues", "RealClass"))
  mr <- mean_roc(pred)
  ma <- mean_AUC(pred, list.element)
  list(predictions = pred,
       ROC = mr,
       AUC = ma
  )
}



#....................................................................................................................
# load saved models
loadModels <- function(model){
  if(model == "M1"){
    readRDS(file = "model_base_seeds.Rds")
  } else if(model == "M2") {
    readRDS(file = "model_miR_seeds.Rds")
  } else if(model == "M3") {
    readRDS(file = "model_best_perf_seeds.Rds")
  } else {
    stop("Please specify model (one of M1 (base model), M2 (miR model), M3 (best performing model))")
  }
}

# set working directory
setwd("Z:/Aktuell/Eigene Dateien/Eigene Dateien_Marc/R/Projekte/Doktorarbeiten_Melanom_Liquid_Biopsies/Daten")

# F?r Ver?ffentlichung Worklfow gut darstellen (Bhalla et al. prediction  and Analysis of Skin Cancer Progression using Genomics Profiles of Patients)
# verschiedene Algorithmen darstellen und vergleichen, verschiedene Feature selection Methoden vergleichen 
# (nur signifikante, RFE, nur signifikante ohne high correlated)

###################################################################################
#                                                                                 #
#                             1. raw Data processing                              #
#                                                                                 #
###################################################################################

# load data
dat_combined <- load_melanoma_data(characterAsFactor = TRUE) %>% filter(!is.na(Responder))


#detect columns containing NAs
NAs <- sapply(dat_combined, function(df){
  sum(is.na(df) ==TRUE)/length(df);
})

#index of columns which contain NAs
which(NAs != 0)

# impute missing values in numeric data (befallene Organe, age,S100) by random forest iteration
dat_imp <- dat_combined %>% select_if(is.numeric) %>% as.data.frame() %>% missForest()

## add NA in factors as factor level and indclude imputed values for NA in numerical data
dat_complete <- dat_combined %>%
  # mutate(Hirnmetastase = addNA(Hirnmetastase),
  #        Stadium = addNA(Stadium),
  #        adjuvant_IFN = addNA(adjuvant_IFN)) %>%
  mutate(age = dat_imp$ximp$age,
         befallen_Organe = dat_imp$ximp$befallen_Organe, Eosinophile = dat_imp$ximp$Eosinophile,
         CRP = dat_imp$ximp$CRP,
         LDH = dat_imp$ximp$LDH,
         S100 = dat_imp$ximp$S100)

#summary statistics
summary(dat_complete)







###################################################################################
#                                                                                 #
#                             2. Data transformation                              #
#                                                                                 #
###################################################################################
#exclude Responder from one-hot encoding
index_Responder <- grep("Responder", colnames(dat_complete))

# one hot encoding for categorical variables
dmy <- dummyVars(" ~ .", data = dat_complete[,-index_Responder])
dat_dummy <- data.frame(predict(dmy, newdata = dat_complete[,-index_Responder]))

#add Responder column 
dat_dummy <- dat_dummy %>% mutate(Responder = dat_complete$Responder)











###################################################################################
#                                                                                 #
#                    3.Feature selection/inspection                               #
#                                                                                 #
###################################################################################

# use only significantly changed parameters in model (+ miR-137 as RFE suggested it to be a good feature, miR-514a 
# lowly expressed, so maybe not?, miR-509-3p dazu? Feature miR-29c-3p zwar nicht signifikant aber geringer p-value und performt gut bei RFE
# allerdings hohe Korrelation mit miR-29b-3p)
# colums containg more than 5% NAs were omitted for the analysis other missing values were imputed by a random forest imputation
# miR-34a und miR-514 lieferten keine gute Performance
dat_signif <- dat_dummy %>% select(X.hsa.mir.132.3p.,X.hsa.mir.137.,X.hsa.mir.197.3p.,X.hsa.mir.214.3p.,X.hsa.mir.514a.3p., BRAFmu, BRAFnras, BRAFwt, Eosinophile, CRP, LDH, S100,Responder)

#plot correlation of miRNAs (scatter plot)
# source_gist("524eade46135f6348140")
# png("miR_corr.png", units="in", width=4.5, height=3, res=1200)
# dat_signif %>% ggplot(aes(X.hsa.mir.29c.3p.,X.hsa.mir.29b.3p.)) +
#   stat_smooth_func(geom="text",method="lm",hjust=0,parse=T)  +
#   geom_smooth(method="lm", se=F,color='skyblue') +
#   geom_point(shape =1) + 
#   xlab("miR-29c-3p expression (a.u.)") +
#   ylab("miR-29b-3p expression (a.u.)") +
#   theme_bw()+
#   theme(legend.key.size = unit(1,"line"),
#         panel.grid.minor=element_blank(),
#         strip.background=element_blank())
# dev.off()

#log transform lab parameters as they seem to follow a log-normal distribution and remove highly correlated miR-29c
dat_model <- dat_signif %>%
  mutate(Eosinophile = ifelse(is.infinite(log(Eosinophile )),0,log(Eosinophile )),
         CRP = log(CRP),
         X.hsa.mir.137. = log(X.hsa.mir.137.)) #%>%
#mutate(Responder = factor(Responder, levels = c("pos", "neg"))) # Je nachdem ob wir Ansprechen oder Nicht-Ansprechen vorhersagen wollen



# die 3 transformierten Variablen zeigen durch Transformation eine deutlich bessere Approximation an die Normalverteilung
# Die Verteilung von miR-29b ist schon untransformatiert sehr nah an einer Normalverteilung 
# Die Transformation von LDH und S100 ?ndert nicht viel an der "Normalit?t" der Verteilung
# -> Model einmal mit log-transformierten und einmal mit untransformierten LDH,S100, miR-29b gefittet
# -> besseres Ergebnis wenn LDH und S100 nicht transformiert werden

# split data into 5 different training and test sets to account for random effects 
n <- c(14,237,1576,2300,4032)
dat_split <- lapply(n, function(x){
  set.seed(x, sample.kind = "Rounding")
  trainIndex <- createDataPartition(dat_model$Responder, p = .7, 
                                    list = FALSE, 
                                    times = 1)
  train <-data.frame(dat_model[ trainIndex,])
  testing <- data.frame(dat_model[-trainIndex,])
  list(train = train,
       testing = testing)
})













###################################################################################
#                                                                                 #
#                                    4. Modelling                                 #
#                                                                                 #
###################################################################################

# activate parallel computing
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)

# preprocess within model as values are back transformed afterwards, z-value scaling used
# to customize the model and change tuning parameters use tuneList and define tuneGrid or tuneLength

# define control for cross validation
set.seed(12,sample.kind = "Rounding")
trainctrl <- trainControl(classProbs= TRUE,
                          verboseIter = TRUE,
                          summaryFunction = twoClassSummary,
                          method = "repeatedcv",
                          number = 10,
                          repeats = 3,
                          returnResamp = "all", # if error: change back to final
                          savePredictions = "all",
                          returnData = TRUE,
                          allowParallel = TRUE)

# define algorithms and specify hyperparamters
algorithmList <- c("xgbTree","svmRadial","lda","glmnet","knn") #svm durch andere Algorithmen ersetzen? lvq zb use more models?
TuneList <- list(rf=caretModelSpec(method="rf", tuneGrid=data.frame(.mtry=1:10)))



# model mit 2 miRNAs, 4 serum Markern und BRAF 
model_miR_seeds <- lapply(1:length(n), function(x){
  set.seed(222,sample.kind = "Rounding")
  tmp <- caretList(
    Responder~., data = dat_split[[x]]$train[,c(1,2,13)],
    trControl = trainctrl,
    methodList = algorithmList,
    tuneList = TuneList,
    continue_on_fail = FALSE,
    preProcess = c("center","scale")
  )
  return(tmp)
})

# load model data
model_list <- loadModels("M3") # M1 = base model, M2 = mir model, M3 = best performing model
dat_split <- readRDS("dat_split.RDS")



# calculate performance scores for cross validated models 
# calculate sd and ci based on the average of the models generated by different seeds
modelMetrics_train <- modelSeed_meanSD(model_list)

# convert list into dataframe
modelMetrics_train_2 <-lapply(names(modelMetrics_train), function(x){
  modelMetrics_train[[x]] %>% rownames_to_column("metric")
})

# rename list elements based on model name 
names(modelMetrics_train_2) <- names(modelMetrics_train)


# plot ROC with confidence intervals averaged over training sets from different seeds 
png("model_resamples.png", units="in", width=3.5, height=2.6, res=1200)
modelMetrics_train_2 %>%
  bind_rows(.id = "model") %>%
  filter(metric == "ROC") %>%
  select(-metric) %>% 
  mutate(model = factor(model, levels=unique(model[order(mean)]), ordered=TRUE)) %>%
  ggplot(aes(model, mean)) +
  geom_point(color = "skyblue", shape = 1, size = 2, stroke = 0.7) + 
  geom_errorbar(aes(ymin = lower, ymax = upper),color = "skyblue",size = 0.7, width = 0.1) +
  coord_flip() +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        axis.title.y = element_blank(),
        panel.border = element_rect(color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x= element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x=unit(0.2,"cm")
  ) +
  ylab("ROC") +
  scale_y_continuous(sec.axis = dup_axis())
dev.off() 


###################################################################################
#                                                                                 #
#   5.            Evaluate models with test data                                  #
#                                                                                 #
###################################################################################


#### Evaluate Algorithm on test set for 5 different seeds 

# calculate ROC and AUC for each model
model_names <- names(model_list[[1]])

pred_res <- lapply(model_names, pred_summary)
names(pred_res) <- model_names

# extract AUC for every model
res_AUC <- map_df(model_names,function(x){
  pred_res[[x]]$AUC
}) 
rownames(res_AUC) <- model_names

# extract ROC for every model
res_ROC <- map_df(model_names,function(x){
  pred_res[[x]]$ROC %>% 
    mutate(model = x)
}) 
rownames(res_AUC) <- model_names

#define optimal cutpoints for eachmodel 
opt_cut <- res_ROC %>% 
  mutate(youden = sensitivity+specificity-1) %>%
  group_by(model) %>%
  filter(youden == max(youden)) 
  
# define x and y coordinates of text
text.pos <- res_AUC %>% 
  rownames_to_column("model") %>%
  mutate(sensitivity = 0.25, specificity = 0.25) %>%
  mutate(mean = paste("AUC:",round(mean.auc, 3))) %>%
  select(c(mean, sensitivity, specificity, model))


# plot ROC curves
png("ROC.png", units="in", width=6, height=5, res=1200)
res_ROC %>% 
  .[-c(42,52,53),] %>% # were plotted as duplicates and wronlgy so for Model 3
  ggplot(aes(x = 1 - specificity, y = sensitivity)) + geom_step() +
  facet_wrap(~model, scales = "free") +
  geom_point(size= 0.8) + geom_text(data = text.pos,label = text.pos$mean, size =3) +
  theme_bw() +   theme(aspect.ratio = 1,
                       axis.ticks.length = unit(0.2,"cm"),
                       strip.background = element_blank(),
                       panel.spacing = unit(1.5, "lines"),
                       strip.text = element_text(face = "bold")) +
  geom_abline(slope = 1, intercept = 0, color = "red", lty = 2) +
  xlab("1-Specificity") +  ylab("Sensitivity") +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,by = 0.2))+
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by = 0.2)) 
dev.off() 


# calculate model metrics for each individual model for each seed
modelMetrics_test <- lapply(1:length(n), function(x){
  sapply(model_names,function(model){
    confusionMatrix(table(predict.train(model_list[[x]][[model]], newdata = dat_split[[x]]$testing, type="prob")[,"neg"] >= 0.5,#opt_cut[opt_cut$model == model,]$cutoff, 
                          dat_split[[x]]$testing$Responder == "neg"), positive = "TRUE")$byClass
  }) %>%
    data.frame() %>% 
    t()
})

# calculate avg model metrics summarized over different seeds for individual models
modelMetrics_test_avg <-  calcTestMetrics(data = modelMetrics_test, models = model_names)







###################################################################################
#                                                                                 #
#   6.  calculate Importance of the final model                                   #
#                                                                                 #
###################################################################################



#extract Importance from each separate model
Importance <- map_df(1:length(n), function(x){
  map_df(model_names, function(model){
    varImp(model_list[[x]][["rf"]])$importance %>%
      select(-contains("pos")) %>%
      setNames("VarImp") %>%
      rownames_to_column("Feature") %>%
      mutate(model = model) %>%
      mutate(Feature = str_replace_all(.$Feature, "X.hsa.","")) %>%
      mutate(Feature = str_replace_all(.$Feature, "\\.$","")) %>%
      mutate(Feature = str_replace_all(.$Feature, "\\.","-"))
    })
  }) %>% 
  group_by(model, Feature) %>%#
  summarize(mean = mean(VarImp))
  
         
#plot results similar to varImp function (not working for stacked models, therefore the workaround)
png("varImp.png", units="in", width=5, height=3.5, res=1200)
Importance %>% filter(model %in% c("rf", "xgbTree")) %>%
ggplot(aes(x=reorder(Feature,mean),y=mean)) +
  facet_wrap(~model) +
  geom_bar(position="dodge",stat="identity",width = 0, color = "black") + 
  coord_flip() + 
  geom_point(color='skyblue') + 
  ylab("variable importance")+
  ylim(0,100) +
  xlab("")+
  theme_bw() +   theme(aspect.ratio = 1,
                      axis.ticks.length = unit(0.2,"cm"),
                      axis.ticks.y = element_blank(),
                      strip.background = element_blank(),
                      panel.spacing = unit(1.5, "lines"),
                      panel.grid.major.x = element_blank(),
                      panel.grid.minor.x= element_blank(),
                      strip.text = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) 
dev.off()   


# turn of parallel computing
stopCluster(cl)




