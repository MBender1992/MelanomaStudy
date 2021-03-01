data(ROCR.hiv)
attach(ROCR.hiv)


out <- cvAUC(hiv.svm$predictions, hiv.svm$labels)

#Plot fold AUCs

plot(out$perf, col="grey82", lty=3, main="10-fold CV AUC")

#Plot CV AUC
plot(out$perf, col="red", avg="vertical", add=TRUE)



cvAUC(model_list$svmRadial$pred$nein, model_list$svmRadial$pred$obs)


predictions <- split(model_list$svmRadial$pred$nein, f = model_list$svmRadial$pred$Resample)
labels <- split(model_list$svmRadial$pred$obs, f = model_list$svmRadial$pred$Resample)
names(predictions) <- NULL 

out <- cvAUC(predictions, labels)

plot(out$perf, col="grey82", lty=3, main="10-fold CV AUC")

#Plot CV AUC
plot(out$perf, col="red", avg="vertical", add=TRUE)

unlist.model(models.lasso.complete, "AUC", "test.metrics") 
unlist.model(models.lasso.complete, "pred.ja", "predictions") 



predictions <- lapply(1:10, function(d){
  test <- sapply(models.lasso.complete[[d]], "[", "predictions")
  lapply(1:10, function(x){
    test[[x]]$pred.ja
  })
})

predictions <- lapply(rapply(predictions, enquote, how="unlist"), eval)

labels <- lapply(1:10, function(d){
  test <- sapply(models.lasso.complete[[d]], "[", "predictions")
  lapply(1:10, function(x){
    test[[x]]$obs
  })
})

labels <- lapply(rapply(labels, enquote, how="unlist"), eval)

out <- cvAUC(predictions, labels)

unlist.model(models.lasso.complete, "AUC", "test.metrics")

data.frame(cvAUC = out$fold.AUC, modelAUC = unlist.model(models.lasso.complete, "AUC", "test.metrics"))





models.lasso.complete$Rep10$Fold7$predictions

pred <- ROCR::prediction(predictions[[98]],labels[[98]])
perf <- ROCR::performance(pred, "tpr", "fpr")

as.numeric(ROCR::performance(pred, measure = "auc", 
                             x.measure = "cutoff")@y.values)

auc(roc(labels[[98]],predictions[[98]], direction = ">", levels = c("nein", "ja")))


test <- ROCR::performance(pred, measure = "auc", 
                          x.measure = "cutoff")

