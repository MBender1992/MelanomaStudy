predictions <- lapply(1:10, function(d){
  test <- sapply(models.lasso.relaxedLasso[[d]], "[", "predictions")
  lapply(1:10, function(x){
    1-test[[x]]$pred.ja
  })
})

predictions <- lapply(rapply(predictions, enquote, how="unlist"), eval)

labels <- lapply(1:10, function(d){
  test <- sapply(models.lasso.relaxedLasso[[d]], "[", "predictions")
  lapply(1:10, function(x){
    test[[x]]$obs
  })
})

labels <- lapply(rapply(labels, enquote, how="unlist"), eval)

out <- cvAUC(predictions, labels)



ci.cvAUC(predictions, labels)

plot(out$perf, col="grey82", lty=3, main="10-fold CV AUC")

#Plot CV AUC
plot(out$perf, col="red", avg="vertical", add=TRUE)
abline(0,1)



# to do
# write function to implement ci calculation with cvAUC
# show cvAUC curve? 
# ADD ci to curve? 