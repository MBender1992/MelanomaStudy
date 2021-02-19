# create bootstrap sample
set.seed(4)


resample <- dat_fct[sample(nrow(dat_fct), replace = TRUE), ]

# define test and training set
set.seed(123)
ind.train <- createDataPartition(resample$Responder, p = 0.7, list = FALSE)

train.data  <- resample[ind.train, ] 
test.data <- resample[-ind.train, ] 

# transform training data
tmp.train <- train.data %>% select(where(is.numeric)) 
fctrs <- train.data %>% select(!where(is.numeric))
train <- data.frame(cbind(log(tmp.train  + 1), fctrs))

# transform testing data
tmp.test  <- test.data %>% select(where(is.numeric))
fctrs <- test.data %>% select(!where(is.numeric))
test <- data.frame(cbind(log(tmp.test  + 1), fctrs))

# dummy encoding of factors with model.matrix
x.train <- model.matrix(Responder~.,data=train)[,-1] 
y.train <- train$Responder

# dummy encoding of factors with model.matrix for test set
x.test <- model.matrix(Responder~.,data=test)[,-1] 
y.test <- test$Responder





