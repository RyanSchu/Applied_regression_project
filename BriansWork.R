# The reason these are written as functions is that we would like to use this for picking the best model for any Gene down the line :). Each gene might have different assumptions, interactions,
# and even different distributions.


remove_singularities <- function(dataset, gene_name){
  dataset_copy <- dataset
  item <- paste(gene_name, "~.", sep="")
  full_formula <- as.formula(item)
  fit <- lm(full_formula, data=dataset)
  singularities <- attributes(alias(fit)$Complete)$dimnames[[1]]
  for (singularity in singularities){
    dataset_copy[singularity] <- NULL
  }
  return(dataset_copy)
}

read_in_pruned_datasets_for_gene_0.8 <- function(gene_name, path){
  full_path0.8 <- paste(path, gene_name, "_for_r_0.8.txt", sep="")
  Data0.8 <- read.table(full_path0.8, header=TRUE, sep=',')
  Data0.8 <- remove_singularities(Data0.8, gene_name)
  return(Data0.8)
}

gene_data <- read_in_pruned_datasets_for_gene_0.8("ENSG00000068097", "D:\\Project\\ALL\\Selected\\chr17\\")


fit_naieve_model_for_gene_dataset <- function(gene_dataset, gene_name){
  second_half <- names(gene_dataset[, names(gene_dataset) != gene_name])
  formula_as_text <- paste(gene_name, "~.", sep="")
  full_formula <- as.formula(formula_as_text)
  model_name <- paste(gene_name, "Model", sep='')
  f <-  as.formula(paste(gene_name, "~", paste(second_half, collapse="+")))
  print(f)
  fit <- eval(bquote(lm(.(f), data=gene_dataset)))
  # Penalize number of parameters that do not add new information first! There's a lot of variables so this will be helpful.
  return(fit)
}
library(caret)
library(tidyverse)
library(lmtest)
library(MASS)
library(car)
library(reshape)
library(plotmo)
orchestrator <- function(gene_name, path){
  
  set.seed(123)
  # We already note there seems to be a large issue with multicolinearity as a priroi we know that SNP data tends to be highly correlated.
  # Thus we go ahead and look at some of these values first.
  
  # Exploratory Data Analysis in R
  gene_data <- read_in_pruned_datasets_for_gene_0.8(gene_name, path)
  gene_data <- as.data.frame(scale(gene_data))
  hist(gene_data[[gene_name]], main=paste(gene_name, "Gene Expression Distribution", sep=""), ylab="Frequency", xlab="Gene Expression Value")
  boxplot(gene_data[[gene_name]])
  
  R <- cor(gene_data)
  R[R == 1] <- NA #drop perfect one's 
  R[abs(R) < 0.5] <- NA
  R <- na.omit(melt(R)) # melt! 
  R[order(-abs(R$value)),] # sort
  print(R)
  
  
  # For the purposes of being thorough we will go ahead and do a naieve multiple regression model and cross validate(K=10) using caret.
  # We anticpate this model will fail to generalize well for any of the genes we plan to test.
  
  train.control <- trainControl(method = "cv", number=10)
  second_half <- names(gene_data[, names(gene_data) != gene_name])
  
  # This allows us to dynamically populate the train variable with a formula. We use bquote here to do variable substitution.
  
  f <-  as.formula(paste(gene_name, "~", paste(second_half, collapse="+")))
  model <- train(bquote(.(f)), data=gene_data, method="lm", trControl=train.control)
  print(model)
  
  # We do diagnostics here.
  naieve_model <- fit_naieve_model_for_gene_dataset(gene_data, gene_name)

  # We can see that the vifs for this gene are very large suggesting severe problems with multicolinearity. 
  vifs <- car::vif(naieve_model)
  print(vifs[order(-vifs)])
  
  # We do a Breusch-Pagan test to see if there is constant variance. We also do a plot of the residuals versus the predictors. 
  print(bptest(naieve_model))
  plot(naieve_model, which=1)
  plot(naieve_model)
  
  # We do not anticipate any issue with autocorrelation.
  print(dwtest(naieve_model))
  
  # We do a shaprio wilks normality test and a QQPlot to see if there is an issue with the normality assumption. 
  print(shapiro.test(studres(naieve_model)))
  qqnorm(studres(naieve_model))
  qqline(studres(naieve_model))

  # We note that there is an issue with normality as indicated by a qqplot of the residuals in this model. 
  

  
  # We look at an added variable plot and note that a lot of the predictors do not add new information when other predictors are in the model as indicated by a horizontal
  # line
  # avPlots(naieve_model)
  

  # We will use something called foba from the Caret package that does variable selection(number of variables specified by k)
  # We also use a lambda value(equivalent to K)
  #grid <- expand.grid(k=c(2,3,4,5,6,7,8,9,10, 15,20, 25, 40), lambda=c(0, 0.001, 0.01))
  
  # There doesn't appear to be an issue with unequal variance but rather there isn't enough representation of the fitted values beyond the range from 2 to 3.
  # Ideally, we would want to get a much larger sample size in a follow-up study to remedy this.
  
  # We notice that there is a ton of predictors that don't have an association with the response when all other variables are accounted for in the model.
  # Doing backward selection here would be challenging as a lot of the predictors are strongly correlated. We will need to consider using ridge regression as 
  # as result but would also like to have the model be as simple as possible. We probably will use ElasticNet to combine LASSO and RIDGE penalty terms.
  
  
  # This was extremely challenging to figure out. Multicolinearity is a nasty problem.
  
  
  # Basically, the model was overfit to the data. On the training data R-squared average was 0.80.
  # However on the validation dataset(80/20 split) the R-squared average was 0.2 suggesting this model probably overfits the data very easily even after backward
  # selection. 
  
  # Upon further evaluation it became clear that the model itself has extreme issues with multicolinearity that led to high variance in the Beta coefficients.
  # We still would like to eliminate some variables however to keep the model simpler so that it can generalize well and be added to the broader population.
  # Thus we wanted a mixture of sparse feature vectors and weight shrinkage to deal with multicolinearity and feature selection without having to do a separate 
  # step. Thus, we chose ElasticNet to help the model generalize well. For being thorough, we used Caret to do 10-fold cross validation using both pure Ridge and 
  # ElasticNet with parameter grids for both. #glmnet automatically scales the variables. 
  # Noba is greedy variable selection and Ridge! This is very useful for these types of issues.
  
  second_half <- names(gene_data[, names(gene_data) != gene_name])
  f <-  as.formula(paste("1/sqrt(",gene_name, "+2)", "~", paste(second_half, collapse="+")))
  print(f)
  grid <- expand.grid(lambda=c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
  grid2 <- expand.grid(alpha=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), lambda=c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3,0.4,0.5,0.6,0.7,0.8))
  model <- train(bquote(.(f)), data=gene_data, method="ridge", trControl=train.control, tuneGrid = grid, standardize=FALSE)
  model_glmnet <- train(bquote(.(f)), data=gene_data, method="glmnet", trControl=train.control, tuneGrid = grid2, standardize=FALSE)
  final_model_ridge <- model$finalModel
  final_model_glmnet <- model_glmnet$finalModel
  final_model <- model$finalModel
  ridge_mse <- model$results[which.min(model$results[, "RMSE"]), ][, "RMSE"]
  glmnet_mse <- model_glmnet$results[which.min(model$results[, "RMSE"]), ][, "RMSE"]
  if (ridge_mse < glmnet_mse){
   final_model <- final_model_ridge
   print(model)
  }
  else{
   final_model <- final_model_glmnet
    print(model_glmnet)
  }
  grid3 <- expand.grid(k=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,26), lambda=c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
  model <- train(bquote(.(f)), data=gene_data, method="foba", trControl=train.control, tuneGrid=grid3)
  print(model)
}

# We notice a huge improvment in Cross Validated R-squared. We also note that RMSE is at a minimum with alpha = 0.7 and lambda = 0.05 considering the grid above.
# We would also still like to notice if the assumptions of the model are met with these parameters. 

# We will fit the final model using all of the data using glmnet and look at the assumptions.


data <- read_in_pruned_datasets_for_gene_0.8("ENSG00000142794", "D:\\Project\\ALL\\Selected\\chr1\\")
fit_naieve_model_for_gene_dataset <- function(gene_dataset, gene_name){
  second_half <- names(gene_dataset[, names(gene_dataset) != gene_name])
  formula_as_text <- paste(gene_name, "~.", sep="")
  full_formula <- as.formula(formula_as_text)
  model_name <- paste(gene_name, "Model", sep='')
  f <-  as.formula(paste(gene_name, "~", paste(second_half, collapse="+")))
  print(f)
  fit <- eval(bquote(lm(.(f), data=gene_dataset)))
  # Penalize number of parameters that do not add new information first! There's a lot of variables so this will be helpful.
  return(fit)
  
}

# Best model performance comes from FOBA(Adaptive Forward Backward Algorithm!:

#  lambda  k   RMSE       Rsquared   MAE   
#
# 0.300   10  0.1407143  0.5241538  0.09770583

# The algorithm is greedy and in R it is implemented and allows for a ridge lamdba coefficient to be specified. 
# This gives us a lot of power to work through multicolinearity and yield sparse solutions. ElasticNet might also make
# sense but it appears this performs slightly better empirically.



