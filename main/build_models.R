library(pamr)
library(randomForest)
library(e1071)
library(foreach)
library(doParallel)

source('main/get_results.R')
source('main/get_eval_mat.R')
source('main/CV/cv_models.R')

#' Gets the trained Shrunken Centroid classifiers for the AFs of a feature selection algorithm
#' 
#' @param tr.data normalised training data containing samples as rows and columns as fea
#' @param fea.list List of AFs of a particular feature selection algorithm
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return The list containing the trained models for the AFs
build.shrunken.classifier <- function(tr.data, fea.list, stages.train, cores = 1)
{
  cla.pr <- levels(as.factor(stages.train))[1]
  cla.rem <- levels(as.factor(stages.train))[2]
  shrunken.train.list <- mclapply(fea.list, function(fea)
  {
    if(is.null(fea) | length(fea) == 1)
      return(NULL)
    mod <- pamr.train(data = list(x = as.matrix(t(tr.data[,fea])), y = stages.train))
    mod$cla.pr <- cla.pr
    mod$cla.rem <- cla.rem
    mod
  }, mc.cores = cores)
  names(shrunken.train.list) <- names(fea.list)
  return(shrunken.train.list)
}

#' Returns the trained Random Forest classifiers for AFs of a feature selection algorithm
#' 
#' @param tr.data normalised training data containing samples as rows and columns as fea
#' @param fea.list List of AFs of a feature selection algorithm
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return The list containing the trained models for the AFs
build.rf.classifier <- function(tr.data, fea.list, stages.train, n.trees, cores = 1)
{
  cla.pr <- levels(as.factor(stages.train))[1]
  cla.rem <- levels(as.factor(stages.train))[2]
  rf.train.list <- mclapply(fea.list, function(fea)
  {
    if(is.null(fea) | length(fea) == 1)
      return(NULL)
    
    rf.cv <- cv.rf(tr.data = tr.data[, fea], tr.model = NULL, folds = 5, stages.train = stages.train, cores = cores,
                   n.trees = n.trees, cla.pr = cla.pr, cla.rem = cla.rem )
    
    n.tree <- rf.cv[['ntree']]
    rf.mod <- randomForest(tr.data[,fea], stages.train, ntree = n.tree)
    rf.mod$cla.pr <- cla.pr
    rf.mod$cla.rem <- cla.rem
    rf.mod
  }, mc.cores = 1)
  names(rf.train.list) = names(fea.list)
  return(rf.train.list)
}

#' Returns the trained SVM classifiers for AFs of a feature selection algorithm
#' 
#' @param tr.data normalised training data containing samples as rows and columns as fea
#' @param fea.list List of AFs of a feature selection algorithm
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return The list containing the trained models for the AFs
build.svm.classifier <- function(tr.data, fea.list, stages.train, cores = 1,
                                 gamma = 0, kernel = 'linear', C_set = 1,
                                 class.weights =if(length(levels(stages.levels)) == 4) 
                                   c('stage i' = 1, 'stage ii' =1, 'stage iii' = 1, 'stage iv' =1)
                                 else c('stage i' = 1, 'stage iv' = 1))
{
  
  svm.train.list <- mclapply(fea.list, function(fea)
  {
    if(is.null(fea) | length(fea) == 1)
      return(NULL)
    
    svm.cv <- cv.svm(tr.data = tr.data[, fea], tr.model = NULL, folds = 5, stages.train = stages.train, 
                     cores = cores, C_set = C_set, cla.pr = levels(as.factor(stages.train))[1], cla.rem = levels(as.factor(stages.train))[2])
    
    cost.req <- svm.cv[['C']]
    svm.cv
    svm.model <- svm(x = tr.data[, fea], y = stages.train, kernel = kernel, gamma = gamma, cost = cost.req,
                     probability = T)
    svm.model$cla.pr <- levels(as.factor(stages.train))[1]
    svm.model$cla.rem <- levels(as.factor(stages.train))[2]
    svm.model
  }, mc.cores = 1)
  names(svm.train.list) <- names(fea.list)
  return(svm.train.list)
}

#' Returns the trained naive bayes classifiers for AFs of a feature selection algorithm
#' 
#' @param tr.data normalised training data containing samples as rows and columns as fea
#' @param fea.list List of AFs of a feature selection algorithm
#' @param cores Number of CPU cores to be used
#' @return The list containing the trained models for the AFs
build.nb.classifier <- function(tr.data, fea.list, stages.train, cores = 1)
{
  cla.pr <- levels(as.factor(stages.train))[1]
  cla.rem <- levels(as.factor(stages.train))[2]
  nb.train.list <- mclapply(seq_along(fea.list), function(i){
  if(is.null(fea.list[[i]]) | length(fea.list[[i]]) == 1)
    return(NULL)
  nb.mod <- naiveBayes(tr.data[, fea.list[[i]]], stages.train)
  nb.mod$cla.pr <- cla.pr
  nb.mod$cla.rem <- cla.rem
  nb.mod
  }, mc.cores = cores)
  names(nb.train.list) <- names(fea.list)
  return(nb.train.list)
}

#' Returns the different training classifiers for the AFs of different feature selection algorithms
#' 
#' @param tr.data normalised training data containing samples as rows and columns as fea
#' @param fea.list List of list of AFs yielded by different feature selection algorithms 
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return The list containing the different trained models for every feature in fea.list
get.train.model <- function(tr.data, fea.list, stages.train, C_set = c(1), n.trees = c(500), cores = 1)
{
  #registerDoParallel(cores = cores)
  train.model <- foreach(i = 1:length(fea.list)) %do%
  {
    fea.name <- names(fea.list)[i]
    model <- list()
    print('Starting')
    model[['shrunken']] <- build.shrunken.classifier(tr.data = tr.data, fea.list = fea.list[[fea.name]], 
                                                     stages.train = stages.train, cores = cores)
    print('Completed Shrunken')
    model[['rf']] <- build.rf.classifier(tr.data = tr.data, fea.list = fea.list[[fea.name]], 
                                         stages.train = stages.train, n.trees = n.trees, cores = cores)
    print('Completed RF')
    model[['nb']] <- build.nb.classifier(tr.data = tr.data, fea.list = fea.list[[fea.name]], 
                                         stages.train = stages.train, cores = cores)
    print('Completed NB')
    model[['svm']] <- build.svm.classifier(tr.data = tr.data, fea.list = fea.list[[fea.name]], 
                                           stages.train = stages.train, C_set = C_set, cores = cores)
    print('Completed SVM')
    return(model)
  }
  names(train.model) = names(fea.list)
  return(train.model)
}
