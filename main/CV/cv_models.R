library(pamr)
library(randomForest)

library(pROC)
library(foreach)
library(doParallel)
source('main/CV/shrunken.R')
source('main/CV/svm.R')
source('main/CV/rf.R')
source('main/CV/nb.R')
source('main/CV/knn.R')

#' Returns the different cross validation models for all the features for the different classifiers
#' 
#' @param tr.data normalised training data containing samples as rows and columns as features
#' @param fea.list List of features yielded by different feature selection algorithms
#' @param folds Number of folds for k-fold cross validation
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return The list containing the different cross validation predictions for the different features for different classifiers
get.cv.model <- function(tr.data, fea.list, folds = 10, stages.train, tr.model, cores = 1)
{
  registerDoParallel(cores = cores)
  cv.model <- foreach(i = 1:length(fea.list)) %dopar%
  {
    model <- list()
    fea.name <- names(fea.list)[i]
    print(fea.name)
    model[['shrunken']] <- cv.shrunken(tr.data = tr.data, folds = folds, fea.list = fea.list[[fea.name]],
                                     train.models.list = tr.model[[fea.name]]$shrunken,
                                     stages.train = stages.train, cores = cores)
    print('shrunken completed')
    # writeLines(c(""), "log.txt")
    # sink('log.txt', append = T)
    # cat(nrow(tr.data))
    model[['rf']] <- cv.rf.list(tr.data = tr.data, folds = folds, fea.list = fea.list[[fea.name]],
                                train.model.list = tr.model[[fea.name]]$rf, stages.train = stages.train, cores = cores)
    print('RF completed')
    model[['nb']] <- cv.nb.list(tr.data = tr.data, folds = folds, fea.list = fea.list[[fea.name]],
                                train.model.list = tr.model[[fea.name]]$nb, stages.train = stages.train, cores = cores)
    print('NB completed')
    model[['svm']] <- cv.svm.list(tr.data = tr.data, folds = folds, fea.list = fea.list[[fea.name]],
                                  train.model.list = tr.model[[fea.name]]$svm, stages.train = stages.train, cores = cores)
    print('SVM completed')
    model[['knn']] <- cv.knn.list(tr.data = tr.data, folds = folds, fea.list = fea.list[[fea.name]],
                                                 stages.train = stages.train, cores = cores)
    return(model)
  }
  
  registerDoSEQ()
  names(cv.model) = names(fea.list)
  return(cv.model)
}
