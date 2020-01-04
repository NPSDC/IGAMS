library(pamr)
library(pROC)
library(mccr)
library(foreach)
library(doParallel)
source('main/helper_func.R')
source('main/get_results.R')
#' Finds the threshold having the maximum AUC
#' 
#' @param cv.models List of cross validation models
#' @param stages of training data
#' @return Index of best threshold
find.best.shrunken.threshold <- function(cv.models, stages)
{
  #print(cv.models)
  cla.pr <- levels(as.factor(stages))[1]
  thr.inds <- sapply(cv.models, function(model)
  {
    if(is.null(model))
      return(NULL)
    else
    {
      pamr.mcc.comb <- sapply(seq_along(model$threshold), function(x)
      {
        mccr(get.order(stages, cla.pr), get.order(model$yhat[,x], cla.pr))
      })
      print(max(pamr.mcc.comb))
      #print(pamr.aucs.comb)
      thr.ind = sort(which(pamr.mcc.comb == max(pamr.mcc.comb)), decreasing = T)[1]
    }
  })
  
  names(thr.inds) <- names(cv.models)
  return(thr.inds) 
}

#' Gets the cross valiation models for the shrunken centroids for the AFs of a feature selection algorithm
#' 
#' @param tr.data normalised training data containing samples as rows and columns as fea
#' @param folds no of folds for k-fold CV
#' @param fea.list List of different groups of a particular feature selection algorithm
#' @param train.models.list List of training models of shrunken for the different fea sets
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return List containing the predicted stages, cv models and threshold indexes for AFs
cv.shrunken <- function(tr.data, folds, fea.list, train.models.list, stages.train, cores = 1)
{
  cv.models <- mclapply(seq_along(fea.list), function(i)
  {
     if(is.null(fea.list[[i]]) | (length(fea.list[[i]]) == 1))
       return(NULL)
     else
      return(pamr.cv(train.models.list[[i]], list(x = t(as.matrix(tr.data[, fea.list[[i]]])),
                                         y = stages.train), folds))
  }, mc.cores = cores)
  names(cv.models) <- names(fea.list)
  thr <- unlist(find.best.shrunken.threshold(cv.models, stages.train))
  pred <- lapply(seq_along(cv.models), function(i)
  {
    if(is.null(cv.models[[i]]))
      return(NULL)
    return(cv.models[[i]]$prob[,1,thr[i]])
  })
  names(pred) <- names(fea.list)
  
  res <- list()
  res[['cv.models']] <- cv.models
  res[['thr']] <- thr
  res[['pred']] <- lapply(seq_along(pred), function(i)
  {
    list(pr=pred[[i]],cla.pr=train.models.list[[i]][['cla.pr']],
         cla.rem=train.models.list[[i]][['cla.rem']])
  })
  names(res[['pred']]) <- names(train.models.list)
  return(res)
}