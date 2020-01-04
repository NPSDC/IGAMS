library(class)
library(foreach)
library(doParallel)
library(pROC)
library(mccr)
source('main/get_results.R')
#' Gets the cross valiation stages for KNN for given k
#' 
#' @param tr.data normalised training data containing the fea set
#' @param folds no of folds for k-fold CV
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return The predicted stages after k fold CV
cv.knn <- function(tr.data, folds, stages.train, k, cores = 1)
{
  total.samp <- nrow(tr.data)
  gr <- build.groups(stages = stages.train, num.group = folds, strat = T)
  registerDoParallel(cores = cores)
  predicted.prob <- foreach(i = 1:folds, .combine = c) %dopar%
  {
    train.index = sort(unlist(gr[-i]))
    test.index = sort(unlist(gr[i]))
    return(get.pr(tr.data = tr.data[train.index,], test.data = tr.data[test.index, ], 
                 k = k, stages.train = stages.train[train.index]))
  }
  predicted.prob[unlist(gr)] <- predicted.prob
  registerDoSEQ()
  return(list(pr.prob = predicted.prob, gr = gr))
}

#' Gets the best K and given predicted stages
#' 
#' @param tr.data normalised training data containing the fea set
#' @param folds no of folds for k-fold CV
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return The predicted stages and best k for k fold CV
find.best.k <- function(tr.data, folds, stages.train, cores = 1)
{
  class.pr <- levels(as.factor(stages.train))[1]
  class.rem <- levels(as.factor(stages.train))[2]
  cv.pred.prob <- mclapply(seq(folds), function(k){
    cv.knn(tr.data, folds, stages.train, k, cores = cores)
  }, mc.cores = cores)
  mccs <- sapply(seq(folds), function(k)
  {
    mccr(get.order(stages.train, class.pr), 
         get.order(get.class.pred(cv.pred.prob[[k]]$pr.prob, 
                                  class.pr, class.rem), class.pr))
  })
  k <- which.max(mccs)
  res <- list()
  res[['best_k']] <- k
  res[['pred']] <- cv.pred.prob[[k]]
  res[['pred']][['cla.pr']] <- levels(as.factor(stages.train))[1]
  res[['pred']][['cla.rem']] <- levels(as.factor(stages.train))[2]
  return(res)
}

#' Gets the predicted stages and best k for KNN for the AFs of a feature selection algorithm
#' 
#' @param tr.data normalised training data containing samples as rows and columns as fea
#' @param folds no of folds for k-fold CV
#' @param fea.list List of AFs of a particular feature selection algorithm
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return List of predicted stages on CV SVM models for AFs
cv.knn.list <- function(tr.data, folds, fea.list, stages.train, cores = 1)
{
  knn.pred <- mclapply(fea.list, function(fea)
  {
    if(is.null(fea) | (length(fea) == 1))
      return(NULL)
    else
      return(find.best.k(tr.data[, fea], folds, stages.train, cores = cores))
  }, mc.cores = cores)
  names(knn.pred) <- names(fea.list)
  return(knn.pred)
}