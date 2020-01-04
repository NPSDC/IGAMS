library(e1071)
library(foreach)
library(doParallel)

#' Gets the cross valiation stages for Naive Bayes for a given data
#' 
#' @param tr.data normalised training data containing the fea set
#' @param tr.model training model
#' @param folds no of folds for k-fold CV
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return The predicted stages after 10 fold cross validation
cv.naiveBayes <- function(tr.data, tr.model, folds, stages.train, cores = 1)
{
  total.samp <- nrow(tr.data)
  gr <- build.groups(stages = stages.train, num.group = folds, strat = T)
  registerDoParallel(cores = cores)
  predicted.prob <- foreach(i = 1:folds, .combine = c) %dopar%
  {
    train.index = sort(unlist(gr[-i]))
    test.index = sort(unlist(gr[i]))
    nb.model <- naiveBayes(x = tr.data[train.index,], y = stages.train[train.index])
    return(get.pr(tr.model = nb.model, test.data = tr.data[test.index,]))
  }
  predicted.prob[unlist(gr)] <- predicted.prob
  registerDoSEQ()
  return(list(pr.prob = predicted.prob, gr = gr, cla.pr = tr.model[['cla.pr']],
              cla.rem = tr.model[['cla.rem']]))
}

#' Gets the predicted stages for Naive Bayes for the AFs of a feature selection algorithm
#' 
#' @param tr.data normalised training data containing samples as rows and columns as fea
#' @param train.model.list List of training models
#' @param folds no of folds for k-fold CV
#' @param fea.list List of AFs of a particular feature selection algorithm
#' @param stages.train Stage of every sample in training data
#' @return List of predicted stages on CV Naive Bayes models for AFs
cv.nb.list <- function(tr.data, train.model.list, folds, fea.list, stages.train, cores = 1)
{
  nb.pred <- mclapply(seq_along(fea.list), function(i)
  {
    if(is.null(fea.list[[i]]) | (length(fea.list[[i]]) == 1))
      return(NULL)
    else
      return(cv.naiveBayes(tr.data = tr.data[, fea.list[[i]]], 
                           tr.model = train.model.list[[i]], folds, stages.train))
  }, mc.cores = cores)
  names(nb.pred) <- names(fea.list)
  return(nb.pred)
}
