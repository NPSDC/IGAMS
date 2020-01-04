library(randomForest)
library(foreach)
library(doParallel)

#' Gets the cross valiation stages for Random Forest for a given data
#' 
#' @param tr.data normalised training data containing the fea set
#' @param tr.model training model
#' @param folds no of folds for k-fold CV
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return The predicted stages after 10 fold cross validation
cv.rf <- function(tr.data, tr.model, folds, stages.train, cores = 1, sampsize = if (replace) nrow(data) else ceiling(.632*nrow(data)), n.trees = NULL, cla.pr = NULL, cla.rem = NULL)
{
  #writeLines(c(""), "log.txt")
  #sink('log.txt', append = T)
  total.samp <- nrow(tr.data)
  gr <- build.groups(stages = stages.train, num.group = folds, strat = T)
  
  pr.c.list <- vector(mode = 'list', length = length(n.trees))
  if(!is.null(tr.model))
  {
    n.trees <- tr.model$ntree
    cla.pr <- tr.model[['cla.pr']]
    cla.rem <- tr.model[['cla.rem']]
  }
  for(j in seq_along(n.trees))
  {
    
    registerDoParallel(cores = cores)
    predicted.prob <- foreach(i = 1:folds, .combine = c) %dopar%
    {
      train.index = sort(unlist(gr[-i]))
      test.index = sort(unlist(gr[[i]]))
      rf.model <- randomForest(x = tr.data[train.index, ], y = stages.train[train.index], ntree = n.trees[j])
      return(get.pr(tr.model = rf.model, test.data = tr.data[test.index,]))
    }
    registerDoSEQ()
    predicted.prob[unlist(gr)] <- predicted.prob
    pr.c.list[[j]] <- list(pr.prob = predicted.prob, gr = gr, cla.pr = cla.pr,
                           cla.rem = cla.rem)
  }
  gc()
  if(!is.null(tr.model))
    return(pr.c.list[[1]])
  
  res.list <- get.eval.list(stages.train, predict.list = pr.c.list, classifier = 'randomForest', thr = 0.5, cv = T)
  tree.ind <- get.best.mcc.pr_aucs(res.list)
  pr.c.list[[tree.ind]][['ntree']] <- n.trees[tree.ind]
  return(pr.c.list[[tree.ind]])
}

#' Gets the cross valiation models for RFs for the AFs of a feature selection algorithm
#' 
#' @param tr.data normalised training data containing samples as rows and columns as fea
#' @param train.model.list List of trained models w.r.t feature selection algorithm
#' @param folds no of folds for k-fold CV
#' @param fea.list List of AFs of a particular feature selection algorithm
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return List of predicted stages on CV SVM models for AFs
cv.rf.list <- function(tr.data, train.model.list, folds, fea.list, stages.train, cores = 1,
                       sampsize = if (replace) nrow(data)
                       else ceiling(.632*nrow(data)))
{
  rf.pred <- mclapply(seq_along(fea.list), function(i)
  {
    if(is.null(fea.list[[i]]) | (length(fea.list[[i]]) == 1))
      return(NULL)
    else
    return(cv.rf(tr.data = tr.data[, fea.list[[i]]], tr.model = train.model.list[[i]],
                 folds =  folds, 
                 stages.train = stages.train, cores = cores))
  }, mc.cores = cores)
  names(rf.pred) <- names(fea.list)
  return(rf.pred)
}
