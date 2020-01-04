library(e1071)
library(foreach)
library(doParallel)
source('main/helper_func.R')
source('main/predict_models.R')
#' Gets the cross valiation stages for SVM for a given data
#' 
#' @param tr.data normalised training data containing the gene set
#' @param tr.model training model
#' @param folds no of folds for k-fold CV
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return The predicted stages after 10 fold cross validation
cv.svm <- function(tr.data, tr.model = NULL, folds, stages.train, cores = 1, gamma = 0, kernel = 'linear', C_set = NULL, cla.pr = NULL, cla.rem = NULL,
                   class.weights = sapply(levels(as.factor(stages.train)), function(x){
                     c(x=1)
                   } ))
{
  #sink("log.txt", append=TRUE)
  names(class.weights) = levels(as.factor(stages.train))
  total.samp <- nrow(tr.data)
  gr <- build.groups(stages = stages.train, num.group = folds, strat = T)
  
  if(length(gr) != folds)
    return('not correct')
  registerDoParallel(cores = cores)
  pr.c.list <- vector(mode = 'list', length = length(C_set))
  if(!is.null(tr.model))
  {
    C_set <- tr.model$cost
    cla.pr <- tr.model[['cla.pr']]
    cla.rem <- tr.model[['cla.rem']]
  }
    
  for(j in seq_along(C_set))
  {
    predicted.prob <- foreach(i = 1:folds, .combine = c) %dopar%
    {
      train.index = sort(unlist(gr[-i]))
      test.index = sort(unlist(gr[[i]]))
      svm.model <- svm(x = tr.data[train.index, ], y = stages.train[train.index], kernel = kernel, cost = C_set[j])
      return(get.pr(tr.model = svm.model, test.data = tr.data[test.index,]))
    }
    predicted.prob[unlist(gr)] <- predicted.prob
    pr.c.list[[j]] <- list(pr.prob = predicted.prob, gr = gr, cla.pr = cla.pr,
                           cla.rem = cla.rem)
  }
  
  registerDoSEQ()
  gc()
  if(!is.null(tr.model))
    return(pr.c.list[[1]])
  
  res.list <- get.eval.list(stages.train, predict.list = pr.c.list, classifier = 'svm', thr = 0, cv = T)
  c.ind <- get.best.mcc.pr_aucs(res.list)
  pr.c.list[[c.ind]][['C']] <- C_set[c.ind]
  return(pr.c.list[[c.ind]])
}

#' Gets the predicted stages for SVM for the AFs of a feature selection algorithm
#' 
#' @param tr.data normalised training data containing samples as rows and columns as genes
#' @param train.model.list List of training models
#' @param folds no of folds for k-fold CV
#' @param fea.list List of AFs of a particular feature selection algorithm
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return List of predicted stages on CV SVM models for AFs
cv.svm.list <- function(tr.data, train.model.list, folds, fea.list, stages.train,
                        gamma = 0, kernel = 'linear', cores = 1, cost =1,
                        class.weights =if(length(levels(stages.train)) == 4) 
                          c('stage i' = 1, 'stage ii' =1, 'stage iii' = 1, 'stage iv' =1)
                        else c('stage i' = 1, 'stage iv' = 1))
{
  #writeLines(c(""), "log.txt")
  
  svm.pred <- mclapply(seq_along(fea.list), function(i)
  {
    if(is.null(fea.list[[i]]) | (length(fea.list[[i]]) == 1))
      return(NULL)
    return(cv.svm(tr.data = tr.data[,fea.list[[i]]], folds = folds,
                     stages.train = stages.train, tr.model = train.model.list[[i]], cores = cores))
    
     
  }, mc.cores = cores)
  names(svm.pred) <- names(fea.list)
  return(svm.pred)
}