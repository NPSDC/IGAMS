library(doParallel)
library(foreach)
library(pROC)
library(PRROC)
library(mccr)
library(caret)

source('main/helper_func.R')
#' Computes the error for a given confusion matrix
#' 
#' @param conf.mat Table/Confusion matrix of the predicted vs actual stages
#' @return Returns the vector of errors of length equal to the number of rows in conf.mat
compute.error.conf.mat <- function(conf.mat)
{
  total.samps <- apply(conf.mat, 1, sum)
  error <- 1-diag(conf.mat)/total.samps
  return(error)
}

#' Binds the error with the given confusion matrix
#' 
#' @param conf.mat Table/Confusion matrix of the predicted vs actual stages
#' @return Returns the table with error bound
create.mat.error <- function(conf.mat)
{
  return(cbind(conf.mat, compute.error.conf.mat(conf.mat)))
}

#' Gets the evaluation parameters given the predicted and actual stages
#' 
#' @param actual.stages Actual labels/stages for the data
#' @param pred.stages.score List containing predicted prob(and gr for CV) with cla.pr and cla.rem
#' @param classifier classifier algorithm used
#' @param thr Threshold for assigning stages
#' @param cv whether evaluating for cross validation or test
#' @return List containing the evaluated parameters and confusion matrix
get.eval <- function(actual.stages, pred.stages.score, classifier,thr = 0.5, cv = F)
{
    res <- list()
    cla.pr <- pred.stages.score$cla.pr
    cla.rem <- pred.stages.score$cla.rem
    pred.stages <- as.factor(get.class.pred(pr = pred.stages.score$pr, 
                                            class.pr = cla.pr, 
                                            class.rem = cla.rem, 
                                            thr = thr))
    if(is.null(pred.stages.score$pr))
      return(NULL)
  # else
  #   pred.stages <- as.factor(get.class.pred(pr = pred.stages.score, 
  #                               class.pr = cla.pr, 
  #                               class.rem = cla.rem, 
  #                               thr = thr))
  
  res[['conf_mat']] <- create.mat.error(table(actual.stages, pred.stages))
  res[['eval']] <- list()
  res[['eval']][['mcc']] <- mccr(get.order(actual.stages, cla.pr), 
                                 get.order(pred.stages, cla.pr))
  res[['eval']][['other']] <- confusionMatrix(pred.stages, actual.stages)

  #prob for 2nd class or the one that is imbalanced
  if(cv)
  {
    res[['eval']][['auc']] <- mean(unlist(sapply(pred.stages.score$gr, function(g)
      {
      if(length(unique(as.character(actual.stages[g]))) == 1)
        return(NULL)
      else
        return(auc(get.order(actual.stages[g], cla.pr), pred.stages.score$pr[g]))
      
    })))
    
     res[['eval']][['pr_auc']] <- mean(unlist(sapply(pred.stages.score$gr, function(g)
    {
      if(length(unique(as.character(actual.stages[g]))) == 1)
        return(NULL)
      else
      {
        if(classifier == 'svm')
          pr.curve(scores.class0 = -1*pred.stages.score$pr[g],
                   weights.class0 = get.order(actual.stages[g], cla.rem))$auc.integral
        else
          pr.curve(scores.class0 = 1-pred.stages.score$pr[g],
                   weights.class0 = get.order(actual.stages[g], cla.rem))$auc.integral
      }
    })))
  }
  else
  {
    res[['eval']][['auc']] <- auc(get.order(actual.stages, cla.pr),
                                  pred.stages.score$pr)
    if(classifier == 'svm')
      res[['eval']][['pr_auc']] <- pr.curve(scores.class0 = -1*pred.stages.score$pr,
                    weights.class0 = get.order(actual.stages, cla.rem))$auc.integral
    else
      res[['eval']][['pr_auc']] <- pr.curve(scores.class0 = 1-pred.stages.score$pr,
          weights.class0 = get.order(actual.stages, cla.rem))$auc.integral
  }
  # if(classifier == 'svm')
  #   res[['eval']][['pr_auc']] <- pr.curve(scores.class0 = -1*pred.stages.score$pr, 
  #                                         weights.class0 = get.order(actual.stages, cla.rem))$auc.integral
  # else
  #   res[['eval']][['pr_auc']] <- pr.curve(scores.class0 = 1-pred.stages.score$pr, 
  #                                         weights.class0 = get.order(actual.stages, cla.rem))$auc.integral
  return(res)
}

#' Gets the predicted classes from probabilites
get.class.pred <- function(pr, class.pr, class.rem, thr = 0.5)
{
  return(ifelse(pr > thr, class.pr, class.rem))
}

#' Gets the list of evaluation parameters of AFs of a feature selection algorithm for a classifier
#' 
#' @param actual.stages Actual labels/stages for the data
#' @param predict.list List of predicted stages for AFs of a feature selection algorithm for a classifier
#' @param clasifier Classifier algorithm used
#' @param thr threshold used for computing classes
#' @return List of AFs with their evaluation metrics
get.eval.list <- function(actual.stages, predict.list, classifier, thr, cv = F)
{
  eval.list <- lapply(predict.list, function(predicted)
  {
    get.eval(actual.stages = actual.stages, pred.stages.score = predicted, 
             classifier = classifier, thr = thr, cv = cv)
  })
  return(eval.list)
}

#' Gets the results object for a list
#' 
#' @param cv.req List containing list/s of cv of various kernels
#' @param stages Vector containing stages the samples on which CV was performed
#' @param child.name Vector containing the children names indicating list to be passed to get.eval function
#' @param classifier Classifier for which CV has to be done
#' @param thr Threshold to be used for computing Labels
#' @return List containing list/s of results object
get.res.cv.req <- function(cv.req, stages, child.name, classifier, thr)
{
  res <- list()
  names.list <- names(cv.req)
  if(!is.null(names.list) && sum(child.name %in% names.list) == length(names.list))
    return(get.eval(actual.stages = stages, pred.stages.score = cv.req, classifier = classifier, 
                    thr = thr, cv = T))
  
  for(i in seq_along(cv.req))
  {
    if(!is.null(names.list))
      res[[names.list[i]]] <- get.res.cv.req(cv.req = cv.req[[names(cv.req)[i]]], stages = stages,
                                             child.name = child.name, classifier = classifier, thr = thr)
    else
      res[[i]] <- get.res.cv.req(cv.req = cv.req[[i]], stages = stages, child.name = child.name, 
                                 classifier = classifier, thr = thr)
  }
  return(res)
}
#' Gets the cross validation results w.r.t various metrics
#' 
#' @param stages.train Vector of stages that was used for training the models
#' @param cv.model object/list output of the function get.cv.model
#' @param cores Number of CPU cores to be used
#' @return List containing the evluation metric w.r.t AFs of all feature selection algorithms for all classifiers
get.cv.res <- function(stages.train, cv.model, cores = 1)
{
  fea.names <- names(cv.model)
  registerDoParallel(cores = cores)
  cv.results <- foreach(i = 1:length(fea.names)) %dopar%
  {
    fea.name <- fea.names[i]
    #print(fea.name)
    res <- list()
    #print(i)
    res[['shrunken']] <- get.eval.list(stages.train, cv.model[[fea.name]]$shrunken$pred, 
                                       classifier = 'shrunken', thr = 0.5, cv = F)
    knn.pred <- lapply(cv.model[[fea.name]][['knn']], function(model)
    {
      model$pred
    })
    res[['knn']] <- get.eval.list(stages.train, knn.pred, classifier = 'knn', thr = 0.5, cv = T)
    
    res[['svm']] <- get.eval.list(stages.train, cv.model[[fea.name]]$svm, classifier = 'svm',
                                  thr = 0, cv = T)
    res[['rf']] <- get.eval.list(stages.train, cv.model[[fea.name]][['rf']], classifier = 'rf',
                                 thr = 0.5, cv = T)
    res[['nb']] <- get.eval.list(stages.train, cv.model[[fea.name]][['nb']], classifier = 'nb',
                                 thr = 0.5, cv = T)
    return(res)
  }
  names(cv.results) <- fea.names
  registerDoSEQ()
  return(cv.results)
}

#' Gets the test results w.r.t various metrics
#' 
#' @param stages.test Vector of actual stages of the testing data
#' @param test.pred object/list output of the function get.test.pred
#' @param cores Number of CPU cores to be used
#' @return List containing the evluation metric w.r.t AFs of all feature selection algorithms for all classifiers
get.test.res <- function(stages.test, test.pred, cores = 1)
{
  fea.names <- names(test.pred)
  #cores = detectCores() - 1
  test.results <- foreach(i = 1:length(fea.names)) %do%
  {
    fea.name <- fea.names[i]
    
    res <- mclapply(names(test.pred[[fea.name]]), function(name){
      thr = 0.5
      if(name == 'svm')
        thr = 0
      
      get.eval.list(actual.stages = stages.test, predict.list = test.pred[[fea.name]][[name]],
                    classifier = name, thr = thr, cv = F)
    }, mc.cores = cores)  
    names(res) <- names(test.pred[[fea.name]])
    return(res)
  }
  names(test.results) <- fea.names
  registerDoSEQ()
  return(test.results)
}

get.best.mcc.pr_aucs <- function(res.list, type = 1)
{
  mcc <- get.mcc(res.list = res.list, type = type)
  pr_aucs <- get.pr_aucs(res.list = res.list, type = type)
  comb.met <- mcc + pr_aucs
  return(which.max(comb.met))
}
