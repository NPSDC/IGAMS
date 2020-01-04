library(pROC)
library(caret)

filter.list <- function(res.list, ret.inds = F)
{
  inds <- seq(length(res.list))[sapply(res.list, function(x) length(x) == 0)]
  res.list[inds] <- NULL
  if(ret.inds)
    return(list(res.list, inds))
  return(res.list)
}

#' Get the AUC of the AFs within a feature selection algorithm for a classifier
#' 
#' @param res.list List within the results object of AFs of a feature selection algorithm for a classifier
#' @type need to be written
#' @ind need to be written
#' @return Vector of AUCs of the AFs of a feature selection algorithm for a classifier
get.aucs <- function(res.list, type = 1, ind = 1)
{
  res.list <- filter.list(res.list)
  if(type == 1)
    aucs <- sapply(res.list, function(x) x$eval$auc )
  else
      aucs <- sapply(res.list, function(x) x[[ind]]$eval$auc)
  return(aucs)
}

#' Get the MCC of the AFs within a feature selection algorithm for a classifier
#' 
#' @param res.list List within the results object of AFs of a feature selection algorithm for a classifier
#' @type need to be written
#' @ind need to be written
#' @return Vector of MCCs of the AFs of a feature selection algorithm for a classifier
get.mcc <- function(res.list, type = 1, ind = 1)
{
  res.list <- filter.list(res.list)
  if(type == 1)
    mccs <- sapply(res.list, function(x) x$eval$mcc)
  else
    mccs <- sapply(res.list, function(x) if(length(x) > 0) x[[ind]]$eval$mcc)
  return(mccs)
}

#' Get the PR_AUC of the AFs within a feature selection algorithm for a classifier
#' 
#' @param res.list List within the results object of AFs of a feature selection algorithm for a classifier
#' @type need to be written
#' @ind need to be written
#' @return Vector of PR_AUCs of the AFs of a feature selection algorithm for a classifier
get.pr_aucs <- function(res.list, type = 1, ind = 1)
{
  res.list <- filter.list(res.list)
  if(type == 1)
    pr_aucs <- sapply(res.list, function(x) x$eval$pr_auc)
  else
    pr_aucs <- sapply(res.list, function(x) if(length(x) > 0) x[[ind]]$eval$pr_auc)
  return(pr_aucs)
}

#' Get the confusion matrix of the AFs within a feature selection algorithm for a classifier
#' 
#' @param res.list List within the results object of AFs of a feature selection algorithm for a classifier
#' @type need to be written
#' @ind need to be written 
#' @return Vector of confusion matrix of the AFs of a feature selection algorithm for a classifier
get.conf.mat <- function(res.list, type = 1, ind = 1)
{
  res.list <- filter.list(res.list)
  if(type == 1)
    conf.mat <- lapply(res.list, function(x) x$conf_mat)
  else
    conf.mat <- lapply(res.list, function(x) x[[ind]]$conf_mat)
  return(conf.mat)
}

#' Get the accuracy of the AFs within a feature selection algorithm for a classifier
#' 
#' @param res.list List within the results object of AFs of a feature selection algorithm for a classifier
#' @type need to be written
#' @ind need to be written 
#' @return Vector of accuracy of the AFs of a feature selection algorithm for a classifier
get.accuracy <- function(res.list, type = 1, ind = 1)
{
  res.list <- filter.list(res.list)
  if(type == 1)
    acc <- lapply(res.list, function(x) x$eval$other$overall[[1]])
  else
    acc <- lapply(res.list, function(x) x[[ind]]$eval$other$overall[[1]])
  return(acc)
}

#' Get the sensitivites of the AFs within a feature selection algorithm for a classifier
#' 
#' @param res.list List within the results object of AFs of a feature selection algorithm for a classifier 
#' @type need to be written
#' @ind need to be written
#' @return Vector of sensitivites of the AFs of a feature selection algorithm for a classifier
get.sens <- function(res.list, type = 1, ind = 1)
{
  res.list <- filter.list(res.list)
  if(type == 1)
    sens <- lapply(res.list, function(x) x$eval$other$byClass[1])
  else
    sens <- lapply(res.list, function(x) x[[ind]]$eval$other$byClass[1])
  return(sens)
}

#' Get the Specificities of the AFs within a feature selection algorithm for a classifier
#' 
#' @param res.list List within the results object of AFs of a feature selection algorithm for a classifier
#' @type need to be written
#' @ind need to be written
#' @return Vector of Specificities of the AFs of a feature selection algorithm for a classifier
get.spec <- function(res.list, type = 1, ind = 1)
{
  res.list <- filter.list(res.list)
  if(type == 1)
    spec <-  lapply(res.list, function(x) x$eval$other$byClass[2])
  else
    spec <-  lapply(res.list, function(x) x[[ind]]$eval$other$byClass[2])
  return(spec)
}

#' Get the F values of the AFs within a feature selection algorithm for a classifier
#' 
#' @param res.list List within the results object of AFs of a feature selection algorithm for a classifier 
#' @type need to be written
#' @ind need to be written
#' @return Vector of F values of the AFs of a feature selection algorithm for a classifier
get.f <- function(res.list, type = 1, ind = 1)
{
  res.list <- filter.list(res.list)
  if(type == 1)
    f.val <-  lapply(res.list, function(x) x$eval$other$byClass[7])  
  else
    f.val <-  lapply(res.list, function(x) x[[ind]]$eval$other$byClass[7])  
  return(f.val)
}