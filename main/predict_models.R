library(pamr)
library(e1071)
library(class)
library(randomForest)
library(doParallel)
library(foreach)

#' Gets the predicted probabilites from a trained model or knn for the higher order
#' 
#' @param test.data
get.pr <- function(test.data,  shr.thresh = NULL, k = NULL, tr.data = NULL,
                   tr.model = NULL, stages.train = NULL)
{
  if(class(tr.model) == 'pamrtrained')  
  {
    if(is.null(shr.thresh))
    {
      print('give threshold index')
      return(-1)
    }
    return(pamr.predict(fit = tr.model, t(test.data), 
                        threshold = tr.model$threshold[shr.thresh],
                        type = 'posterior'))
  }
  else if(class(tr.model) == 'randomForest')  
    return(predict(tr.model, test.data, type = 'prob')[,1])
  else if(class(tr.model) == 'svm')
  {
    p.svm <- predict(tr.model, test.data, decision.values = T)
    p.svm <- attr(p.svm, 'decision.values')[,1]
    if(tr.model$labels[1] == 2)
      p.svm <- p.svm*-1  
    return(p.svm)
  }
  else if(class(tr.model) == 'naiveBayes')
    return(predict(tr.model, test.data, type = 'raw')[,1])
           
  else
  {
    if(is.null(tr.data) | is.null(stages.train) | is.null(k))
    {
      print('Missing arguments for KNN')
      return(-1)
    }      
    k.pred <- knn(train = tr.data, test = test.data, cl = stages.train, k = k, prob = T)
    k.prob <- attr(k.pred, 'prob')
    type <- levels(as.factor(stages.train))[1]
    return(ifelse(k.pred == levels(as.factor(stages.train))[1], k.prob, 1-k.prob))
  }
}


#' Gets the predicted stages of test data for the shrunken centroid classifiers trained on AFs of 
#' a feature selection algorithm
#' 
#' @param train.models List of trained shrunken centroid models on AFs
#' @param genes.list AFs of the feature selection algorithm
#' @param test.data normalised data to be tested with samples in rows and genes as columns
#' @param thresholds.list Vector of thresholds indexes picked by cross validation to be used for 
#'                  selecting the threshold in the trained model
#' @param cores Number of CPU cores to be used
#' @return List of predicted stages of test data for the AFs
predict.shrunken <- function(train.models, genes.list, test.data, thresholds.list, cores = 1)
{
  if(length(train.models) != length(genes.list))
      print("Length of training models not same as length of gene lists")
  
  shrunken.pred <- mclapply(seq_along(genes.list), function(i)
  {
    if(is.null(train.models[[i]]))
      return(NULL)
    pr <- pamr.predict(train.models[[i]], t(test.data[,genes.list[[i]]]), 
                 threshold = train.models[[i]]$threshold[thresholds.list[i]], type = 'posterior')[,1]
    names(pr) <- rownames(test.data)
    list(pr=pr, cla.pr = train.models[[i]][['cla.pr']], 
         cla.rem = train.models[[i]][['cla.rem']])
  }, mc.cores = cores)
  names(shrunken.pred) <- names(genes.list)
  return(shrunken.pred)
}

#' Gets the predicted stages of test data for KNN on AFs of a feature selection algorithm
#' 
#' @param cv.models Cross validation models output of get.cv.model
#' @param genes.list AFs of a feature selection algorithm
#' @param train.data normalised training data with rows as samples and columns as genes
#' @param test.data normalised test data with rows as samples and columns as genes
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return List of predicted stages of test data for the AFs
predict.knn <- function(cv.models, genes.list, train.data, test.data, stages.train, cores = 1)
{
  if(length(cv.models) != length(genes.list))
    print("Length of training models not same as length of gene lists")
  pred.knn <- mclapply(seq_along(genes.list), function(i)
  {
    if(is.null(genes.list[[i]]) | length(genes.list[[i]]) == 1)
      return(NULL)
    k.pr <- knn(train.data[,genes.list[[i]]], test.data[,genes.list[[i]]], stages.train, 
        cv.models[[i]]$best_k, prob = T)
    k.prob <- attr(k.pr, 'prob')
    type <- levels(as.factor(stages.train))[1]
    pr <- ifelse(k.pr == type, k.prob, 1-k.prob)
    names(pr) <- rownames(test.data)
    list(pr=pr, cla.pr = levels(as.factor(stages.train))[1], 
         cla.rem = levels(as.factor(stages.train))[2])
  }, mc.cores = cores)
  names(pred.knn) <- names(genes.list)
  return(pred.knn)
}

#' Gets the predicted stages of test data for different trained models such as RF, 
#' Naive Bayes and SVM on AFs of a feature selection algorithm
#' 
#' @param train.models List of models trained on AFs
#' @param genes.list AFs of the feature selection algorithm
#' @param test.data normalised data to be tested with samples in rows and genes as columns
#' @param cores No of different CPU cores
#' @return List of predicted stages for the models trained on AFs
predict.model <- function(train.models, genes.list, test.data, cores = 1)
{
  if(length(train.models) != length(genes.list))
    print("Length of training models not same as length of gene lists")
  pred <- mclapply(seq_along(genes.list), function(i)
  {
    if(is.null(train.models[[i]]))
      return(NULL)
    pr <- get.pr(tr.model = train.models[[i]], test.data = test.data[,genes.list[[i]]])
    names(pr) <- rownames(test.data)
    list(pr=pr, cla.pr = train.models[[i]][['cla.pr']], 
         cla.rem = train.models[[i]][['cla.rem']])
  }, mc.cores = cores)
  names(pred) <- names(genes.list)
  return(pred)
}

#' Gets the predicted stages of test data for different classifier models on AFs of all 
#' feature selection algorithm
#' 
#' @param tr.data normalised training data with rows as samples and columns as genes
#' @param te.data normalised test data with rows as samples and columns as genes
#' @param fea.list List of AFs of every feature selection algorithm
#' @param stages.tr Stage of every sample in training data
#' @param tr.model Output of get.train.model containing trained models for all AFs
#' @param cv.models Cross validation models output of get.cv.model
#' @param cores No of different CPU cores
get.test.pred <- function(tr.data, te.data, fea.list, stages.tr, tr.model, cv.model, cores = 1)
{
  registerDoParallel(cores = cores)
  test.pred <- foreach(i = 1:length(fea.list)) %dopar%
  {
    fea.name <- names(fea.list)[i]
    print(fea.name)
    pred.list <- list()
    pred.list[['shrunken']] <- predict.shrunken(train.models = tr.model[[fea.name]][['shrunken']],
                                                genes.list = fea.list[[fea.name]], test.data = te.data,
                                                thresholds.list = cv.model[[fea.name]][['shrunken']]$thr,
                                                cores = cores)
    pred.list[['rf']] <- predict.model(train.models = tr.model[[fea.name]][['rf']],
                                      genes.list = fea.list[[fea.name]], test.data = te.data, cores = cores)

    pred.list[['nb']] <- predict.model(train.models = tr.model[[fea.name]][['nb']],
                                        genes.list = fea.list[[fea.name]], test.data = te.data, cores = cores)
    pred.list[['svm']] <- predict.model(train.models = tr.model[[fea.name]][['svm']],
                                        genes.list = fea.list[[fea.name]], test.data = te.data, cores = cores)
    pred.list[['knn']] <- predict.knn(cv.models = cv.model[[fea.name]][['knn']],
                                      genes.list = fea.list[[fea.name]], train.data = tr.data,
                                      test.data = te.data, stages.train = stages.tr, cores = cores)
    return(pred.list)
  }
  names(test.pred) <- names(fea.list)
  registerDoSEQ()
  return(test.pred)
}