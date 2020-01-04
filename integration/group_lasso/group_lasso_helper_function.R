source('integration/group_lasso/group_lasso_multiple_kernel_classification_train.R')
source('integration/group_lasso/group_lasso_multiple_kernel_classification_test.R')
source('integration/group_lasso/solve_classification_svm_mosek.R')
source('integration/group_lasso/classification_helper.R')
source('integration/kernel_helper.R')
source('main/helper_func.R')
source('main/get_eval_mat.R')
source('main/get_results.R')
library(doParallel)
#' Gets the predicted scores on n-fold CV for GL given all the parameters
#'
#' @param ker.mat.list List containing list of training and test matrices w.r.t each group
#' @param train.data.list List of matrix/data.frame training data with samples in rows
#' @param folds Folds for n-fold CV
#' @param stages.train Stages of samples in training data
#' @param features.list List of list of features on whose each element individual kernel would be formed
#' @param kern.func.list List of list containing different kernel functions
#' @param kern.hash Vector that acts as a hash table associating index to each kernel function name
#' @param kern.ind.list List equal to length of features.list containing indexes of every kernel function to be applied to each partition
#' @param params List of parameters to train group lasso svm
#' @param gr.list NULL or List containing distribution for n-fold validation for each iteration
#' @param n Number indicating the times k-fold CV has to be repeated 
#' @param scale Vector indicating whether individual training data needs to be scaled
#' @param unit.norm Logical whether to normalise kernel matrix or not
#' @param cores Number of cores
#' @param new Whether to ksvm for training kernels
#' @return List containing scores for each group of partitions
cv.group.lasso <- function(ker.mat.list = NULL, train.data.list, folds, stages.train, features.list, kern.func.list, kern.hash, 
                   kern.ind.list, params, gr.list = NULL, n = 1, scale = rep(T, length(train.data.list)), unit.norm = T, cores = 1)
{
  if(!is.null(ker.mat.list) & (length(ker.mat.list) != n | sum(sapply(ker.mat.list, length) == folds) != n))
    stop("length of ker.mat.list should be equal to folds")
  if(!is.null(gr.list))
  {
    if(typeof(gr.list) != 'list')
      stop('gr has to be a list')
    if(length(gr.list) != n)
      stop('N is not equal to the number of times k-fold CV has to be repeated')
    if(sum(sapply(gr.list, length) == folds) != n)
      stop('gr length not equal to folds')
  }
  y <- get.order(stages = stages.train, type = levels(as.factor(stages.train))[1], neg = -1)
  registerDoParallel(cores = cores)
  
  
  cv.prob.pred <- foreach(iter = 1:n) %dopar%
  {
  #  total.samp <- nrow(train.data)
    gr <- list()
    if(is.null(gr.list))
      gr <- build.groups(stages = stages.train, num.group = folds, strat = T)
    else
      gr <- gr.list[[iter]]
    
    pr.prob <- foreach(i = 1:folds, .combine = c) %dopar%
    {
      train.index = sort(unlist(gr[-i]))
      test.index = sort(unlist(gr[i]))
      if(is.null(ker.mat.list))
        kern.mat <- get.kernel.matrices(train.data.list = lapply(train.data.list, function(x) x[train.index,]), 
                                        test.data.list = lapply(train.data.list, function(x) x[test.index,]), 
                                        features.list = features.list, kern.func.list = kern.func.list,
                                        kern.hash = kern.hash, kern.ind.list = kern.ind.list, test = F,
                                        scale = scale, unit.norm = unit.norm)
      else
        kern.mat <- ker.mat.list[[iter]][[i]]
      
      pred.prob <- rep(c(1), nrow(ker.mat.list[[iter]][[i]][['test.train.mat']]))
      
      #return(list(kern.mat[['train.mat']],y[train.index]))
      #print(pred.prob)
      # if(new)
      #   state <- group_lasso_multiple_kernel_classification_train_new(Km = kern.mat[['train.mat']],
      #                                                             y[train.index], parameters = params)
      #else
      
      state <- group_lasso_multiple_kernel_classification_train(Km = kern.mat[['train.mat']],
                                                y[train.index], parameters = params)
      
      pred.prob <- group_lasso_multiple_kernel_classification_test(Km = kern.mat[['test.train.mat']],
                                                                   state = state)$f[,1]
      #print(pred.prob)
    }
    #print("Ss")
    if(length(pr.prob) != length(unlist(gr)))
      pr.prob <- NULL
    else
      pr.prob[unlist(gr)] <- pr.prob
    
    return(list('pr' = pr.prob, gr = gr, cla.pr = levels(as.factor(stages.train))[1],
                cla.rem = levels(as.factor(stages.train))[2]))
  }
  registerDoSEQ()
  gc()
  return(cv.prob.pred)
}

#' Gets the predicted scores on n-fold CV for GL for a given combination of kernels and scaling parameters, optimising the costs
cv.group.lasso.cost <- function(ker.mat.list = NULL, train.data.list, stages.train, features.list,
                               costs, kern.func.list, kern.hash, kern.ind.list, params, folds, 
                               gr.list, n = 1, scale = rep(T, length(train.data.list)), unit.norm = T, cores = 1)
  
{
  if(sum(sapply(train.data.list, function(x) dim(x)[1] == dim(train.data.list[[1]])[1])) != length(train.data.list))
    stop('Training data rows not equal')
  stages.train <- as.factor(stages.train)
  ##Cross Validation
  if(is.null(ker.mat.list))
  {
    ker.mat.list <- list()
    for(i in 1:n)
    {
      ker.mat.list[[i]] <- list()
      for(j in 1:folds)
      {
        train.index = sort(unlist(gr.list[[i]][-j]))
        test.index = sort(unlist(gr.list[[i]][j]))
        ker.mat.list[[i]][[j]] <- get.kernel.matrices(train.data.list = lapply(train.data.list, function(x) x[train.index, ]), 
                                                      test.data.list = lapply(train.data.list, function(x) x[test.index, ]), 
                                                      features.list = features.list, 
                                                      kern.func.list = kern.func.list, 
                                                      kern.hash = kern.hash, kern.ind.list = kern.ind.list, 
                                                      test = F, scale = scale, unit.norm = unit.norm)
      }
    }
  }
  
  registerDoParallel(cores = cores)
  
  cv.costs.res <- foreach(i=1:length(costs), .combine = 'comb', .multicombine = T, 
                          .init = list(list(), list())) %dopar%
  {
    params$C = costs[i]
    cv.pred <- cv.group.lasso(ker.mat.list = ker.mat.list, train.data.list = train.data.list, folds = folds, stages.train = stages.train,
                              features.list = features.list, kern.func.list = kernels, gr.list = gr.list,
                              kern.hash = kernels.hash.table, kern.ind.list = kern.ind.list,
                              params = params, n = n, scale = scale, unit.norm = unit.norm, cores = 1)
    cv.res <- get.res.cv.req(cv.req = cv.pred, stages = as.factor(stages.train), child.name = c('pr', 'gr', 'cla.pr', 'cla.rem'),
                             classifier = 'svm', thr = 0)
    return(list(cv.pred, cv.res))
  }
  
  names(cv.costs.res[[1]]) <- as.character(costs)
  names(cv.costs.res[[2]]) <- as.character(costs)
  gc()
  registerDoSEQ()
  return(list('pr.list' = cv.costs.res[[1]], 'res.list' = cv.costs.res[[2]]))
}
train.group.lasso <- function(train.data.list, test.data.list, stages.train, stages.test, features.list, cv.costs.res = NULL,
                              costs, kern.func.list, kern.hash, kern.ind.list, params, folds, 
                              gr.list, n = 1, scale = rep(T, length(train.data.list)), unit.norm = T, cores = 1)
{
  if(length(train.data.list) != length(test.data.list))
    stop('Training data not equal to test data')
  if(sum(sapply(train.data.list, function(x) dim(x)[1] == dim(train.data.list[[1]])[1])) != length(train.data.list))
    stop('Training data rows not equal')
  if(sum(sapply(test.data.list, function(x) dim(x)[1] == dim(test.data.list[[1]])[1])) != length(test.data.list))
    stop('Test data rows not equal')
  if(!is.null(cv.costs.res) & length(kern.ind.list) != length(cv.costs.res))
    stop('Length of kernel index list should be equal to cv costs res')
  
  stages.train <- as.factor(stages.train)
  stages.test <- as.factor(stages.test)
  ##Cross Validation
  if(is.null(cv.costs.res))
  {
    cv.costs.res <- vector(mode = 'list', length = length(kern.ind.list))
    for(i in seq_along(kern.ind.list))
    {
      cv.costs.res[[i]] <- cv.group.lasso.cost(train.data.list = train.data.list, stages.train = stages.train, 
                                   features.list = features.list, costs = costs, kern.func.list = kern.func.list, kern.hash = kern.hash, 
                                   kern.ind.list = kern.ind.list[[i]], params = params, folds = folds, gr.list = gr.list, n = n, scale = scale, unit.norm = unit.norm, 
                                   cores = cores)  
    }
  }
  
  ker.best <- get.best.cv.kern(res.kern.list = cv.costs.res, atts = costs)
  params$C <- ker.best$atts[1]
  ker <- ker.best$k.inds[1]
  print(ker.best)
  #ker <- 5
  #Training
  ker.mat.req <- get.kernel.matrices(train.data.list = train.data.list, test.data.list = test.data.list,
                                     features.list = features.list, kern.func.list = kern.func.list,
                                     kern.hash = kern.hash, kern.ind.list = kern.ind.list[[ker]], test = F,
                                     scale = scale, unit.norm = unit.norm)
  y.train <- get.order(stages = stages.train, type = levels(as.factor(stages.train))[1], neg = -1)
  # if(new)
  #   state <- group_lasso_multiple_kernel_classification_train_new(Km = ker.mat.req[['train.mat']],
  #                                                           y = y.train, parameters = params)
  #else
    state <- group_lasso_multiple_kernel_classification_train(Km = ker.mat.req[['train.mat']],
                                                              y = y.train, parameters = params)
  pr.test <- group_lasso_multiple_kernel_classification_test(Km = ker.mat.req[['test.train.mat']],
                                                             state = state)$f[,1]

  res <- get.eval(actual.stages = as.factor(stages.test),
                  pred.stages.score = list('pr' = pr.test, cla.pr = levels(as.factor(stages.test))[1],
                                            cla.rem = levels(as.factor(stages.test))[2]),
                  classifier = 'svm', thr = 0, cv = F)
  gc()
  return(list('res' = res, 'pr.test' = pr.test, 'state' = state, 'cost' = params$C, 'ker.ob' = ker.best))
}

  
train.group.lasso.trial <- function(train.data.list, test.data.list, stages.train, stages.test, features.list, cv.costs.res = NULL,
                                    costs, kern.func.list, kern.hash, kern.ind.list, params, folds, 
                                    gr.list, n = 1, unit.norm = T, cores = 1)
{
  if(length(train.data.list) != length(test.data.list))
    stop('Training data not equal to test data')
  if(sum(sapply(train.data.list, function(x) dim(x)[1] == dim(train.data.list[[1]])[1])) != length(train.data.list))
    stop('Training data rows not equal')
  if(sum(sapply(test.data.list, function(x) dim(x)[1] == dim(test.data.list[[1]])[1])) != length(test.data.list))
    stop('Test data rows not equal')
  if(!is.null(cv.costs.res) & length(kern.ind.list) != length(cv.costs.res))
    stop('Length of kernel index list should be equal to cv costs res')
  if(length(features.list) != length(train.data.list))
    stop("Length of features list not equal to training data")
  
  stages.train <- as.factor(stages.train)
  stages.test <- as.factor(stages.test)
  
  k.inds <- set.same.kernels(features.list = features.list, k.list = kern.ind.list)
  
  ##Cross Validation
  
  if(is.null(cv.costs.res))
  {
    grid.list <- create.grid.params(n = length(train.data.list), k.list = k.inds)
    cv.costs.res <- vector(mode = 'list', length = length(grid.list))
    registerDoParallel(cores = cores)
    cv.costs.res <- foreach(i = 1:length(grid.list)) %dopar%
    {
      scale <- grid.list[[i]]$scale
      k.ind <- grid.list[[i]]$k.ind
      cv.costs.res[[i]] <- cv.group.lasso.cost(train.data.list = train.data.list, stages.train = stages.train, 
                                              features.list = features.list, costs = costs, kern.func.list = kern.func.list, kern.hash = kern.hash, 
                                              kern.ind.list = k.ind, params = params, folds = folds, gr.list = gr.list, n = n, scale = scale, unit.norm = unit.norm, 
                                              cores = cores)  
    }
    gc()
    registerDoSEQ()
    inds.best <- get.best.cv.kern(res.kern.list = cv.costs.res, atts = costs)
    params$C <- inds.best$atts[1]
    ind <- inds.best$k.inds[1]
    scale <- grid.list[[ind]][['scale']]
    k.ind <- grid.list[[ind]][['k.ind']]
    print(scale)
    ker.mat.req <- get.kernel.matrices(train.data.list = train.data.list, test.data.list = test.data.list,
                                       features.list = features.list, kern.func.list = kern.func.list,
                                       kern.hash = kern.hash, kern.ind.list = k.ind, test = F,
                                       scale = scale, unit.norm = unit.norm)
    y.train <- get.order(stages = stages.train, type = levels(as.factor(stages.train))[1], neg = -1)
    
    state <- group_lasso_multiple_kernel_classification_train(Km = ker.mat.req[['train.mat']],
                                                                y = y.train, parameters = params)
    pr.test <- group_lasso_multiple_kernel_classification_test(Km = ker.mat.req[['test.train.mat']],
                                                               state = state)$f[,1]
    
    res <- get.eval(actual.stages = as.factor(stages.test),
                    pred.stages.score = list('pr' = pr.test, cla.pr = levels(as.factor(stages.test))[1],
                                             cla.rem = levels(as.factor(stages.test))[2]),
                    classifier = 'svm', thr = 0, cv = F)
    gc()
    return(list('res' = res, 'pr.test' = pr.test, 'state' = state, 'cost' = params$C, "inds.best" = inds.best, "grid.best" = grid.list[[ind]]))
  }
}