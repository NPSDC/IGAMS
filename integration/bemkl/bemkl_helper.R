source('integration/bemkl/bemkl_supervised_classification_variational_train.R')
source('integration/bemkl/bemkl_supervised_classification_variational_test.R')
source('integration/kernel_helper.R')
source('main/helper_func.R')
source('main/get_eval_mat.R')
source('main/get_results.R')
library(doParallel)
#' Gets the predicted scores on n-fold CV for BEMKL for the combination of kernels
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
#' @param scale Logical whether to scale or not
#' @param unit.norm Logical whether to normalise kernel matrix or not
#' @param cores Number of cores
#' @return List containing scores for each group of partitions
cv.bemkl <- function(ker.mat.list = NULL, train.data.list, folds, stages.train, features.list, kern.func.list, kern.hash,
                           kern.ind.list, params, gr.list = NULL, n = 1, scale = T, unit.norm = T, cores = 1)
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
    #    total.samp <- nrow(train.data)
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
                                        be = T, scale = scale, unit.norm = unit.norm)
      else
        kern.mat <- ker.mat.list[[iter]][[i]]
      
      
      state <- bemkl_supervised_classification_variational_train(Km = kern.mat[['train.mat']], 
                                                                 y = y[train.index], parameters = params)
      pred.prob <- bemkl_supervised_classification_variational_test(Km = kern.mat[['train.test.mat']],
                                                                   state = state)$p[,1]
      print(pred.prob)
    }
    
    pr.prob[unlist(gr)] <- pr.prob
    
    return(list('pr' = pr.prob, gr = gr, cla.pr = levels(as.factor(stages.train))[1],
                cla.rem = levels(as.factor(stages.train))[2]))
  }
  registerDoSEQ()
  return(cv.prob.pred)
}

cv.bemkl.set <- function(train.data.list, stages.train, features.list,
                         gamma.prior, kern.func.list, kern.hash, kern.ind.list, params, folds,
                         gr.list, n = 1, scale = T, unit.norm = T, cores = 1)
{
  if(typeof(gamma.prior) != 'list')
    stop("Gamma prior not list")
  if(sum(sapply(gamma.prior, length) == 2) != length(gamma.prior))
    stop("Gamma prior element should be of length 2")
  
  stages.train <- as.factor(stages.train)
  ##Cross Validation
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
                                                    features.list = features.list, kern.func.list = kern.func.list, 
                                                    kern.hash = kern.hash, kern.ind.list = kern.ind.list, 
                                                    test = F, be = T, scale = scale, unit.norm = unit.norm)
    }
  }
  print('Kernel matrices computed')
  #registerDoMC(cl)
  registerDoParallel(cores = cores)
  
  cv.param.res <- foreach(i=1:length(gamma.prior), .combine = 'comb', .multicombine = T, .init = list(list(), list())) %dopar%
  {
    params$alpha_gamma = gamma.prior[[i]][1]
    params$beta_gamma = gamma.prior[[i]][2]
    
    cv.pred <- cv.bemkl(ker.mat.list = ker.mat.list, train.data.list = train.data.list, folds = folds, stages.train = stages.train,
                        features.list = features.list, kern.func.list = kern.func.list, gr.list = gr.list,
                        kern.hash = kern.hash, kern.ind.list = kern.ind.list,
                        params = params, n = n, scale = scale, unit.norm = unit.norm, cores = 5)
    cv.res <- get.res.cv.req(cv.req = cv.pred, stages = as.factor(stages.train), child.name = c('pr', 'gr', 'cla.pr', 'cla.rem'),
                   classifier = 'bemkl', thr = 0.5)
    return(list(cv.pred, cv.res))
  }
  gc()
  names(cv.param.res[[1]]) <- as.character(gamma.prior)
  names(cv.param.res[[2]]) <- as.character(gamma.prior)
  return(list('pr.list' = cv.param.res[[1]], 'res.list' = cv.param.res[[2]]))
}

train.bemkl <- function(train.data.list, test.data.list, stages.train, stages.test, features.list, cv.costs.res = NULL,
                              gamma.prior, kern.func.list, kern.hash, kern.ind.list, params, folds, 
                              gr.list, n = 1, scale = T, unit.norm = T, cores = 1)
{
  if(typeof(gamma.prior) != 'list')
    stop("Gamma prior not list")
  if(sum(sapply(gamma.prior, length) == 2) != length(gamma.prior))
    stop("Gamma prior element should be of length 2")
  if(!is.null(cv.costs.res) & length(cv.costs.res) != length(kern.ind.list))
    stop("Kernel index length not equal tp cv.costs.res")
  
  
  stages.train <- as.factor(stages.train)
  stages.test <- as.factor(stages.test)
  
  if(is.null(cv.costs.res))
  {
    cv.costs.res <- vector(mode = 'list', length = length(kern.ind.list))
    for(i in seq_along(kern.ind.list))
    {
      cv.costs.res[[i]] <- cv.bemkl.set(train.data.list = train.data.list, stages.train = stages.train, 
                             features.list = features.list, gamma.prior = gamma.prior, kern.func.list = kern.func.list, kern.hash = kern.hash, 
                             kern.ind.list = kern.ind.list[[i]], params = params, folds = folds, gr.list = gr.list, n = n, scale = scale, unit.norm = unit.norm, 
                             cores = cores)
    }
  }
  
  ker.best <- get.best.cv.kern(res.kern.list = cv.costs.res, atts = gamma.prior)
  ker <- ker.best$k.inds[1]
  
  params$alpha_gamma <- ker.best$atts[[1]][1]
  params$beta_gamma <- ker.best$atts[[1]][2]
  
  
  ##Training
  ker.mat.req <- get.kernel.matrices(train.data.list = train.data.list, test.data.list = test.data.list, 
                                     features.list = features.list, kern.func.list = kern.func.list, 
                                     #kern.hash = kern.hash, kern.ind.list = kern.ind.list[[5]], test = F,
                                     kern.hash = kern.hash, kern.ind.list = kern.ind.list[[ker]], test = F,
                                     scale = scale, unit.norm = unit.norm, be = T)
  
  y.train <- get.order(stages = stages.train, type = levels(as.factor(stages.train))[1], neg = -1)
  state <- bemkl_supervised_classification_variational_train(Km = ker.mat.req[['train.mat']], 
                                                            y = y.train, parameters = params)
  pr.test <- bemkl_supervised_classification_variational_test(Km = ker.mat.req[['train.test.mat']], 
                                                             state = state)$p[,1]
  
  res <- get.eval(actual.stages = stages.test, 
                  pred.stages.score = list('pr' = pr.test, cla.pr = levels(as.factor(stages.test))[1],
                                           cla.rem = levels(as.factor(stages.test))[2]), 
                  classifier = 'bemkl', thr = 0.5, cv = F)
  gc()
  return(list('res' = res, 'pr.test' = pr.test, 'state' = state, 'gamma' = ker.best$atts[[1]], 'ker.ob' = ker.best))
}

train.bemkl.trial <- function(train.data.list, test.data.list, stages.train, stages.test, features.list, cv.costs.res = NULL,
                        gamma.prior, kern.func.list, kern.hash, kern.ind.list, params, folds, 
                        gr.list, n = 1, unit.norm = T, cores = 1)
{
  if(typeof(gamma.prior) != 'list')
    stop("Gamma prior not list")
  if(sum(sapply(gamma.prior, length) == 2) != length(gamma.prior))
    stop("Gamma prior element should be of length 2")
  if(!is.null(cv.costs.res) & length(cv.costs.res) != length(kern.ind.list))
    stop("Kernel index length not equal tp cv.costs.res")
  
  
  stages.train <- as.factor(stages.train)
  stages.test <- as.factor(stages.test)
  k.ind <- set.same.kernels(features.list = features.list, k.list = kern.ind.list)
  
  if(is.null(cv.costs.res))
  {
    
    grid.list <- create.grid.params(n = length(train.data.list), k.list = k.ind)
    registerDoParallel(cores = cores)
    cv.costs.res <- foreach(i = 1:length(grid.list)) %dopar%
    {
      scale <- grid.list[[i]]$scale
      k.ind <- grid.list[[i]]$k.ind
      cv.costs.res[[i]] <- cv.bemkl.set(train.data.list = train.data.list, stages.train = stages.train, 
                                        features.list = features.list, gamma.prior = gamma.prior, kern.func.list = kern.func.list, kern.hash = kern.hash, 
                                        kern.ind.list = k.ind, params = params, folds = folds, gr.list = gr.list, n = n, scale = scale, unit.norm = unit.norm, 
                                        cores = cores)
    }
    gc()
    registerDoSEQ()
  }
  
  inds.best <- get.best.cv.kern(res.kern.list = cv.costs.res, atts = gamma.prior)

  ind <- inds.best$k.inds[1]
  
  scale <- grid.list[[ind]][['scale']]
  k.ind <- grid.list[[ind]][['k.ind']]
  print(scale)
  params$alpha_gamma <- inds.best$atts[[1]][1]
  params$beta_gamma <- inds.best$atts[[1]][2]
  
  
  ##Training
  ker.mat.req <- get.kernel.matrices(train.data.list = train.data.list, test.data.list = test.data.list, 
                                     features.list = features.list, kern.func.list = kern.func.list, 
                                     #kern.hash = kern.hash, kern.ind.list = kern.ind.list[[5]], test = F,
                                     kern.hash = kern.hash, kern.ind.list = k.ind, test = F,
                                     scale = scale, unit.norm = unit.norm, be = T)
  
  y.train <- get.order(stages = stages.train, type = levels(as.factor(stages.train))[1], neg = -1)
  state <- bemkl_supervised_classification_variational_train(Km = ker.mat.req[['train.mat']], 
                                                             y = y.train, parameters = params)
  pr.test <- bemkl_supervised_classification_variational_test(Km = ker.mat.req[['train.test.mat']], 
                                                              state = state)$p[,1]
  
  res <- get.eval(actual.stages = stages.test, 
                  pred.stages.score = list('pr' = pr.test, cla.pr = levels(as.factor(stages.test))[1],
                                           cla.rem = levels(as.factor(stages.test))[2]), 
                  classifier = 'bemkl', thr = 0.5, cv = F)
  gc()
  return(list('res' = res, 'pr.test' = pr.test, 'state' = state, 'gamma' = inds.best$atts[[1]], "grid.best" = grid.list[[ind]]))
}