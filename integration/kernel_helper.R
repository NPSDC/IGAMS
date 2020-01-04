library(kernlab)
library(genefilter)
library(abind)
library(doParallel)
library(stats)
library(philentropy)

#' Creates partitions of data given the features list
#' 
#' @param train.data Matrix/Data.frame with rows as samples and features as columns
#' @param test.data Matrix/Data.frame with rows as samples and features as columns
#' @return list containing partitions of test and training data
create.data.partition <- function(train.data, test.data, features.list)
{
  if(ncol(train.data) != ncol(test.data))
    stop('Different number of features in train and test data')
  if(is.null(names(features.list)))
    stop('Enter valid feature list names')
  train.data.list <- lapply(features.list, function(feature) train.data[,feature])
  test.data.list <- lapply(features.list, function(feature) test.data[,feature])
  names(train.data.list) = names(features.list)
  names(test.data.list) = names(features.list)
  return(list('train.data' = train.data.list, 'test.data' = test.data.list))
}

#' Generates kernel function required with parameters (RBF, Polynomial, Jaccard, Gower)
#' 
#' @param kernel.name - Name of the kernel, (rbf, polynomial, Jaccard, Gower)
#' @param kernel.param - Parameter of the kernel function - sigma for RBF, degree for polynomial
#' @return Kernel function computed
create.kernel.function <- function(kernel.name, kernel.param)
{
  if(kernel.name == 'rbf')
  {
    if(!('sigma' %in% names(kernel.param)))
      stop('Provide sigma')
    else
      return(rbfdot(sigma = kernel.param[['sigma']]))
  }
  else if(kernel.name == 'polynomial')
  {
    if(!('degree' %in% names(kernel.param)))
      stop('Provide degree of polynomial')
    else
      return(polydot(degree = kernel.param[['degree']]))
  }
  else if(kernel.name == 'gower')
  {
    f <- function(x1, x2)
      return(1-(distance(rbind(x1, x2), method = 'gower')))
    
    class(f) <- "kernel"
    return(f)
  }
  else if(kernel.name == 'jaccard')
  {
    f <- function(x1, x2)
      return(1-(distance(rbind(x1, x2), method = 'jaccard')))
    
    class(f) <- "kernel"
    return(f)
  }
  else
    stop('Only rbf, polynomial, jaccard, kernel')
}

#' Scales the data in which features are mean centred and divided by standard deviation
#' 
#' @param train.data Matrix/data.frame with rows as samples with columns as features
#' @param test.data Matrix/data.frame with rows as samples with columns as features
#' @param mean Mean value that has to be subtracted
#' @param sd Standard Deviation to be divided
#' @return List containing centred training and test data
scale.data <- function(train.data, test.data, means = NULL, sds = NULL)
{
  if(ncol(train.data) != ncol(test.data))
    stop('Different number of features in train and test data')
  if(is.null(means))
    means = colMeans(train.data)
  if(is.null(sds))
    sds = rowSds(t(train.data))
  
  zero.sds = which(sds <= 1e-10)
  
  if(length(zero.sds) != 0)
  {
    train.data <- train.data[,-zero.sds]
    test.data <- test.data[,-zero.sds]
    sds <- sds[-zero.sds]
    means <- means[-zero.sds]
  }
  train.data.scale <- scale(train.data, center = T, scale = T)
  test.data.scale <- test.data - (rep(1,nrow(test.data)) %*% t(means)) 
  
  for(i in seq_along(sds))
    test.data.scale[,i] <- test.data.scale[,i]/sds[i]
  
  return(list('train.data' = train.data.scale, 'test.data' = test.data.scale))
}

#' Computes kernel matrix given a kernel function
#' 
#' @param train.data Matrix/data.frame with rows as samples with columns as features
#' @param test.data Matrix/data.frame with rows as samples with columns as features
#' @param kern.func Kernel function to be used for computing kernel matrix
#' @param be Logical indicating whether train.test should be returned or test.train
#' @param test Logical indicating whether kernel between test and train to be computed
#' @param unit.norm Logical whether to normalise kernel matrix or not
#' @return List of kernel matrix containing test and training data
create.kernel.matrix <- function(train.data, test.data, kern.func, be = F, test = F, unit.norm = T)
{
  if(ncol(train.data) != ncol(test.data))
    stop('Different number of features in train and test data')
  kern.train <- kernelMatrix(kernel = kern.func, train.data)
  kern.test.train <- NULL
  kern.train.test <- NULL
  if(be)
    kern.train.test <- kernelMatrix(kernel = kern.func, train.data, test.data)
  else
    kern.test.train <- kernelMatrix(kernel = kern.func, test.data, train.data)
  kern.test <- NULL
  if(test)
    kern.test <- kernelMatrix(kernel = kern.func, test.data)
  if(unit.norm)
  {
    w = sum(diag(kern.train))
    if(w > 1e-10)
    {
      kern.train <- kern.train/w
      kern.test <- kern.test/w  
      kern.test.train <- kern.test.train/w
      kern.train.test <- kern.train.test/w
    }
  }
  if(test)
  {
    if(be)
      return(list('train.data' = kern.train, 'train.test.data' = kern.train.test, 'test.data' = kern.test))
    else
      return(list('train.data' = kern.train, 'test.train.data' = kern.test.train, 'test.data' = kern.test))
  }
    
  else
  {
    if(be)
      return(list('train.data' = kern.train, 'train.test.data' = kern.train.test))
    else
      return(list('train.data' = kern.train, 'test.train.data' = kern.test.train))
  }
    
}

#' Gets the final kernel matrices doing scaling and normalisation on raw data
#' 
#' @param train.data.list List of matrix/data.frame with rows as samples with columns as features
#' @param test.data.list List of matrix/data.frame with rows as samples with columns as features
#' @param features.list List of list of features on whose each element individual kernel would be formed coresponding to each training data
#' @param kern.func.list List containing different kernel functions
#' @param kern.hash Vector that acts as a hash table associating index to each kernel function name
#' @param kern.ind.list List of list equal to length of features.list containing indexes of every kernel function to be applied to each partition in each training data
#' @param be Logical indicating whether train.test should be returned or test.train
#' @param test Logical to keep kernel matrix between test and training data
#' @param scale Vector indicating whether individual training data needs to be scaled
#' @param unit.norm Logical whether to normalise kernel matrix or not
#' @return List of Kernel matrices in test and training data 
get.kernel.matrices <- function(train.data.list, test.data.list, features.list, kern.func.list, kern.hash, kern.ind.list, 
                              be = F, test = T, scale = rep(T, length(train.data.list)), unit.norm =  T)
{
  if(typeof(train.data.list) != 'list' | typeof(test.data.list) != 'list')
    stop('Training and test data should be list')
  if(sum(sapply(train.data.list, nrow) == nrow(train.data.list[[1]])) != length(train.data.list))
    stop('Rows not consistent across train data')
  if(length(train.data.list) != length(test.data.list))
    stop('Lists of training and test data should be equal')
  if(length(train.data.list) != length(scale))
    stop('Length of scale not equal to training data list')
  for(i in seq_along(train.data.list))
  {
    if(ncol(train.data.list[[i]]) != ncol(test.data.list[[i]]))
      stop(paste('Columns should be equal in both test and training data in ', i))
    if(sum(colnames(train.data.list[[i]]) == colnames(test.data.list[[i]])) != ncol(train.data.list[[i]]))
      stop('Columns should be equal in test and training in ', i)
  }
  if((length(train.data.list) != length(features.list)) | length(train.data.list) != length(kern.ind.list))
    stop("Length of training data should be same as that of features list or kern.ind.list")
  for(i in seq_along(train.data.list))
  {
    if(is.null(names(features.list[[i]])))
      stop(paste('Enter valid feature list names ', i))
    if(length(features.list[[i]]) != length(kern.ind.list[[i]]))
      stop(paste('Number of partitions and kernels provided not equal ', i))
    if(sum(names(features.list[[i]]) == names(kern.ind.list[[i]])) != length(features.list[[i]]))
      stop(paste('Names in kernel list and features list not equal ', i))
  }
  
  mat.names <- c()
  for(i in seq_along(features.list))
  {
    mat.names <- c(mat.names, unlist(sapply(names(features.list[[i]]), function(n)
    {
      paste("data",i,n, names(kern.hash)[kern.ind.list[[i]][[n]]], sep = '_')
    }, simplify = F)))
  }
  
  # #print(kernels.ind.list)
  req.kernel.mat.train <- array(, dim = c(nrow(train.data.list[[1]]), nrow(train.data.list[[1]]), length(unlist(kern.ind.list))),
                          dimnames = list(NULL, NULL, mat.names))
  if(test)
    req.kernel.mat.test <- array(, dim = c(nrow(test.data.list[[1]]), nrow(test.data.list[[1]]), length(unlist(kern.ind.list))),
                                dimnames = list(NULL, NULL, mat.names))
  if(be)
    req.kernel.mat.train.test <- array(, dim = c(nrow(train.data.list[[1]]), nrow(test.data.list[[1]]), length(unlist(kern.ind.list))),
                               dimnames = list(NULL, NULL, mat.names))
  else
    req.kernel.mat.test.train <- array(, dim = c(nrow(test.data.list[[1]]), nrow(train.data.list[[1]]), length(unlist(kern.ind.list))),
                                       dimnames = list(NULL, NULL, mat.names))
  
  for(ind in seq_along(train.data.list))
  {
    for(fea in names(features.list[[ind]]))
    {
      temp <- list()
      temp[['train.data']] <- train.data.list[[ind]][,features.list[[ind]][[fea]]]
      temp[['test.data']] <- test.data.list[[ind]][,features.list[[ind]][[fea]]]
      if(scale[ind])
      {
        temp <- scale.data(train.data = temp[['train.data']],
                           test.data = temp[['test.data']])
      }
      #print(nrow(train.data.list[[ind]][,features.list[[fea]]]))
      for(kern in names(kern.hash)[kern.ind.list[[ind]][[fea]]])
      {
        req.name <- paste("data", ind, fea, kern, sep = '_')
        #print(req.name)
        kern.mat <- create.kernel.matrix(train.data = temp[['train.data']],
                                                           test.data = temp[['test.data']],
                                                           kern.func = kern.func.list[[kern]],
                                         be = be, unit.norm = unit.norm)
        req.kernel.mat.train[,,req.name] <- kern.mat[['train.data']]
        if(be)
          req.kernel.mat.train.test[,,req.name] <- kern.mat[['train.test.data']]
        else  
          req.kernel.mat.test.train[,,req.name] <- kern.mat[['test.train.data']]
  
        if(test)
          req.kernel.mat.test[,,req.name] <- kern.mat[['test.data']]
      }
    }
  }
  if(test)
  {
    if(be)
      return(list('train.mat' = req.kernel.mat.train, 'train.test.mat' = req.kernel.mat.train.test,
                  'test.mat' = req.kernel.mat.test))
    else
      return(list('train.mat' = req.kernel.mat.train, 'test.train.mat' = req.kernel.mat.test.train, 
                  'test.mat' = req.kernel.mat.test))
  }
    
  else
  {
    if(be)
      return(list('train.mat' = req.kernel.mat.train, 'train.test.mat' = req.kernel.mat.train.test)) 
    else
      return(list('train.mat' = req.kernel.mat.train, 'test.train.mat' = req.kernel.mat.test.train)) 
  }
}

#' Gets the final kernel matrices doing scaling and normalisation on raw data
#' 
#' @param train.data Matrix/data.frame with rows as samples with columns as features
#' @param test.data Matrix/data.frame with rows as samples with columns as features
#' @param features.list List of features on whose each element individual kernel would be formed
#' @param kern.func.list List containing different kernel functions
#' @param kern.hash Vector that acts as a hash table associating index to each kernel function name
#' @param kern.ind.list List equal to length of features.list containing indexes of every kernel function to be applied to each partition
#' @param kern.test.train Logical to keep kernel matrix between test and training data
#' @param scale Logical whether to scale or not
#' @param unit.norm Logical whether to normalise kernel matrix or not
#' @return List of Kernel matrices in test and training data 

get.kernel.matrices.temp <- function(train.data, test.data, features.list, kern.func.list, kern.hash, kern.ind.list, 
                                test = T, scale = T, unit.norm =  T, cores = 1)
{
  if(ncol(train.data) != ncol(test.data))
    stop('Different number of features in train and test data')
  if(is.null(names(features.list)))
    stop('Enter valid feature list names')
  if(length(features.list) != length(kern.ind.list))
    stop('Number of paritions and kernels provided not equal')
  if(sum(names(features.list) == names(kern.ind.list)) != length(features.list))
    stop('Names in kernel list and features list not equal')
  req.kernel.train.mat <- array(,dim=c(nrow(train.data), nrow(train.data), length(unlist(kern.ind.list))))
  req.kernel.test.train.mat <- array(,dim=c(nrow(test.data), nrow(train.data), length(unlist(kern.ind.list))))
  registerDoParallel(cores = cores)  
  mat.list <- foreach(i = 1:length(features.list), .combine = 'c') %dopar%
  {
    fea <- names(features.list)[i]
    fea.req <- intersect(features.list[[fea]], colnames(train.data))
    
    if(length(fea.req) <= 1)
      stop('Incorrect features entered')
      
    temp <- list()
    temp[['train.data']] <- train.data[,fea.req]
    temp[['test.data']] <- test.data[,fea.req]
    
    if(scale)
    {
      temp <- scale.data(train.data = temp[['train.data']],
                         test.data = temp[['test.data']])
    }
    req.list <- list()
    m <- foreach(j=1:length(kern.ind.list[[fea]])) %dopar%
    {
      kern = names(kern.hash)[kern.ind.list[[fea]]][j]
      #req.name <- paste(fea, kern, sep = '_')
      kern.mat <- create.kernel.matrix(train.data = temp[['train.data']], 
                           test.data = temp[['test.data']],
                           kern.func = kern.func.list[[kern]], test = test,
                           unit.norm = unit.norm)
      return(kern.mat)
    }
    return(m)
  }
  registerDoSEQ()
  n.dim <- length(mat.list)
  return(mat.list)
  d <- n.dim
  while(d > 0)
  {
    req.kernel.train.mat[,,d] <- mat.list[[d]][[1]]
    req.kernel.test.train.mat[,,d] <- mat.list[[d]][[2]]
    mat.list[[d]] <- NULL
    d = d - 1
  }
  return(list(req.kernel.train.mat, req.kernel.test.train.mat))
}

##Would be used in cv.group.lasso
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}


get.best.cv.kern <- function(res.kern.list, atts)
{
  if(typeof(atts) == 'list')
    k.atts <- vector(mode = 'list', length = length(res.kern.list))
  else
    k.atts <- vector(mode = 'numeric', length = length(res.kern.list))
  k.vals <- vector(mode = 'numeric', length = length(res.kern.list))
  if(sum(names(res.kern.list[[1]]$res.list) == as.character(atts)) != length(res.kern.list[[1]]$res.list))
    stop('Length of given parameters not equal')
  for(i in seq_along(res.kern.list))
  {
    res.req <- res.kern.list[[i]]$res.list
    res <- filter.list(res.list = res.req, ret.inds = T)
    res.req <- res[[1]]
    
    known.ind <- setdiff(seq_along(atts), res[[2]])
    
    atts.req <- atts[known.ind]
    mcc <- get.mcc(res.req, type = 2)
    pr.auc <- get.pr_aucs(res.req, type = 2)
    #print(paste(length(atts), length(mcc), length(known.ind)))
    
    if(sum(names(mcc) == names(pr.auc)) != length(atts.req))
      stop('MCC not equal to PR AUC')
    
    comb.met <- mcc + pr.auc
    
    if(typeof(atts.req) == 'list')
      k.atts[[i]] <- atts.req[[which.max(comb.met)]]
    else
      k.atts[i] <- atts.req[which.max(comb.met)]
    k.vals[i] <- max(comb.met)
  }
  k.req <- which(k.vals == max(k.vals))
  return(list(k.inds = k.req, atts = k.atts[k.req]))
}

set.same.kernels <- function(features.list, k.list)
{
  kern.list <- vector(mode = 'list', length = length(k.list))
  for(i in seq_along(k.list))
  {
    kern.list[[i]] <- vector(mode = 'list', length = length(features.list))
    for(j in seq_along(features.list))
    {
      kern.list[[i]][[j]] <- vector(mode = 'list', length = length(features.list[[j]]))
      names(kern.list[[i]][[j]]) <- names(features.list[[j]])
      for(fea in names(features.list[[j]]))
        kern.list[[i]][[j]][[fea]] <- k.list[[i]]
    }
  }
  return(kern.list)
}

#' Creates all the permuations/subsets/combinations of T,F for a given n
#' 
#' @param n number of base elements
#' @return List of length 2^n, where each element contains a given combination
get.scale.perm <- function(n)
{
  if(n == 1)
  {
    vals <- list()
    vals[[1]] <- T
    vals[[2]] <- F
    return(vals)
  }
  
  else
  {
    vals <- get.scale.perm(n-1)
    cur.list <- list()
    for(i in seq_along(vals))
      cur.list[[i]] <- c(vals[[i]], T)
    for(j in seq_along(vals))
      cur.list[[i+j]] <- c(vals[[j]], F)
    return(cur.list)
  }
}

create.grid.params <- function(n, k.list)
{
  grid.list <- list()
  count <- 1
  scale.list <- get.scale.perm(n)
  change.params <- c()
  for(i in seq(scale.list))
  {
    for(j in seq_along(k.list))
    {
        grid.list[[count]] <- list("scale" = scale.list[[i]], "k.ind" = k.list[[j]])
        count = count + 1
    }
  }
  return(grid.list)
}
