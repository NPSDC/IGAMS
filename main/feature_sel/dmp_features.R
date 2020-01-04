library(minfi)

#' Gives the samseq result object for given data w.r.t each group
#' 
#' @param train.data Data.frame/Matrix containing methylation beta data containing samples as rows and columns as genes/cpg probes
#' @param train.ind.list List containing indexes of the folds withing the training data
#' @param stages.train Vector with stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @param type Vector indicates whether there are 1 or multiple(type = 1) folds
#' @return List of DMP dfs for all groups
get.dmp.df <- function(train.data, train.ind.list, stages.train, cores = 1, type = 1)
{
  if(length(which(train.data > 1 | train.data < 0)) != 0)
    stop('Input should be beta values between 0 and 1')
  res.dmp <- foreach(i = 1:length(train.ind.list)) %dopar%
  {
    if(type == 1)
      train.ind <- unlist(train.ind.list[-i])
    else
      train.ind <- unlist(train.ind.list[[i]])
    req.data <- t(train.data[train.ind, ])
    
    dmp <- dmpFinder(dat = req.data, 
                     pheno = stages.train[train.ind],
                     type = 'categorical', qCutoff = 0.05)
    ##Calculating mean beta difference
    early.ind <- which(stages.train[train.ind] == 'early')
    late.ind <- which(stages.train[train.ind] == 'late')
    beta.diff <- rowMeans(req.data[rownames(dmp), early.ind]) -
      rowMeans(req.data[rownames(dmp), late.ind])
    dmp <- cbind(dmp, beta.diff)
    return(dmp)
  }
  return(res.dmp)
}

#' Gets the sig
get.genes.dmp <- function(dmp.list, q.val = 0.05, beta.val = 0)
{
  dmp.genes <- lapply(dmp.list, function(df)
  {
    rownames(df)[which((df[,4] <= q.val) & (abs(df[,5]) >= abs(beta.val)))]  
  })
  return(dmp.genes)
}