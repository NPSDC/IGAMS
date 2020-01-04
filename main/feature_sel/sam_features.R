library(foreach)
library(doParallel)
library(samr)
source('main/feature_sel/sam_func_cop.R')
source('main/feature_sel/samr.morefuns.R')

#' Gives the samseq result object for given data w.r.t each group
#' 
#' @param data normalised data containing samples as rows and columns as genes
#' @param train.ind.list List containing indexes of the folds
#' @param stages Stage of every sample in data
#' @param cores Number of CPU cores to be used
#' @param type Indicates whether there are 1 or multiple(type = 1) folds
#' @return List of SAMseq gene objects for all groups
get.sam.object <- function(data, train.ind.list, stages, cores = 1, type = 1)
{
  stages = as.factor(stages)
  y.req = c()
  for(stage in stages)
  {
    if(stage == levels(stages)[1])    
      y.req = c(y.req, 1)
    else
      y.req = c(y.req, 2)
  }
  registerDoParallel(cores = cores)
  res.sam <- foreach(i = 1:length(train.ind.list)) %dopar%
  {
     print(length(train.ind.list))
     if(type == 1)
       train.ind <- unlist(train.ind.list[-i])
     else
       train.ind <- unlist(train.ind.list[[i]])
    req.data <- t(data[train.ind, ])

    ind.to.rem <- which(rowSums(req.data) < 2)
    if(length(ind.to.rem) > 0)
      req.data <- req.data[-ind.to.rem,]  
    res <- SAMseq(x=req.data, y = y.req[train.ind], resp.type = 'Two class unpaired', 
                           genenames = rownames(req.data))
    return(res)
  }
  return(res.sam)
}

#' Creates the data.frame from SAM result object
#' 
#' @param res.sam List of result object for the groups after passing through SAMSeq
#' @return List containing the data frame converted result objects for each group
get.sam.df <- function(res.sam)
{
  sam.df.list <- lapply(res.sam, function(sam.obj)
  {
    res.df <- c()
    if(sam.obj$siggenes.table$ngenes.up == 0)
      res.df <- sam.obj$siggenes.table$genes.lo
    else if(sam.obj$siggenes.table$ngenes.lo == 0)
      res.df <- sam.obj$siggenes.table$genes.up
    else
      res.df <- rbind(sam.obj$siggenes.table$genes.up, sam.obj$siggenes.table$genes.lo)
    data.frame(GeneID = as.character(res.df[,1]), Gene.Name = as.character(res.df[,2]), 
               Score = as.numeric(res.df[,3]), Fold.Change = as.numeric(res.df[,4]),
               q.val = as.numeric(res.df[,5]), stringsAsFactors = F)
  })
  return(sam.df.list)
}

#' Gets the genes for the fold changes and groups
#' 
#' @param res.sam.df List of data frame containing sam res object
#' @param folds log2FC for comparison
#' @return List containing for every group within a fold
get.sam.genes <- function(res.sam.df, folds)
{
  req.genes <- lapply(folds, function(fold)
  {
    genes <- list()
    for(i in seq_along(res.sam.df))
    {
      genes[[i]] <- res.sam.df[[i]]$GeneID[which(abs(log2(res.sam.df[[i]]$Fold.Change)) > fold & 
                                                  res.sam.df[[i]]$q.val < 0.05)]
    }
    genes
  })
  names(req.genes) <- paste0(folds, ' fold')
  return(req.genes)
}
