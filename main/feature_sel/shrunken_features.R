library(pamr)
library(pROC)
library(mccr)
library(foreach)
library(doParallel)
source('main/feature_sel/pamr.listgenes.R')
source('main/get_results.R')
source('main/helper_func.R')
#' Extracts genes from the shrunken object
#' 
#' @param shrunken.genes.df.list list containing data frames of output pamr.listgenes() w.r.t each group
#' @return List containing gene names for each group
get.genes.shrunken <- function(shrunken.genes.df.list)
{
  genes = list()
  for(i in seq_along(shrunken.genes.df.list))
  {
    genes[[i]] = shrunken.genes.df.list[[i]][['genes.list']][,2]
  }
  return(genes)
}

#' Gives the shrunken gene object for given data w.r.t each group
#' 
#' @param data normalised data containing samples as rows and columns as genes
#' @param train.ind.list List containing indexes of the folds
#' @param stages Stage of every sample in data
#' @param cores Number of CPU cores to be used
#' @param type Indicates whether there are 1 or multiple(type = 1) folds
#' @return List containing gene object for all groups and AFs
get.shrunken.object <- function(data, train.ind.list, stages, cores, type = 1, min.genes = 1, min.range = 0.02)
{
  registerDoParallel(cores = cores)
  pamr.genes.list <- foreach(i = 1:length(train.ind.list)) %dopar%
  {
    if(type == 1)
      train.ind <- sort(unlist(train.ind.list[-i]))
    else
      train.ind <- sort(unlist(train.ind.list[[i]]))
    train.model <- pamr.train(list(x = as.matrix(t(data[train.ind,])), 
                                            y = stages[train.ind]))
    cv.model <- pamr.cv(train.model, data = list(x = t(as.matrix(data[train.ind,])),
                                             y = stages[train.ind]), nfold = 5)
    type <- as.factor(stages[train.ind])[1]
    mccs <- sapply(seq_along(cv.model$threshold), function(x)
    {
      mccr(get.order(stages[train.ind], type), get.order(cv.model$yhat[,x], type))
    })
    print(mccs)
    thr.ind = sort(which(mccs == max(mccs)), decreasing = T)[1]
    
    if(min.genes != 1)
      thr.ind <- max(which(mccs > (mccs[thr.ind] - min.range)))  
    
    
    genes.list <- pamr.listgene(train.model, 
                                          data = list(x=as.matrix(t(data[train.ind,])),
                                                      y=stages[train.ind]), 
                                          threshold = cv.model$threshold[thr.ind],
                                          fitcv = cv.model, genenames = T)
    return(list(genes.list = genes.list, tr.model = train.model, thr.ind = thr.ind, cv.mod = cv.model))
  }
  return(pamr.genes.list)
}

##Help in changing shrunken features
see.mccs.shrunken <- function(genes.object, i, tr.list, stage.train, type = 1)
{
  cv.mod <- genes.object[[i]]$cv.mod
  mccs <- sapply(seq_along(cv.mod$threshold), function(x)
  {
    if(type == 1)
      train.ind <- sort(unlist(tr.list[-i]))
    else
      train.ind <- sort(unlist(tr.list[i]))
    mccr(get.order(stage.train[train.ind], 'early'), get.order(cv.mod$yhat[,x], 'early'))
  })
  return(mccs)
}

change.shrunken <- function(shr.obj, new.thrs, tr.list, tr.data, stages.train, type = 1)
{
  no.gr <- length(shr.obj$genes.object)
  for(j in seq(no.gr))
  {
    shr.obj$genes.object[[j]]$thr.ind <- new.thrs[j]
    if(type == 1)
      train.ind <- sort(unlist(tr.list[-j]))
    else
      train.ind <- sort(unlist(tr.list[j]))
    cv.mod <- shr.obj$genes.object[[j]]$cv.mod
    shr.obj$genes.object[[j]]$genes.list <- pamr.listgene(shr.obj$genes.object[[j]]$tr.model, 
                                                        data = list(x=as.matrix(t(tr.data[train.ind,])),
                                                        y=stages.train[train.ind]), 
                                                        threshold = cv.mod$threshold[new.thrs[j]],
                                                        fitcv = cv.mod, genenames = T)
  }
  shr.obj[['genes']] <- get.genes.shrunken(shr.obj$genes.object)
  groups.name <- paste(rep('atleast', no.gr), seq(no.gr), sep = '_')
  shr.obj[groups.name] <- lapply(seq(no.gr), function(group)
  {
    get.genes.common(shr.obj[['genes']], group)
  })
  return(shr.obj)
}
