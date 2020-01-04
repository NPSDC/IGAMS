source('main/helper_func.R')
library(DESeq2)
#' Gets the DESeq2 object w.r.t each group
#' 
#' @param counts.data count data containing samples as rows and columns as genes
#' @param train.ind.list List containing indexes of the folds
#' @param stages Stage of every sample in data
#' @param cores Number of CPU cores to be used
#' @param type Indicates whether there are 1 or multiple(type = 1) folds in train.ind.list
#' @return List containing the deseq2 object as a result of deseq function
create.Deseq2 <- function(counts.data, train.ind.list, stages, cores, type = 1)
{
  sampData <- data.frame(stage=stages, names=rownames(counts.data), row.names = rownames(counts.data))
  deseq.list <- mclapply(seq(length(train.ind.list)), function(x)
  {
    if(type == 1)
      train.ind  <- sort(unlist(train.ind.list[-x]))
    else
      train.ind <- sort(unlist(train.ind.list[[x]]))
    #View(sampData[train.ind,])
    dds.obj <- DESeqDataSetFromMatrix(t(counts.data[train.ind,]),
                                      colData = sampData[train.ind,], design = ~stage)
    dds.obj <- dds.obj[rowSums(assay(dds.obj)) > 2]
    dds.obj <- DESeq(dds.obj)
  }, mc.cores = cores)
  return(deseq.list)
}


#' gets the deseq2 result object for every group
#' 
#' @param dds.list List containing deseq2 objects for every group
#' @param cores Number of CPU cores to be used
#' @return List containing the result data.frame of comparision of stages for every group of deseq2 object 
do.Deseq2 <- function(dds.list, cores)
{
  deseq.res <- mclapply(dds.list, function(dds)
  {
    stages <- levels(as.factor(colData(dds)[,'stage']))
    if(stages %in% c('early', 'late'))
      comp.res(dds, 'stage', 'early', 'late')[['late']]
    else
      comp.res(dds, 'stage', 'stage i', 'stage iv')[['stage iv']]
  }, mc.cores = cores)
  return(deseq.res)
}

#' Gets the deseq2 genes for the respective folds for the result object of every group
#' 
#' @param res.list List containing the result of deseq2 object for every group
#' @param fold.list List containing the threshold log2FC fold changes
#' @return List containing the genes w.r.t every fold and every group
get.deseq2.genes <- function(res.list, fold.list)
{
  genes.list <- lapply(fold.list, function(fold)
  {
    lapply(res.list, function(x)
    {
      get.genes(x, fold, 0.05, 0.05)
    })
  })
  names(genes.list) <- paste(as.character(fold.list), rep(c('fold'), length(fold.list)))
  return(genes.list)
}

