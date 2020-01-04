source('main/feature_sel/shrunken_features.R')
source('main/feature_sel/varselrf_features.R')
source('main/feature_sel/deseq_features.R')
source('main/feature_sel/sam_features.R')
source('main/feature_sel/dmp_features.R')
source('main/helper_func.R')
#' Gets the features w.r.t each group for shrunken centroid
#' 
#' @param train.data normalised data from which features will be extracted. Rows are samples and columns are features
#' @param train.stages stage of each sample in data
#' @param train.ind.list List containing indexes w.r.t each group in training data
#' @param cores number of cores passed by user for computing groups in parallel
#' @param type Indicates whether train.ind.list consists of single or multiple elements(type = 1)
#' @return List containing the shrunken object as a result of pamr function and combined features
get.shrunken.features.group <- function(train.data, train.stages, train.ind.list, cores = 1, type = 1, 
                                        min.genes = 1, min.range = 0.02)
{
  shrunken.features <- list()
  no.groups <- length(train.ind.list)
  shrunken.features[['genes.object']] <- get.shrunken.object(data = train.data, train.ind.list = train.ind.list,
                                                             stages = train.stages, cores = cores, type = type,
                                                             min.genes = min.genes, min.range = min.range)
  shrunken.features[['genes']] <- get.genes.shrunken(shrunken.features$genes.object)
  groups.name <- paste(rep('atleast', no.groups), seq(no.groups), sep = '_')
  shrunken.features[groups.name] <- lapply(seq(no.groups), function(group)
    {
    get.genes.common(shrunken.features[['genes']], group)
  })
  print('Shrunken Completed')
  return(shrunken.features)
}

#' Gets the features w.r.t each group for VarSelRF
#' 
#' @param train.data normalised data containing samples as rows and columns as genes
#' @param train.stages Stage of every sample in data
#' @param train.ind.list List containing indexes of the folds
#' @param min minimum number of features that would be selected by varSelRF
#' @param type Indicates whether there are 1 or multiple(type = 1) folds in train.ind.list
#' @return List if gene object for all groups
get.varSelRf.features.group <- function(train.data, train.stages, train.ind.list, min = 10, cores = 1, type = 1)
{
  varselRF.features <- list()
  #writeLines(c(""), "log.txt")
  #sink("log.txt", append = T)
  varselRF.features[['genes.object']] <- get.varselRF.object(data = train.data, 
                                                             train.ind.list = train.ind.list, 
                                                             stages = train.stages, cores = cores,
                                                             type = type)
  varselRF.features[['genes.list']] <- get.min.oob.varselRf(varselRf.ob.list = varselRF.features$genes.object, min.fea = min)
  no.groups <- length(train.ind.list)
  groups.name <- paste(rep('atleast', no.groups), seq(no.groups), sep = '_')
  varselRF.features[groups.name] <- lapply(seq(no.groups), function(group)
  {
    get.genes.common(varselRF.features[['genes.list']], group)
  })
  print('VarSelRF Completed')
  return(varselRF.features)
}

#' Gets the features w.r.t each group for DESeq2
#' 
#' @param data count data containing samples as rows and columns as genes
#' @param train.ind.list List containing indexes of the folds
#' @param stages Stage of every sample in data
#' @param cores Number of CPU cores to be used
#' @param type Indicates whether there are 1 or multiple(type = 1) folds in train.ind.list
#' @return List containing the result object as a result of deseq2 function, and AFs
get.deseq2.features.group <- function(train.data, train.ind.list, train.stages, cores = 1, type = 1)
{
  dds.obj <- create.Deseq2(counts.data = train.data, train.ind.list, stages = train.stages,
                           cores = cores, type = type)
  
  features.deseq2 <- list()
  features.deseq2[['genes.object']] <- do.Deseq2(dds.list = dds.obj, cores = cores)
  features.deseq2[['genes.list']] <- get.deseq2.genes(res.list = features.deseq2[['genes.object']], 
                                                      fold.list = c(1,1.5,2))
  
  no.groups <- length(train.ind.list)
  groups.name <- paste(rep('atleast', no.groups), seq(no.groups), sep = '_')
  features.deseq2[groups.name] <- lapply(seq(no.groups), function(group)
  {
    sapply(features.deseq2[['genes.list']], function(x)
      {
      get.genes.common(x,group)
    })
  })
  print('Finished DESeq2')
  return(features.deseq2)
  
}

#' Gets the features w.r.t each group for SAMSeq
#' 
#' @param data count data containing samples as rows and columns as genes
#' @param train.ind.list List containing indexes of the folds
#' @param stages Stage of every sample in data
#' @param cores Number of CPU cores to be used
#' @param type Indicates whether there are 1 or multiple(type = 1) folds in train.ind.list
#' @return The list containing the result object as a result of SAMseq function and AFs
get.sam.features.group <- function(train.data, train.ind.list, train.stages, cores = 1, type = 1)
{
  res.sam <- get.sam.object(data = train.data, train.ind.list = train.ind.list, 
                            stages = train.stages, cores = cores, type = type)
  res.sam.df <- get.sam.df(res.sam = res.sam)
  features.sam <- list()
  features.sam[['sam.obj']] <- res.sam
  features.sam[['genes.list']] <- get.sam.genes(res.sam.df = res.sam.df, folds = c(1, 1.5, 2))
  
  no.groups <- length(train.ind.list)
  groups.name <- paste(rep('atleast', no.groups), seq(no.groups), sep = '_')
  features.sam[groups.name] <- lapply(seq(no.groups), function(group)
  {
    sapply(features.sam[['genes.list']], function(x)
    {
      get.genes.common(x,group)
    })
  })
  print('Finished SAMSeq')
  return(features.sam)
}

#' Gets the features w.r.t each group for SAMSeq
#' 
#' @param data count data containing samples as rows and columns as genes
#' @param train.ind.list List containing indexes of the folds
#' @param stages Stage of every sample in data
#' @param cores Number of CPU cores to be used
#' @param type Indicates whether there are 1 or multiple(type = 1) folds in train.ind.list
#' @return The list containing the result object as a result of SAMseq function and AFs
get.dmp.features.group <- function(train.data, train.ind.list, train.stages, cores = 1, type = 1)
{
  features.dmp <- list()
  features.dmp[['genes.object']] <- get.dmp.df(train.data = train.data, train.ind.list = train.ind.list, 
                                               stages.train = train.stages, cores = cores, type = type)
  features.dmp[['genes.list']] <- get.genes.dmp(dmp.list = features.dmp[['genes.object']], q.val = 1e-10, beta.val = 0)
  
  no.groups <- length(train.ind.list)
  groups.name <- paste(rep('atleast', no.groups), seq(no.groups), sep = '_')
  features.dmp[groups.name] <- lapply(seq(no.groups), function(group)
  {
      get.genes.common(features.dmp[['genes']],group)
    
  })
  print('Finished DMP')
  return(features.dmp)
}

#' Gets the features w.r.t each group for all the four selection algorithms namely DESeq2, SAMSeq, Shrunken and VarSelrf
#' 
#' @param train.count.data count data containing samples as rows and columns as genes
#' @param train.norm.data normalised data containing samples as rows and columns as genes
#' @param train.stages Stage of every sample in data
#' @param train.ind.list List containing indexes of the folds
#' @param cores Number of CPU cores to be used
#' @param type Indicates whether there are 1 or multiple(type = 1) folds in train.ind.list
#' @return List containing the resulting features from individual feature selection algorithms
get.features.object <- function(train.count.data, train.norm.data, train.stages, 
                                train.ind.list, cores.list = c(1,1,1,1) , type = 1)
{
  net.features <- list()
  ######Shrunken####
  net.features[['shrunken']] <- get.shrunken.features.group(train.data = train.norm.data, train.ind.list = train.ind.list,
                                                            train.stages = train.stages,
                                                            cores = cores.list[1], type = type)
  #####VarSelRF####
  net.features[['varSelRF']] <- get.varSelRf.features.group(train.data = train.norm.data, train.ind.list = train.ind.list,
                                                            train.stages = train.stages,
                                                            cores = cores.list[2], type = type)
  ######DESeq2
    net.features[['deseq2']] <- get.deseq2.features.group(train.data = train.count.data, 
                                                          train.ind.list = train.ind.list,
                                                          train.stages = train.stages,
                                                          cores = cores.list[3], type = type)
  ####SAMSeq
  net.features[['sam']] <- get.sam.features.group(train.data = train.count.data,
                                                  train.ind.list = train.ind.list,
                                                  train.stages = train.stages,
                                                  cores = cores.list[4], type = type)
  return(net.features)
}

#' Transforms the nested fold specific features into individual categories
#' 
#' @param features.list List containing only the features that need to be transformed
#' @param folds.list Folds which need to be transformed into individual groups
#' @return List containing folds which need to be transformed
get.filter.fea <- function(features.list, folds.list)
{
  filter.fea <- list()
  for(i in seq_along(folds.list))
  {
    filter.fea[[folds.list[i]]] <- lapply(features.list, function(genes.list)
    {
      if(class(genes.list) == 'matrix')
        genes.list[, folds.list[[i]]]
      else
        genes.list[[folds.list[i]]]  
    })
  }
  return(filter.fea)
}

#' Gets only the features of each feature selection algorithm removing additional objects
#' 
#' @param net.fea List which is output of get.features.object
#' @return List which contains only the features w.r.t each AF
get.class.fea <- function(net.fea)
{
  fea.list <- list()
  for(i in seq_along(net.fea))
  {
    fea.name <- names(net.fea)[i]
    print(fea.name)
    if(fea.name %in% c('shrunken','varSelRF','dmp'))
      fea.list[[fea.name]] <- net.fea[[fea.name]][c(3:length(net.fea[[fea.name]]))]
    else
    {
      temp.fea <- get.filter.fea(net.fea[[fea.name]][c(3:length(net.fea[[fea.name]]))], paste0(c(1,1.5,2), ' fold'))
      for(j in seq_along(temp.fea))
        fea.list[[paste0(fea.name, names(temp.fea)[j])]] <- temp.fea[[names(temp.fea)[j]]]
    }
  }
  return(fea.list)
}
