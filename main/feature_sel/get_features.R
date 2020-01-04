source('main/feature_sel/feature_sel_list.R')

get.meth.fea <- function(train.data, train.stages, train.ind.list, cores = 1, type = 2, min.genes = 2,
                         min.range = 0.02)
{
  features <- list()
  features[['shrunken']] <- get.shrunken.features.group(train.data = train.data, train.stages = train.stages,
                                                train.ind.list = train.ind.list, cores = cores,
                                                type = type, min.genes = min.genes, min.range = min.range)
  features[['varSelRF']] <- get.varSelRf.features.group(train.data = train.data, train.stages = train.stages, 
                                                train.ind.list = train.ind.list, cores = cores, type = type)
  return(features)
}

get.meth.fea.list <- function(train.data.list, train.stages, train.ind.list, cores = 1, type = 2, 
                              min.genes = 2, min.range = 0.02)
{
 features.list <- lapply(train.data.list, function(train.data)
   {
    get.meth.fea(train.data = train.data, train.stages = train.stages, train.ind.list = train.ind.list,
                 cores = cores, type = type, min.genes = min.genes, min.range = min.range)
 })
 return(features.list)
}