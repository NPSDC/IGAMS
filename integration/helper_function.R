#' Gets the indexes for list/vector containing list/vector of feature set(s) for
#'  different feature selection algorithm(s)
#'
#' @param features List/Multi-level List/Vector containing bump feature(s) 
#' @return List/Multi-level List/Vector of indexes w.r.t features
get.features.inds <- function(features)
{
  if(typeof(features) == 'list')
    return(lapply(features, function(fea) get.features.inds(features = fea)))
  return(as.numeric(sapply(strsplit(features, 'V'), function(ind.list) ind.list[2])))
}

#' Extracts CpGs from list containing cpgs for different bumps given list of different indexes of 
#' bumps extracted by different feature selection algorithms or just a 1 way list or a simple vector
#' 
#' @param bmp.cpgs List containing cpgs for each bump
#' @param inds Vector or list(1  or mult-level) of indexes indicating which bmps indexes to be extracted
#' @return List/Vector containing CpGs for each ind
get.inds.list.cpgs <- function(bmp.cpgs, inds)
{
  if(typeof(inds) == 'list')
    return(lapply(inds, function(ind) get.inds.list.cpgs(bmp = bmp.cpgs, inds = ind)))
  return(bmp.cpgs[inds])
}

get.cpgs.genes <- function(cpgs.data)
{
  genes.names <- lapply(strsplit(cpgs.data[,2], ';', fixed = T), unique)
  names(genes.names) <- rownames(cpgs.data)
  return(genes.names)
}

# Breaks down the features into list
reduce.features <- function(fea.list)
{
  fea.req.list <- list()
  for(f in names(fea.list))
  {
    for(f_type in names(fea.list[[f]]))
     fea.req.list[[paste(f, f_type, sep = '_')]] <- fea.list[[f]][[f_type]]
  }
  return(fea.req.list)
}
