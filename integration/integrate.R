source('../methylation/helper_func.R')
source('build_models.R')
source('CV/cv_models.R')
source('predict_models.R')
source('get_results.R')
source('get_eval_mat.R')

#' Combines the features of mRNA and methylation

#' @param mrna.fea.list List containing mRNA features
#' @param meth.fea.list List containing methylation features
#' @return List of list containing mRNA and methylation features merged
comb.fea <- function(mrna.fea.list, meth.fea.list)
{
  fea.comb.list <- list()
  if(class(mrna.fea.list) != 'list' | class(meth.fea.list) != 'list')
    stop('mrna or methylation features is not list')
  if(is.null(names(mrna.fea.list)) | is.null(names(meth.fea.list)))
    stop('List of mRNA or methylation features is unnamed')
  
  for(mrna in names(mrna.fea.list))
  {
    fea.comb.list[[mrna]] <- list()
    for(meth in names(meth.fea.list))
      fea.comb.list[[mrna]][[meth]] <- union(mrna.fea.list[[mrna]], 
                                                       meth.fea.list[[meth]])
  }
  return(fea.comb.list)
}

merge.data <- function(gene.exp.data, meth.beta.data, gene.fea.list, meth.fea.list, logit = T)
{
  total.rna.fea <- unique(unlist(gene.fea.list))
  total.meth.fea <- unique(unlist(meth.fea.list))
  #print(length(get.class.fea(meth.fea.list)))
  if(!(sum(total.meth.fea %in% colnames(meth.beta.data)) == length(total.meth.fea)))
  {
    print('methylation features not found in the methylation array')
    return(-1)
  }
  if(!(sum(total.rna.fea %in% colnames(gene.exp.data)) == length(total.rna.fea)))
  {
    print('gene features not found in the gene expression array')
    return(-1)
  }
  rna.samp <- get.case.from.sample(rownames(gene.exp.data))
  meth.samp <- get.case.from.sample(rownames(meth.beta.data))
  if(length(rna.samp) != length(meth.samp))
  {
    print("lengths not same")
    return(-1)
  }
  if(length(intersect(rna.samp,meth.samp)) != length(rna.samp))
  {
    print("samples not match")
    return(-1)
  }
  if(logit & (min(meth.beta.data) < 0 | max(meth.beta.data) > 1))
    stop('Beta values not provided')
  
  match.ids <- match(meth.samp, rna.samp)
  if(logit)
    comb.df <- cbind(gene.exp.data[match.ids,total.rna.fea], logit2(meth.beta.data[,total.meth.fea]))
  else
    comb.df <- cbind(gene.exp.data[match.ids,total.rna.fea], meth.beta.data[,total.meth.fea])
  return(comb.df)
}

get.int.cv <- function(merged.train.data, merge.fea, stages.train, folds = 5, cores = 1, top = 2)
{
  merged.tr.mod <- get.train.model(tr.data = merged.train.data, 
                                   fea.list = merge.fea,
                                   stages.train = stages.train, cores = cores)
  merged.cv.mod <- get.cv.model(tr.data = merged.train.data, fea.list = merge.fea, 
                                folds = folds, stages.train = stages.train, 
                                tr.model = merged.tr.mod, cores = cores)
  cv.res <- get.cv.res(stages.train = stages.train, cv.model = merged.cv.mod, cores = cores)
  if(length(cv.res) != length(merge.fea))
    stop('CV is wrong with length')
  compute.top.sum <- function(res, merge.fea)
  {
    meth.fea.length <- length(merge.fea[[1]])
    res.calc <- lapply(res, function(r)
    {
      acc.sum <- list()
      for(i in seq(meth.fea.length))
        acc.sum[[names(merge.fea[[1]])[i]]] <- sum(sort(mapply(sum, get.mcc(r, type = 2, ind = i), get.pr_aucs(r, type = 2, ind = i)), decreasing = T)[1:2])
      acc.sum
    })
    return(res.calc)
  }
  compute.max <- function(res.ca, top = top)
  {
    inds = order(unlist(res.ca), decreasing = T)[1:top]
    fea.names <- list()
    names.res <- names(unlist(res.ca))
    for(i in inds)
    {
      mrna.fea <- strsplit(names.res[i], '.', fixed = T)[[1]][1]
      meth.fea <- strsplit(names.res[i], '.', fixed = T)[[1]][2]
      if(!(mrna.fea %in% names(fea.names)))
        fea.names[[mrna.fea]] <- c()
      fea.names[[mrna.fea]] <- c(fea.names[[mrna.fea]], meth.fea)
    }
    return(fea.names)
  }
  res.calc <- compute.top.sum(res = cv.res, merge.fea = merge.fea)
  fea.names <- compute.max(res.ca = res.calc, top = top)
  return(list(cv.res = cv.res, res.calc = res.calc, fea.names = fea.names, tr.mod = merged.tr.mod))
}

#' Gets the classification results of the merged data for a given gene feature and list of CpG features
#' 
#' @param gene.exp.data List containing training and test datasets of gene expression data
#' @param meth.beta.data List containing training and test datasets of methylation beta data
#' @param gene.fea.list Vector of gene features
#' @param meth.fea.list object obtained by running the net feature selection algorithms
#' @param stages List containing stages in test and training data set
#' @param cores Number of cores
#' @param return List containing CV and test performances
simple.integrate <- function(gene.exp.data, meth.beta.data, gene.fea.list, meth.fea.list, stages, scale = T, top = 2, cv = F, cores = 1)
{
  train.comb <- merge.data(gene.exp.data = gene.exp.data[['train']], meth.beta.data = meth.beta.data[['train']],
                           gene.fea.list = gene.fea.list, meth.fea.list = meth.fea.list)
  test.comb <- merge.data(gene.exp.data = gene.exp.data[['test']], meth.beta.data = meth.beta.data[['test']],
                           gene.fea.list = gene.fea.list, meth.fea.list = meth.fea.list)
  remove(gene.exp.data, meth.beta.data)
  gc()
  if(scale)
  {
    scaled.data <- scale.data(train.data = train.comb, test.data = test.comb)
    train.comb <- scaled.data[['train.data']]
    test.comb <- scaled.data[['test.data']]
  }
  
  get.best.ob <- function(res.ob, fea.names)
  {
    req.res.ob <- list()
    for(f in names(fea.names))
      req.res.ob[[f]] <- lapply(res.ob[[f]], function(x) x[fea.names[[f]]])
    return(req.res.ob)
  }
  
  fea.comb <- comb.fea(meth.fea.list = meth.fea.list, mrna.fea.list = gene.fea.list)
  res.after.cv <- NULL
  if(cv)
  {
    res.after.cv <- get.int.cv(merged.train.data = train.comb, merge.fea = fea.comb, stages.train = stages[['train']], folds = 5, 
                               cores = cores, top = 4)
    
    req.train.mod <- get.best.ob(res.ob = res.after.cv$tr.mod, fea.names = res.after.cv$fea.names)
    fea.comb.req <- lapply(names(res.after.cv[['fea.names']]), function(n)
    {
      fea.comb[[n]][res.after.cv[['fea.names']][[n]]]
    })
    names(fea.comb.req) <- names(res.after.cv[['fea.names']])
  }
  else
  {
    req.train.mod <- get.train.model(tr.data = train.comb, 
                                    fea.list = fea.comb,
                                    stages.train = stages[['train']], cores = cores)
    fea.comb.req <- fea.comb
  }
  
  req.cv.mod <- get.cv.model(tr.data = train.comb, fea.list = fea.comb.req, folds = 5, 
                     stages.train = stages[['train']], tr.model = req.train.mod, cores = cores)
  pr <- get.test.pred(tr.data = train.comb, te.data = test.comb, 
                            fea.list = fea.comb.req, stages.tr = stages[['train']],
                            tr.model = req.train.mod, cv.model = req.cv.mod, cores = cores)
  cv.res <- get.cv.res(stages.train = stages[['train']], cv.model = req.cv.mod, cores = cores)
  test.res <- get.test.res(stages.test = stages[['test']], test.pred = pr, cores = cores)
  return(list(train.model = req.train.mod, cv.model = req.cv.mod, test.pred = pr, cv.res = cv.res, test.res = test.res, 
              res.sel = res.after.cv))
}

merge.data.t <- function(train.data.list, test.data.list, fea.list, scale = rep(T, length(train.data.list)))
{
  for(i in seq_along(train.data.list))
  {
    train.data.list[[i]] <- train.data.list[[i]][,unlist(fea.list[[i]])]
    test.data.list[[i]] <- test.data.list[[i]][,unlist(fea.list[[i]])]
    if(scale[i])
    {
      temp <- scale.data(train.data = train.data.list[[i]], test.data = test.data.list[[i]])
      train.data.list[[i]] <- temp[['train.data']]
      test.data.list[[i]] <- temp[['test.data']]
    }
  }
  train.data <- Reduce(cbind, train.data.list)
  test.data <- Reduce(cbind, test.data.list)
  return(list(tr.data = train.data, te.data = test.data))
}

merge.features <- function(lis)
{
  perms <- expand.grid(lapply(lis, seq))
  req.fea <- list()
  for(i in seq(nrow(perms)))
  {
    names.fea <- c()
    fea <- c()
    for(j in seq(ncol(perms)))
    {
      cur.name <- names(lis[[j]])[perms[i,j]]
      names.fea <- c(names.fea, cur.name)
      fea <- c(fea, lis[[j]][[cur.name]])
    }
    names.fea <- paste(names.fea, collapse = '_')
    req.fea[[names.fea]] <- fea
  }
  return(req.fea)
}
#' Gets the classification results of the merged data for a given gene feature and list of CpG features
#' 
#' @param gene.exp.data List containing training and test datasets of gene expression data
#' @param meth.beta.data List containing training and test datasets of methylation beta data
#' @param gene.fea.list Vector of gene features
#' @param meth.fea.list object obtained by running the net feature selection algorithms
#' @param stages List containing stages in test and training data set
#' @param cores Number of cores
#' @param return List containing CV and test performances
simple.integrate.t <- function(train.data.list, test.data.list, fea.list, stages, scale = rep(T, length(train.data.list)), top = 2, cv = F, C_set = NULL, n.trees = NULL, cores = 1)
{
  if(sum(sapply(list(train.data.list, test.data.list, fea.list), class) == 'list') != 3)
    stop('Data or input features are not list')
  if((length(train.data.list) != length(test.data.list)) | (length(train.data.list) != length(fea.list)))
    stop('Length of training data not equal to test data or features')
  if(sum(sapply(train.data.list, function(x) dim(x)[1]) == dim(train.data.list[[1]])[1]) != length(train.data.list))
    stop('Training samples not equal across the datasets')
  if(sum(sapply(test.data.list, function(x) dim(x)[1]) == dim(test.data.list[[1]])[1]) != length(test.data.list))
    stop('Test samples not equal across the datasets')
  if(sum(mapply(function(x,y){
    sum(colnames(x) %in% colnames(y)) == ncol(y) & sum(colnames(x) %in% colnames(y)) == ncol(x)
  }, train.data.list, test.data.list)) != length(train.data.list))
    stop('Colnames not equal')
  
  for(i in seq(fea.list))
  {
    if(sum(unlist(fea.list[[i]]) %in% colnames(train.data.list[[i]])) != length(unlist(fea.list[[i]])))
      stop(paste("Features in",i,'missing'))
  }
  
  merged.data <- merge.data.t(train.data.list = train.data.list, test.data.list = test.data.list, fea.list = fea.list, scale = scale)
  remove(train.data.list, test.data.list)
  gc()
  
  get.best.ob <- function(res.ob, fea.names)
  {
    req.res.ob <- list()
    for(f in names(fea.names))
      req.res.ob[[f]] <- lapply(res.ob[[f]], function(x) x[fea.names[[f]]])
    return(req.res.ob)
  }

  fea.comb <- merge.features(lis = fea.list)
  res.after.cv <- NULL
  if(cv)
  {
    res.after.cv <- get.int.cv(merged.train.data = merged.data[['tr.data']], merge.fea = fea.comb, stages.train = stages[['train']], folds = 5,
                               cores = cores, top = 4)
    req.train.mod <- get.best.ob(res.ob = res.after.cv$tr.mod, fea.names = res.after.cv$fea.names)
    fea.comb.req <- lapply(names(res.after.cv[['fea.names']]), function(n)
    {
      fea.comb[[n]][res.after.cv[['fea.names']][[n]]]
    })
    names(fea.comb.req) <- names(res.after.cv[['fea.names']])
  }
  else
  {
    fea.comb <- list(f=fea.comb)
    print(length(fea.comb))
    print(names(fea.comb))
    req.train.mod <- get.train.model(tr.data = merged.data[['tr.data']],
                                     fea.list = fea.comb, C_set = C_set, n.trees = n.trees,
                                     stages.train = stages[['train']], cores = cores)
    fea.comb.req <- fea.comb
  }

  req.cv.mod <- get.cv.model(tr.data = merged.data[['tr.data']], fea.list = fea.comb.req, folds = 5,
                             stages.train = stages[['train']], tr.model = req.train.mod, cores = cores)
  pr <- get.test.pred(tr.data = merged.data[['tr.data']], te.data = merged.data[['te.data']],
                      fea.list = fea.comb.req, stages.tr = stages[['train']],
                      tr.model = req.train.mod, cv.model = req.cv.mod, cores = cores)
  cv.res <- get.cv.res(stages.train = stages[['train']], cv.model = req.cv.mod, cores = cores)
  test.res <- get.test.res(stages.test = stages[['test']], test.pred = pr, cores = cores)
  return(list(train.model = req.train.mod, cv.model = req.cv.mod, test.pred = pr, cv.res = cv.res, test.res = test.res,
              res.sel = res.after.cv))
}