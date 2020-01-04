source("main/get_eval_mat.R")
library(matrixStats)
#' Creates random partitions based on total number required and total samples
#' 
#' @param stages Vector with stage distribution
#' @param  num.group The number of partitions that need to be created
#' @param strat Logical determining whether to stratify or not
#' @return a list containing the indexes within with each partition
build.groups <- function(stages, num.group, strat = F, seed = 1000)
{
  gr <- lapply(seq(1:num.group), function(x) x <- c())
  if(strat)
  {
    if(length(levels(as.factor(stages))) == 1)
      stop('Stages cannot be NULL or should contain more than 1 level')
  }
  else
    stages <- rep(1,length(stages))
  i = 1
  for(stage in levels(as.factor(stages)))
  {
    set.seed(i*seed)
    stage.ind = which(stages == stage)
    total.stage.samples <- length(stage.ind)
    stage.dist <- sample(rep(1:num.group, ceiling(total.stage.samples/num.group)), total.stage.samples)
    gr <- lapply(seq(1:num.group), function(i) c(gr[[i]], stage.ind[which(stage.dist == i)]))
    i = i + 1
  }
  gr <- lapply(gr, sort)
  return(gr)
}

#' Gives the stage distribution of a partition/s
#' @param  gr List containing the partitions
#' @param stages Vector containing the overall stages
#' @return A list containing the stage distribution of each sample within the partition 
get.stage.distribution <- function(gr, stages)
{
  stage.dist <- lapply(gr, function(x)
  {
    table(stages[x])
  })
  return(stage.dist)
}

#' Gives the genes present in atleast max no of models from the given genes list
#' 
#' @param genes.list List containing genes across the different groups
#' @param max.no.of.models No of groups should the genes atleast be present in
#' @return final genes that are atleast present in the asked number of groups
get.genes.common <- function(genes.list, max.no.of.models)
{
  total.genes <- Reduce(union, genes.list)
  genes.req <- c()
  for(i in total.genes)
  {
    count = 0
    for(j in seq_along(genes.list))
    {
      if(i %in%  genes.list[[j]])
        count = count + 1
      if(count == max.no.of.models)
      {
        genes.req <- c( genes.req, i)
        break
      }
    }
  }
  return(genes.req)
}

#'Creates the deseq2 res object
#'
#' @param dds_comp deseq object that would be used for differential analysis
#' @param contrast column in colData of dds object that would be used for comparison
#' @param stage.cont  base stage that would be used for comparision
#' @param stages.tum vector of stages that would be compared  against stage.cont
#' @return The list of resut matrix for comparision of stages.tum against stage.cont
comp.res <- function(dds.comp, contrast, stage.cont, stages.tum)
{
  library(qvalue)
  stages.comp.cont <- lapply(stages.tum, function(x)
  {
    res = results(dds.comp, contrast = c(contrast, x, stage.cont))
    indexes = is.na(res[,'pvalue'])
    res = res[!indexes,]
    res$qvalue = qvalue(res$pvalue)$qvalues
    res
  })
  names(stages.comp.cont) = stages.tum
  return(stages.comp.cont)
}

#' Gets the stages in 1 or 0 form with 1 being the dominant one
#' 
#' @param stages Vector that needs to be converted in the above form
#' @param type Level of stages that would be the dominant category
#' @param neg Value for negative 0 or 1
#' @return Vector with stages converted into 1 or 0
get.order <- function(stages, type, neg = 0)
{
  return(ifelse(stages == type, 1, neg))
}
#' Gets the genes from result object passing the threshold of log2FC and pvalue
#' 
#' @param res The result object created by deseq2
#' @param logfc log2FC threshold
#' @param adj.pval adjusted pvalue cutoff
#' @return the genes passing the threshold criteria
#' @param pval pvalue cutoff  
get.genes <- function(res, logfc, adj.pval, pval)
{
  inds <- which(abs(res[,2]) > logfc & res[,6] < adj.pval & res[,5] < pval)
  if(length(inds) > 0)
    return(rownames(res)[inds])
  else
    return(NULL)
#  return(rownames(res)[abs(res[,2]) > logfc & res[,6] < adj.pval & res[,5] < pval])
}

#' Gets the Capitalised full Classifier Name
#' 
#' @param cla classifier name to be converted
#' @param return updated classifier name
get.class.name <- function(cla)
{
  cla <- tolower(cla)
  if(cla == 'knn')
    return(c('KNN'))
  else if(cla == 'nb')
    return(c('NB'))
  else if(cla == 'rf')
    return(c('RF'))
  else if(cla == 'shrunken')
    return(c('SC'))
  else if(cla == 'svm')
    return(c('SVM'))
  else if(cla == 'gr.lasso' | cla == 'gr' | cla == 'lasso' | cla == 'gl')
    return(c('GL'))
  else if(cla == 'bemkl')
    return(c('BEMKL'))
}

#' Creates a consolidated data frame containing the value of an evaluation metric w.r.t all AFs for every feature selection w.r.t every classifier
#' 
#' @param test.ac List returned by get.cv.res or get.test.res
#' @return data.frame that can be used by ggplot for drawing plots w.r.t each metric
create.net.df <- function(test.ac)
{
  req.df <- list()
  req.df$roc_auc <- c()
  req.df$accuracy <- c()
  req.df$sens <- c()
  req.df$spec <- c()
  req.df$mcc <- c()
  req.df$pr_auc <- c()
  req.df$f_val <- c()
  req.df$classifier <- c()
  req.df$feature_sel <- c()
  
  for(j in seq_along(test.ac))
  {
    test.ob <- test.ac[[j]]
    for(i in seq_along(test.ob))
    {
      req.df$roc_auc <- c(req.df$roc_auc, unlist(get.aucs(test.ob[[i]])))
      req.df$accuracy <- c(req.df$accuracy, unlist(get.accuracy(test.ob[[i]])))
      req.df$sens <- c(req.df$sens, unlist(get.sens(test.ob[[i]])))
      req.df$spec <- c(req.df$spec, unlist(get.spec(test.ob[[i]])))
      req.df$f_val <- c(req.df$f_val, unlist(get.f(test.ob[[i]])))
      req.df$mcc <- c(req.df$mcc, unlist(get.mcc(test.ob[[i]])))
      req.df$pr_auc <- c(req.df$pr_auc, unlist(get.pr_aucs(test.ob[[i]])))
      req.df$classifier <- c(req.df$classifier,rep(get.class.name(names(test.ob)[i]), 
                                                   length(unlist(get.f(test.ob[[i]])))))
      req.df$feature_sel <- c(req.df$feature_sel,rep(names(test.ac)[j],
                                                     length(unlist(get.f(test.ob[[i]])))))
    }
  }

  df <- data.frame(ROC_AUC <- req.df$roc_auc,  Accuracy <- req.df$accuracy, 
                   Sensitivity <- req.df$sens, Specificity <- req.df$spec,
                   F_Value <- req.df$f_val, MCC <- req.df$mcc, PR_AUC <- req.df$pr_auc,
                   Classifier <- req.df$classifier,
                   Feature_Selection <- req.df$feature_sel)
  colnames(df) <- c('ROC-AUC', 'Accuracy', 'Sensitivity', 'Specificity', 'F1 Score', 'MCC', 'PR-AUC',
                    'Classifier', 'Feature_Selection')
  return(df)
}

create.net.df.one <- function(res.ob, ind = 1, rep, platform)
{
  req.df <- list()
  req.df$auc <- c()
  req.df$accuracy <- c()
  req.df$sens <- c()
  req.df$spec <- c()
  req.df$f_val <- c()
  req.df$mcc <- c()
  req.df$pr_auc <- c()
  req.df$classifier <- c()
  req.df$feature_sel <- c()
  req.df$rep <- c()
  req.df$platform <- c()
  for(j in seq_along(res.ob))
  {
      class.name <- sapply(names(get.aucs(res.ob[[j]])), get.class.name)
      #print(class.name)
      fea.name <- names(res.ob)[j]
      temp.ind = ind
      # if(fea.name != 'varSelRF')
      #   temp.ind = 1
      req.df$auc <- c(req.df$auc, get.aucs(res.ob[[j]], type = 2, ind = temp.ind))
      req.df$accuracy <- c(req.df$accuracy, unlist(get.accuracy(res.ob[[j]], 
                                                                type = 2, ind = temp.ind)))
      req.df$sens <- c(req.df$sens, unlist(get.sens(res.ob[[j]], type = 2, ind = temp.ind)))
      req.df$spec <- c(req.df$spec, unlist(get.spec(res.ob[[j]], type = 2, ind = temp.ind)))
      req.df$f_val <- c(req.df$f_val, unlist(get.f(res.ob[[j]], type = 2, ind = temp.ind)))
      req.df$mcc <- c(req.df$mcc, unlist(get.mcc(res.ob[[j]], type = 2, ind = temp.ind)))
      req.df$pr_auc <- c(req.df$pr_auc, unlist(get.pr_aucs(res.ob[[j]], type = 2, ind = temp.ind)))
      req.df$classifier <- c(req.df$classifier, class.name)
      req.df$feature_sel <- c(req.df$feature_sel, rep(fea.name, length(class.name)))
      req.df$rep <- c(req.df$rep, rep(rep, length(class.name)))
      req.df$platform <- c(req.df$platform, rep(platform, length(class.name)))
  }
  #print(req.df$classifier)
  df <- data.frame(AUC <- req.df$auc,  Accuracy <- req.df$accuracy, 
                   Sensitivity <- req.df$sens, Specificity <- req.df$spec,
                   F_Value <- req.df$f_val, MCC <- req.df$mcc, PR_AUC <- req.df$pr_auc, Classifier <- req.df$classifier,
                   Feature_Selection <- req.df$feature_sel, Representation <- req.df$rep, Platform <- req.df$platform)
  
  colnames(df) <- c('ROC-AUC', 'Accuracy', 'Sensitivity', 'Specificity', 'F_val', 'MCC', 'PR_AUC',
                   'Classifier', 'Feature_Selection', 'Representation', 'Platform')
  
  if(tolower(names(res.ob)[1]) %in% c('sc', 'shrunken', 'shrunken centroids'))
    levels(df[['Feature_Selection']])[1] <- c('Shrunken Centroids')
  return(df)
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
  library(genefilter)
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

merge.group.lasso.res <- function(int.res, gr.res, req.int, mkl.name)
{
  int.res[[length(int.res)+1]] <- list()
  names(int.res)[length(int.res)] <- mkl.name
  for(i in seq_along(req.int))
  {
    if(!is.null(req.int[[i]]))
      int.res[[mkl.name]][names(int.res$shrunken)[i]] <- gr.res[[req.int[[i]]]]
    else
      int.res[[mkl.name]][names(int.res$shrunken)[i]] <- list(NULL)
  }
  return(int.res)
}

remove.dots <- function(ens.ids.all)
{
  ###ens.ids.all <- gets the ids returned from get.genes.files
  
  ##The ens ids contain symbols after dots making them as invalid ensembl ids for using for enrichment
  ##analysis, so stripping the same ids removing the unwanted things after dot
  
  g = sapply(ens.ids.all, function(x) 
  {
    unlist(strsplit(x, split = '.', fixed = T))[1]
  }) ##removing the symbols after .
  return(g)
}

get.list.results.df <- function(res.list, feaSelName, evalQuant, int = F, round = 2)
{
    nrows <- length(res.list)
    if(int)
      names.class <- names(res.list[[1]][["test.res"]][[feaSelName]])
    else
      names.class <- names(res.list[[1]][[feaSelName]])
    df <- vector(mode = "list", length = length(names.class))
    names(df) <- names.class
    for(i in seq_along(df))
      df[[i]] <- rep(0, nrows)
    #print(feaSelName)
    
    df <- data.frame(df)
    for(i in seq(nrows))
    {
      if(int)
        res.list[[i]] <- res.list[[i]][["test.res"]]
      if(evalQuant == "mcc")
        df[i,] <- round(get.mcc(res.list[[i]][[feaSelName]], type = 2), round)
      if(evalQuant == "pr_auc")
        df[i,] <- round(get.pr_aucs(res.list[[i]][[feaSelName]], type = 2), round)
    }
    df <- rbind(df, round(colMeans(df), round))
    df <- rbind(df, round(colSds(as.matrix(df[1:nrows,])), round))
    rownames(df) <- c(paste("G", seq(nrows), sep = ""), "mean", "sd")
    return(df)
}