#' Gets the patient specific id from sample id
#' @param samples Vector of case ids
#' @param type 1 to include tumor determing code, 2 exclude 
#' @return Vector containing the patient ids
get.case.from.sample <- function(samples, type = 1)
{
  cases <- sapply(strsplit(samples, split = '-', fixed = T), function(case)
  {
    if(type == 1)
      paste(case[1], case[2], case[3], case[4], sep = '-')
    else
      paste(case[1], case[2], case[3], sep = '-')
  })
  return(cases)
}


#'Gets the list of data.frame of cpgs ordered by Chromosomes and sorted by position

#' @param ann.df Ilumina array of locations of Cpgs
#' @return List of data frame with cpgs corresponding to each chromosome
create.chr.cpg <- function(ann.df)
{
  chroms <- levels(as.factor(ann.df$chr))
  chr.cpgs.list <- list()
  for(chrom in chroms)
  {
        chr.cpgs.list[[chrom]] <- ann.df[which(ann.df$chr == chrom), -1]
        chr.cpgs.list[[chrom]] <- chr.cpgs.list[[chrom]][order(chr.cpgs.list[[chrom]]$pos),]
  }
  return(chr.cpgs.list)
}

#'Gets the cpgs lying within certain locations for a chromosome
#'
#' @param chr.cpgs.list List containing 
#' @param chrom chromosome for which the cpgs have to be extracted
#' @param start.pos starting position for the cpg sites in location
#' @param end.pos ending position for the cpg sites in location
#' @return A vector of cpgs lying within the mentioned location
get.cpgs <- function(chr.cpgs.list, chrom, start.pos, end.pos)
{
  return(rownames(chr.cpgs.list[[chrom]])[chr.cpgs.list[[chrom]]$pos >= start.pos & 
                                            chr.cpgs.list[[chrom]]$pos <= end.pos])
}

#' Gets the genes that map to CpGs given annotation data frame
#' @param cpgs Vector of CpGs that have to be mapped
#' @param ann.df rowData containing information about probes
#' @return Vector of genes mapping to CpGs 
get.genes.cpgs <- function(cpgs, ann.df)
{
  return(unique(unlist(strsplit(ann.df[,2][match(cpgs, ann.df[,1])], ';', fixed = T))))
}

#' Gets the length of specific index of a list within lists
#' @fea.list List of lists
#' @param ind index within list
#' @return length of list
get.fea.length <- function(fea.list, ind)
{
  return(sapply(fea.list, function(fea) length(fea[[ind]]), USE.NAMES = T))
}

#' Gets genes from varselrf select object for a particular index
#' @param var.gene.ob list returned by varSelRF run whose object has to be used
#' @param ind index within select object from which genes need to be extracted
#' @return Genes/features from that particular index
get.varSelRF.gene <- function(var.gene.ob, ind)
{
  return(unlist(strsplit(as.character(var.gene.ob$selec.history[ind,2]), ' + ', fixed = T)))
}

get.dmp <- function(GRSet, stages, q.val)
{
  dmp <- dmpFinder(dat = assay(GRSet), pheno = stages, type = 'categorical', qCutoff = q.val)
  early.ind <- which(stages == 'early')
  late.ind <- which(stages == 'late')
  beta.diff <- rowMeans(assay(GRSet)[rownames(dmp), early.ind]) -
    rowMeans(assay(GRSet)[rownames(dmp), late.ind])
  dmp <- cbind(dmp, beta.diff)
  return(dmp)
}

get.cpgs.df <- function(paired, cpgs = NULL, genes = NULL, ens = NULL, rem = F, dmp = NULL, cores = 1)
{
  library(parallel)
  if(class(paired$cpgs) == 'factor')
    paired$cpgs <- as.character(paired$cpgs)
  if(!is.null(dmp))
    paired[,'beta.diff'] <- dmp[paired$cpgs,'beta.diff']
  inds.req <- c()
  if(!is.null(cpgs))
  {
    if(length(cpgs) != length(intersect(cpgs, paired$cpgs)))
      print('Cpgs length not equal')
    cpgs = intersect(cpgs, paired$cpgs)
    print(length(cpgs))
    inds.req <- unlist(mclapply(cpgs, function(x) which(paired$cpgs == x), mc.cores = cores))
    if(length(inds.req) == 0)
      inds.req <- NULL
    #print(length(inds.cpgs))
  }
  if(!is.null(genes))
  {
    if(length(genes) != length(intersect(genes,paired$gene.ids)))
      print('Genes length not equal')
    genes = intersect(genes, paired$gene.ids)
    inds.genes <- unlist(mclapply(genes, function(x) which(paired$gene.ids == x), mc.cores = cores))
    if(length(inds.req) == 0)
      inds.req <- inds.genes
    else
      inds.req <- intersect(inds.req, inds.genes)
  }
  if(!is.null(ens))
  {
    if(length(ens) != length(intersect(ens,paired$ens)))
      print('Ensembl length not equal')
    ens.ids = intersect(ens, paired$ens)
    inds.ens <- unlist(mclapply(ens.ids, function(x) which(paired$ens == x), mc.cores = cores))
    if(length(inds.req) == 0)
      inds.req <- inds.ens
    else
      inds.req <- intersect(inds.req, inds.ens)
  }
  if(length(inds.req) == 0)
    return(NULL)
  
  
  if(rem)
    return(paired[-inds.req,])
  return(paired[inds.req,])
}

get.pairings.all <- function(row.df, dmp.T = NULL, dmp.N.T = NULL)
{
  pairings <- list()
  pairings[['ens']] <- c()
  pairings[['cpgs']] <- c()
  pairings[['gene']] <- c()
  pairings[['pos']] <- c()
  for(i in seq(nrow(row.df)))
  {
    print(i)
    cpg <- rownames(row.df)[i]
    gene <- unlist(row.df[i,'gene.ids'])
    g.dot <- which(gene == '.')
    ens <- unlist(row.df[i,'ens'])
    pos <- unlist(row.df[i,'gene.pos'])
    
    if(sum(sapply(list(gene,ens,pos), length) == length(gene)) != 3)
      stop('length not same')
    if(length(g.dot) != 0)
    {
      if(length(g.dot) == length(gene))
        next()
      gene <- gene[-g.dot]
      ens <- ens[-g.dot]
      pos <- pos[-g.dot]
    }
    na.inds <- which(is.na(gene))
    if(length(na.inds) != 0)
    {
      gene <- gene[-na.inds]
      ens <- ens[-na.inds]
      pos <- pos[-na.inds]
      if(length(gene) == 0)
        next()
    }
    pairings[['ens']] <- c(pairings[['ens']], ens)
    pairings[['cpgs']] <- c(pairings[['cpgs']], rep(cpg, length(gene)))
    pairings[['gene']] <- c(pairings[['gene']], gene)
    pairings[['pos']] <- c(pairings[['pos']], pos)
  }
  if(!is.null(dmp.N.T))
    pairings[,'N.T.beta.diff'] <- dmp[pairings$cpgs,'beta.diff']
  if(!is.null(dmp.T))
    pairings[,'T.beta.diff'] <- dmp[pairings$cpgs,'beta.diff']
  return(data.frame(pairings, stringsAsFactors = F))
}

get.cpg.ens <- function(p.df)
{
  cpg.list <- list()
  for(i in seq(nrow(p.df)))
  {
    cpg <- p.df$cpgs[i]
    ens <- p.df$ens[i]
    gene <- p.df$gene[i]
    if(!(cpg %in% names(cpg.list)))
    {
      cpg.list[[cpg]] <- list()
      cpg.list[[cpg]][['gene']] <- c()
      cpg.list[[cpg]][['ens']] <- c()
    }
    cpg.list[[cpg]][['gene']] <- c(cpg.list[[cpg]][['gene']], gene)
    cpg.list[[cpg]][['ens']] <- c(cpg.list[[cpg]][['ens']], ens)
  }
  return(cpg.list)
}

create.count.df <- function(paired.df, prom, unique = F)
{
  if(unique)
  {
    cpgs.unique <- names(table(paired.df$cpgs))[table(paired.df$cpgs) > 1]
    paired.df <- get.cpgs.df(paired = paired.df, cpgs = cpgs.unique, rem = T)
  }
  
  l <- length(unique(paired.df$gene))
  df <- data.frame('Gene' = unique(paired.df$gene), "Promoters" = rep(0,l), "3'UTR" = rep(0,l),
                   "5'UTR" = rep(0,l), "Body" = rep(0,l), "Not Mapped" = rep(0,l), 'Total' = rep(0,l),
                   stringsAsFactors = F, row.names = unique(paired.df$gene))
  colnames(df)[c(3,4,6)] = c("3'UTR", "5'UTR", "Not Mapped")
  for(i in seq(nrow(paired.df)))
  {
    pos <- paired.df[i,'pos']
    gene <- paired.df[i,'gene.ids']
    if(is.na(pos))
    {
      df[gene,'Not Mapped'] = df[gene,'Not Mapped']+1
      #print(gene)
    }
    else if(pos %in% prom)
      df[gene,'Promoters'] = df[gene,'Promoters']+1
    else if(pos %in% colnames(df)[3:5])
      df[gene,pos] = df[gene,pos]+1
    else
      stop(paste(pos, 'Not Found'))
  }
  df[,'Total'] <- rowSums(df[,c(2:6)])
  return(df = df)
}

get.chroms <- function(row.df, old.ann = F, cpgs = NULL)
{
  if(!is.null(cpgs))
    row.df <- row.df[cpgs,]
  if(old.ann)
  {
    chr <- row.df[,'chr']
    chr <- ifelse(chr == 'chrX' | chr == 'chrY', NA, chr)
  }
  else
  {
    l <- strsplit(row.df$CGI_Coordinate, ':', fixed = T)
    chr <- sapply(l, function(x) if(length(x) == 0) return(NA) else if(x == 'chrX' | x == 'chrY') return(NA) else return(x[2]))
  }
  return(chr)
}

get.ens <- function(map.genes, genes, cores = 1)
{
  library(parallel)
  ens.req <- mclapply(genes, function(g) 
  {
    g=tolower(g)
    id = c(which(tolower(map.genes[,2]) == g), which(tolower(map.genes[,3]) == g), which(tolower(map.genes[,4]) == g))
    if(length(id) > 0)
      return(map.genes[id[1],'ensembl_gene_id'])
    else
      return(NA)
  }, mc.cores =  cores)
  return(unlist(ens.req))
}

##rna.data <- Summarized Experiment
##meth.data <- df with genes in columns
get.beta.log.mean.diff <- function(meth.data, rna.data = NULL, meth.matched.samples, meth.all.samples, rna.matched.samples, rna.all.samples, res.rna, ens.map)
{
  get.means <- function(data)
  {
    samp.names <- colnames(data)
    types <- sapply(strsplit(samp.names, '-', fixed = T), function(x) substr(x = x[[4]], start = 1, stop = 2))
    tum.ind <- which(types != '11')
    norm.ind <- which(types == '11')
    
    mean.tum <- rowMeans(data[,tum.ind], na.rm = T)
    mean.norm <- rowMeans(data[,norm.ind], na.rm = T)
    names(mean.tum) <- rownames(data)
    names(mean.norm) <- rownames(data)
    mean.diff <-  mean.tum - mean.norm
    return(list('mean.tum' = mean.tum, 'mean.norm' = mean.norm, 'mean.diff' = mean.diff))
  }
  ##Filtering ens.map
  if(sum(is.na(ens.map)) > 0)
    ens.map <- ens.map[!is.na(ens.map)]
  
  ##Customising ens.map to our needs
  ens.map.req <- ens.map[intersect(rownames(meth.data), names(ens.map))]
  rna.all.means <- NULL
  rna.matched.means <- NULL
  rna.log.fc <- NULL
  if(!is.null(rna.data))
  {
    ens.map.temp <- intersect(ens.map.req, rownames(rna.data))
    inds <- c()
    for(i in seq_along(ens.map.req))
    {
      if(ens.map.req[i] %in% ens.map.temp)
        inds <- c(inds,i)
    }
    ens.map.req <- ens.map.req[inds]
    rna.data <- rna.data[ens.map.req,]
    rna.all.means <- get.means(rna.data[,rna.all.samples])
    rna.matched.means <- get.means(rna.data[,rna.matched.samples])
  }
  if(!is.null(res.rna))
  {
    ens.map.temp <- intersect(ens.map.req, rownames(res.rna))
    inds <- c()
    for(i in seq_along(ens.map.req))
    {
      if(ens.map.req[i] %in% ens.map.temp)
        inds <- c(inds,i)
    }
    ens.map.req <- ens.map.req[inds]
    res.req <- res.rna[ens.map.req,]
    rna.log.fc <- res.req[,'log2FoldChange']
    names(rna.log.fc) <- rownames(res.req)
  }
  
  meth.data <- meth.data[names(ens.map.req),]
  meth.all.means <- get.means(meth.data[,meth.all.samples])
  meth.matched.means <- get.means(meth.data[,meth.matched.samples])
  
  return(list('meth.all.means' = meth.all.means, 'meth.matched.means' = meth.matched.means, 'rna.all.means' = rna.all.means, 'rna.matched.means' = rna.matched.means, 
              'rna.log.fc' = rna.log.fc, 'ens.map.req' = ens.map.req))
}

create.df.meth.rna.plot <- function(ob.list, ens.genes.map, max.beta = NULL, min.beta = NULL)
{
  library(ggplot2)
  meth.mean.diff <- ob.list$meth.matched.mean$mean.diff
  if(is.null(min.beta))
    min.beta <- round(min(meth.mean.diff) + 0.05, 1)
  if(is.null(max.beta))
    max.beta <- round(max(meth.mean.diff) - 0.05, 1)
  
  vals.range <- round(seq(min.beta, max.beta, 0.1),1)
  vals.names <- c()
  vals.genes <- list()
  for(i in seq_along(vals.range))
  {
    if(i == 1)
    {
      vals.names <- c(vals.names, paste('MeanBeta < ', vals.range[i]))
      vals.genes[[i]] <- names(meth.mean.diff)[which(meth.mean.diff < vals.range[i])]
    }
    else
    {
      vals.names <- c(vals.names, paste(vals.range[i-1], '< MeanBeta <', vals.range[i]))
      vals.genes[[i]] <- names(meth.mean.diff)[which(meth.mean.diff >= vals.range[i-1] & meth.mean.diff < vals.range[i])]
    }
  }
  vals.names <- c(vals.names, paste('MeanBeta >', vals.range[i]))
  vals.genes[[i+1]] <- names(meth.mean.diff)[which(meth.mean.diff >= vals.range[i])]
  
  ob.list$meth.matched.means <- lapply(ob.list$meth.matched.means, function(x) x[unlist(vals.genes)])
  df.req <- data.frame(genes = unlist(vals.genes), log2FC = ob.list$rna.log.fc[ens.genes.map[unlist(vals.genes)]], 
            type = unlist(mapply(function(x,y) rep(x,y), vals.names, sapply(vals.genes,length))),  
            meth.mean.tum = ob.list$meth.matched.means$mean.tum, meth.mean.norm = ob.list$meth.matched.means$mean.norm,
            meth.mean.diff = ob.list$meth.matched.means$mean.diff, row.names = unlist(vals.genes))
  df.req$type <- factor(df.req$type, levels = vals.names)
  if(!is.null(ob.list$rna.matched.means))
  {
    ob.list$rna.matched.means <- lapply(ob.list$rna.matched.means, function(x) x[ens.genes.map[unlist(vals.genes)]])
    df.req <- cbind(df.req, rna.mean.tum = ob.list$rna.matched.means$mean.tum, rna.mean.norm = ob.list$rna.matched.means$mean.norm)
  }
  return(df.req)
  #vals.table <- vector(mode = 'list', length = length(vals.range))
  #for(i in seq_along(vals.range))
}

create.paper.plot <- function(df.req, log2FC = T, mean.norm = F, sc.plot = T)
{
  library(ggplot2)
  library(ggpubr)
  if(log2FC)
    p.bp <- ggplot(df.req, aes(x=type, y=log2FC))
  else
    p.bp <- ggplot(df.req, aes(x=type, y=rna.mean.tum))
  
  p.bp <- p.bp + geom_boxplot() + theme_classic()

  if(mean.norm)
  {
    p.norm <- ggplot(df.req, aes(x=type, y=rna.mean.norm))
    p.norm <- p.norm + geom_boxplot() + theme_classic()
    p.bp <- ggarrange(p.norm, p.bp, ncol = 2)
  }
  p.sc <- NULL
  if(sc.plot)
  {
    p.sc <- ggplot(df.req, aes(x=log2FC, y = meth.mean.diff)) +
      geom_point(aes(color = type)) + stat_smooth(method = 'lm', se = F)
  }
  return(list(p.bp, p.sc))
}

get.type <- function(samples)
{
  types <- sapply(strsplit(samples, '-', fixed = T), function(x) substr(x[4],1,2))
  types.req <- types
  types.req[types == '11'] <- 'N'
  types.req[types != '11'] <- 'T'
  return(types.req)
}