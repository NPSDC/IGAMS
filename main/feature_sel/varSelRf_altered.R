library(randomForest)
varSelRF <- function(xdata, Class,
                     c.sd = 1,
                     mtryFactor = 1,
                     ntree = 5000,
                     ntreeIterat = 2000,                     
                     vars.drop.num = NULL,
                     vars.drop.frac = 0.2,
                     whole.range = TRUE,
                     recompute.var.imp = FALSE,
                     verbose = FALSE,
                     returnFirstForest = TRUE,
                     fitted.rf = NULL,
                     keep.forest = FALSE,
                     strata=NULL, sampsize = if (replace) nrow(xdata) else ceiling(.632*nrow(x)),
                     classwt = NULL, nodesize = 1, replace = T) {
  
  if(!is.factor(Class))
    stop("Class should be a factor")
  if( (is.null(vars.drop.num) & is.null(vars.drop.frac)) |
      (!is.null(vars.drop.num) & !is.null(vars.drop.frac)))
    stop("One (and only one) of vars.drop.frac and vars.drop.num must be NULL and the other set")
  
  max.num.steps <- dim(xdata)[2]
  num.subjects <- dim(xdata)[1]
  
  if(is.null(colnames(xdata)))
    colnames(xdata) <- paste("v", 1:dim(xdata)[2], sep ="")
  
  ##oversize the vectors; will prune later.
  n.vars <- vars <- OOB.rf <- OOB.sd <- mcc.rf <- rep(NA, max.num.steps)
  
  oobError <- function(rf) { ## should not be exported in the namespace.
    ## The out of bag error
    ooo <- rf$confusion[, -dim(rf$confusion)[2]]
    s.ooo <- sum(ooo)
    diag(ooo) <- 0
    sum(ooo)/s.ooo
  }
  
  compute.mcc <- function(conf.mat)
  {
    TP=conf.mat[1,1]
    TN=conf.mat[2,2]
    FP=conf.mat[2,1]
    FN=conf.mat[1,2]
    return((TP*TN-FP*FN)/sqrt(((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  }
  
  if(!is.null(fitted.rf)) {
    if(ncol(fitted.rf$importance) < 2)
      stop("The fitted rf was not fitted with importance = TRUE")
    n.ntree <- fitted.rf$ntree
    mtry <- fitted.rf$mtry
    n.mtryFactor <- mtry/sqrt(ncol(xdata))
    if((n.ntree != ntree) | (n.mtryFactor != mtryFactor))
      warning("Using as ntree and mtry the parameters obtained from fitted.rf",
              immediate.= TRUE)
    ntree <- n.ntree
    mtryFactor <- n.mtryFactor
    rm(n.ntree, n.mtryFactor)
    rf <- fitted.rf
  } else {
    mtry <- floor(sqrt(ncol(xdata)) * mtryFactor)
    rf <- randomForest(x = xdata, y = Class,
                       ntree = ntree, mtry = mtry,
                       importance = TRUE,
                       keep.forest = keep.forest,
                       classwt = NULL, sampsize = sampsize, strata = strata, nodesize = nodesize, replace = replace)
  }
  
  if(returnFirstForest)
    FirstForest <- rf
  else
    FirstForest <- NULL
  m.iterated.ob.error <- m.initial.ob.error <- oobError(rf)
  m.iterated.mcc <- m.initial.mcc <- compute.mcc(rf$confusion)
  sd.iterated.ob.error <- sd.initial.ob.error <-
    sqrt(m.iterated.ob.error * (1 - m.iterated.ob.error) *
           (1/num.subjects))
  
  if(verbose) {
    print(paste("Initial MCC error: mean = ",
                round(m.initial.mcc, 4),
                "; sd = ", round(sd.initial.ob.error, 4), sep = ""))
  }
  
  importances <- importance(rf, type = 1, scale = FALSE)
  selected.vars <- order(importances, decreasing = TRUE)
  ordered.importances <- importances[selected.vars]
  
  initialImportances <- importances
  initialOrderedImportances <- ordered.importances
  
  j <- 1
  n.vars[j] <- dim(xdata)[2] 
  vars[j] <- paste(colnames(xdata), collapse = " + ")
  OOB.rf[j] <- m.iterated.ob.error
  OOB.sd[j] <- sd.iterated.ob.error
  mcc.rf[j] <- m.iterated.mcc
  var.simplify <- TRUE
  
  while(var.simplify) {
    if (verbose){
      print("gc inside loop of varSelRF")
      print(gc())
    } else {
      gc()
    }
    
    last.rf <- rf 
    last.vars <- selected.vars
    previous.m.error <- m.iterated.ob.error
    previous.sd.error <- sd.iterated.ob.error
    
    ## If this is left as only
    ## "if(length(selected.vars) <= 2) var.simplify <- FALSE"
    ## under the if((length(selected.vars) < 2) | (any(selected.vars < 1))),
    ## as it used to be, then we fit a 2 model var, which might be
    ## better or as good as others, but as we never re-enter,
    ## we cannot return it, even if we see it in the history.
    
    ## Alternatively, we cannot just convert
    ## "if((length(selected.vars) < 2) | (any(selected.vars < 1))"
    ## to the <= 2, as we then would re-enter many times because
    ## of the way selected.vars <- selected.vars[1:2] when
    ## num.vars < (vars.drop + 2)
    
    ## This way, we enter just to set last.rf, last.vars and
    ## we bail out
    
    if(length(selected.vars) <= 2) {
      var.simplify <- FALSE
      break
    }
    
    
    if(recompute.var.imp & (j > 1)) {
      importances <- importance(rf, type = 1, scale = FALSE)
      tmp.order <- order(importances, decreasing = TRUE)
      selected.vars <- selected.vars[tmp.order]
      ordered.importances <- importances[tmp.order]
    }
    
    num.vars <- length(selected.vars)
    
    if(is.null(vars.drop.num))
      vars.drop <- round(num.vars * vars.drop.frac)
    else vars.drop <- vars.drop.num
    
    if(num.vars >= (vars.drop + 2)) {
      selected.vars <- selected.vars[1: (num.vars - vars.drop)]
      ordered.importances <- ordered.importances[1: (num.vars - vars.drop)]
    }
    else {
      selected.vars <- selected.vars[1:2]
      ordered.importances <- ordered.importances[1:2]
    }
    
    ## couldn't we eliminate the following?
    if((length(selected.vars) < 2) | (any(selected.vars < 1))) {
      var.simplify <- FALSE
      break
    }
    
    
    
    mtry <- floor(sqrt(length(selected.vars)) * mtryFactor)
    if(mtry > length(selected.vars)) mtry <- length(selected.vars)
    
    if(recompute.var.imp) 
      rf <- randomForest(x = xdata[, selected.vars],
                         y = Class, importance= TRUE,
                         ntree = ntree, mtry = mtry,
                         keep.forest = keep.forest,
                         classwt = NULL, sampsize = sampsize, strata = strata, nodesize = nodesize, replace = replace)
    else
      rf <- randomForest(x = xdata[, selected.vars],
                         y = Class, importance= FALSE,
                         ntree = ntreeIterat, mtry = mtry,
                         keep.forest = keep.forest,
                         classwt = NULL, sampsize = sampsize, strata = strata, nodesize = nodesize, replace = replace)
    
    m.iterated.ob.error <- oobError(rf)
    m.iterated.mcc <- compute.mcc(rf$confusion)
    sd.iterated.ob.error <-
      sqrt(m.iterated.ob.error * (1 - m.iterated.ob.error) *
             (1/num.subjects))
    
    if(verbose) {
      print(paste("..... iteration ", j, "; OOB error: mean = ",
                  round(m.iterated.ob.error, 4),
                  "; sd = ", round(sd.iterated.ob.error, 4),
                  "; num. vars = ", length(selected.vars), 
                  sep = ""))
    }
    j <- j + 1
    
    
    n.vars[j] <- length(selected.vars)
    vars[j] <- paste(colnames(xdata)[selected.vars],
                     collapse = " + ")
    OOB.rf[j] <- m.iterated.ob.error
    OOB.sd[j] <- sd.iterated.ob.error
    mcc.rf[j] <- m.iterated.mcc
    
    if(!whole.range &
       (
         (m.iterated.ob.error >
          (m.initial.ob.error + c.sd*sd.initial.ob.error))
         |
         (m.iterated.ob.error >
          (previous.m.error + c.sd*previous.sd.error)))
    )
      var.simplify <- FALSE
  }
  
  if (!whole.range) {
    if(!is.null(colnames(xdata)))
      selected.vars <- sort(colnames(xdata)[last.vars])
    else
      selected.vars <- last.vars
    
    out <- list(selec.history = data.frame(
      Number.Variables = n.vars,
      Vars.in.Forest = vars,
      OOB = OOB.rf,
      sd.OOB = OOB.sd, MCC = mcc.rf)[1:j,],
      rf.model = last.rf,
      selected.vars = selected.vars,
      selected.model =  paste(selected.vars, collapse = " + "),
      best.model.nvars = length(selected.vars),
      initialImportances = initialImportances,
      initialOrderedImportances = initialOrderedImportances,
      ntree = ntree,
      ntreeIterat = ntreeIterat,
      mtryFactor = mtryFactor,
      #                    mtry = mtry,
      firstForest = FirstForest)
    class(out) <- "varSelRF"
    return(out)
  }
  else {
    ## Prune the too long vectors created at begin.
    ## not needed above, because we select the 1:j rows
    ## of the return matrix selec.history.
    n.vars <- n.vars[1:j]
    vars <- vars[1:j]
    OOB.rf<- OOB.rf[1:j]
    OOB.sd <- OOB.sd[1:j]
    mcc.rf <- mcc.rf[1:j]
    min.oob.ci <- min(OOB.rf) #+ c.sd * OOB.sd[which.min(OOB.rf)]
    best.pos <-
      which(OOB.rf == min.oob.ci)[which.min(n.vars[which(OOB.rf == min.oob.ci)])]
    
    selected.vars <- sort(unlist(strsplit(vars[best.pos],
                                          " + ", fixed = TRUE)))
    out <- list(selec.history = data.frame(
      Number.Variables = n.vars,
      Vars.in.Forest = vars,
      OOB = OOB.rf,
      sd.OOB = OOB.sd, MCC = mcc.rf),
      rf.model = NA,
      selected.vars = selected.vars,
      selected.model = paste(selected.vars, collapse = " + "),
      best.model.nvars = n.vars[best.pos],
      initialImportances = initialImportances,
      initialOrderedImportances = initialOrderedImportances,
      ntree = ntree,
      ntreeIterat = ntreeIterat,
      mtryFactor = mtryFactor,
      ##mtry = mtry,
      firstForest = FirstForest)
    class(out) <- "varSelRF"
    return(out)
  }
}
