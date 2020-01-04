source("integration/integrate.R")
source("integration/load_basic.R")
library(minfi)

setwd("../Stage-Prediction-of-Cancer/Organised/")
load("../methylation/environment/prcc/req_meth.RData") ##data.frame of methylation CpGs probess for the 250 samples common to RNASeq and methylation
load("../methylation/environment/prcc/comb_stages.RData") ##Stage information for the samples
load("../methylation/environment/prcc/train_common_ind.RData") ##Training index within 250 samples
load("../methylation/environment/prcc/test_common_ind.RData") ##Test index within the 250 samples
load("../methylation/environment/prcc/meth_fea_list.RData") ##Methylation features
load("environment/prcc/vst_rna_meth_common.RData") ##VST normalised RNASeq Data
load("environment/prcc/rna_fea.RData") ##RNASeq features
setwd("../../IGAMS/")

meth.req.data <- meth.req.data[unique(unlist(meth.fea.list)),]
vst.rna.req <- vst.rna.req[,unique(unlist(rna.fea))]
gc()
int.res <- simple.integrate.t(train.data.list = list(vst.rna.req[train.common.index,], t(logit2(meth.req.data[,train.common.index]))),
                           test.data.list = list(vst.rna.req[test.common.index,], t(logit2(meth.req.data[,test.common.index]))),
                           fea.list = list(list("varGene" = rna.fea$varSelRF$atleast_1), 
                                           list("varMeth" = meth.fea.list$varSelRF$atleast_1)),
                           stages = list('train' = comb.stage[train.common.index], 'test' = comb.stage[test.common.index]),
                           n.trees = n.trees, C_set = C_set, scale = c(T,T), cores = 2, cv = F)

registerDoSEQ()
save(int.res, file = "environment/integration/int_res.RData")
gr <- build.groups(stages = comb.stage[train.common.index], num.group = 5, strat = T)

tr.rna.meth.gr <- list()
tr.rna.meth.gr[["shrunken"]] <- list()
tr.rna.meth.gr[["varSelRF"]] <- list()

tr.rna.meth.gr[["varSelRF"]][["atleast_1"]] <- train.group.lasso.trial(train.data.list = list(vst.rna.req[train.common.index,], 
                                                                                                   t(logit2(meth.req.data[,train.common.index]))), 
                                                                            test.data.list = list(vst.rna.req[test.common.index,],
                                                                                                  t(logit2(meth.req.data[,test.common.index]))), 
                                                                            stages.train = comb.stage[train.common.index], stages.test = comb.stage[test.common.index], 
                                                                            features.list = list(list("f" = rna.fea$varSelRF$atleast_1), 
                                                                                                 list("v" = meth.fea.list$varSelRF$atleast_1)),
                                                                            cv.costs.res = NULL, costs = C_set, kern.func.list = kernels, 
                                                                            kern.hash = kernels.hash.table, kern.ind.list = k.gr, params = parameters,
                                                                            folds = 5, gr.list = list(gr), n = 1, unit.norm = T, cores = 2)
tr.rna.meth.gr[["shrunken"]][["atleast_1"]] <- train.group.lasso.trial(train.data.list = list(vst.rna.req[train.common.index,], 
                                                                                                   t(logit2(meth.req.data[,train.common.index]))),
                                                                            test.data.list = list(vst.rna.req[test.common.index,],
                                                                                                  t(logit2(meth.req.data[,test.common.index]))),
                                                                            stages.train = comb.stage[train.common.index], stages.test = comb.stage[test.common.index],
                                                                            features.list = list(list("f" = rna.fea$varSelRF$atleast_1), 
                                                                                                 list("v" = meth.fea.list$shrunken$atleast_1)),
                                                                            cv.costs.res = NULL, costs = C_set, kern.func.list = kernels,
                                                                            kern.hash = kernels.hash.table, kern.ind.list = k.gr, params = parameters,
                                                                            folds = 5, gr.list = list(gr), n = 1, unit.norm = T, cores = 2)

save(tr.rna.meth.gr, file = "integration/tr_rna_meth_gr.RData")

tr.rna.meth.bemkl <- list()


tr.rna.meth.bemkl[["shrunken"]] <- list()
tr.rna.meth.bemkl[["varSelRF"]] <- list()
    
tr.rna.meth.bemkl[["varSelRF"]][["atleast_1"]] <- train.bemkl.trial(train.data.list = list(vst.rna.req[train.common.index,], 
                                                                                                t(logit2(meth.req.data[,train.common.index]))), 
                                                                         test.data.list = list(vst.rna.req[test.common.index,], 
                                                                                               t(logit2(meth.req.data[,test.common.index]))), 
                                                                         stages.train = comb.stage[train.common.index], stages.test = comb.stage[test.common.index], 
                                                                         features.list = list(list("f" = rna.fea$varSelRF$atleast_1), 
                                                                                              list("v" = meth.fea.list$varSelRF$atleast_1)),
                                                                         cv.costs.res = NULL, gamma.prior = list(c(1,1), c(1e-10, 1e+10)), 
                                                                         kern.func.list = kernels, kern.hash = kernels.hash.table, kern.ind.list = k.gr,
                                                                         params = params, folds = 5, gr.list = list(gr), n = 1, unit.norm = T, cores = 2)
tr.rna.meth.bemkl[["shrunken"]][["atleast_1"]] <- train.bemkl.trial(train.data.list = list(vst.rna.req[train.common.index,],
                                                                                                t(logit2(meth.req.data[,train.common.index]))), 
                                                                         test.data.list = list(vst.rna.req[test.common.index,], 
                                                                                               t(logit2(meth.req.data[,test.common.index]))),
                                                                         stages.train = comb.stage[train.common.index], stages.test = comb.stage[test.common.index], 
                                                                         features.list = list(list("f" = rna.fea$varSelRF$atleast_1), 
                                                                                              list("v" = meth.fea.list$shrunken$atleast_1)),
                                                                         cv.costs.res = NULL, gamma.prior = list(c(1,1), c(1e-10, 1e+10)), kern.func.list = kernels,
                                                                         kern.hash = kernels.hash.table, kern.ind.list = k.gr, params = params, folds = 5,
                                                                         gr.list = list(gr), n = 1, unit.norm = T, cores = 2)

save(tr.rna.meth.bemkl, file = "environment/integration/tr_rna_meth_bemkl.RData")
