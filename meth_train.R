source("main/feature_sel/get_features.R")
source("main/build_models.R")
source("integration/load_basic.R")
source("integration/kernel_init.R")

setwd("../Stage-Prediction-of-Cancer/Organised/")
load("../methylation/environment/prcc/req_meth.RData") ##data.frame of methylation CpGs probess for the 250 samples common to RNASeq and methylation
load("../methylation/environment/prcc/comb_stages.RData") ##Stage information for the samples
load("../methylation/environment/prcc/train_common_ind.RData") ##Training index within 250 samples
load("../methylation/environment/prcc/test_common_ind.RData") ##Test index within the 250 samples
setwd("../../IGAMS/")

##Extracting methylation features
meth.req.data <- meth.req.data[1:100,]
meth.fea.ob <- get.meth.fea(train.data = t(logit2(meth.req.data[, train.common.index])), train.stages = comb.stage[train.common.index],
             train.ind.list = list(seq(length(train.common.index))), cores = 4,
             type = 2, min.genes = 2, min.range = 0.02)
meth.fea.list <- get.class.fea(meth.fea.ob)
save(meth.fea.ob, file = "environment/methylation/meth_fea_ob.RData")
save(meth.fea.list, file = "environment/methylation/meth_fea_list.RData")

##Training methylation model
pap.meth.train.model <- get.train.model(tr.data = t(logit2(meth.req.data[,train.common.index])), fea.list = meth.fea.list, 
                                        stages.train = comb.stage[train.common.index], C_set = C_set, n.trees = n.trees, cores = 4)
save(pap.meth.train.model, file = "environment/methylation/pap_meth_train_model.RData")

##Cross Validation
pap.meth.cv.model <- get.cv.model(tr.data = t(logit2(meth.req.data[,train.common.index])), 
                               fea.list = meth.fea.list, folds = 5,
                               stages.train = comb.stage[train.common.index],
                               tr.model = pap.meth.train.model, cores = 4)
save(pap.meth.cv.model, file = "environment/methylation/pap_meth_cv_model.RData")

##Predicting test stages
pap.meth.test.pred <- get.test.pred(tr.data = t(logit2(meth.req.data[,train.common.index])), 
                                    te.data = t(logit2(meth.req.data[,test.common.index])),
                                    fea.list = meth.fea.list, stages.tr = comb.stage[train.common.index],
                                    tr.model = pap.meth.train.model, cv.model = pap.meth.cv.model, cores = 4)
save(pap.meth.test.pred, file = "environment/methylation/pap_meth_test_pred.RData")

##Evaluating models on test data
pap.meth.test.res <- get.test.res(stages.test = comb.stage[test.common.index], test.pred = pap.meth.test.pred)
save(pap.meth.test.res, file = "environment/methylation/pap_meth_test_res.RData")

##MKL
gr <- build.groups(stages = comb.stage[train.common.index], num.group = 5, strat = T) 

##Group Lasso
tr.gr.lasso <- list()
tr.gr.lasso[["shrunken"]] <- list()
tr.gr.lasso[["varSelRF"]] <- list()

tr.gr.lasso[["varSelRF"]][["atleast_1"]] <- train.group.lasso.trial(train.data.list = list(t(logit2(meth.req.data[,train.common.index]))), 
                                            test.data.list = list(t(logit2(meth.req.data[,test.common.index]))), 
                                            stages.train = comb.stage[train.common.index], stages.test = comb.stage[test.common.index], 
                                            features.list = list(list("f" = meth.fea.list$varSelRF$atleast_1)),
                                            cv.costs.res = NULL, costs = C_set, kern.func.list = kernels, 
                                            kern.hash = kernels.hash.table, kern.ind.list = k.gr, params = parameters, folds = 5, 
                                            gr.list = list(gr), n = 1, unit.norm = T, cores = 4)

tr.gr.lasso[["shrunken"]][["atleast_1"]] <- train.group.lasso.trial(train.data.list = list(t(logit2(meth.req.data[,train.common.index]))),
                                            test.data.list = list(t(logit2(meth.req.data[,test.common.index]))),
                                            stages.train = comb.stage[train.common.index], stages.test = comb.stage[test.common.index],
                                            features.list = list(list("f" = meth.fea.list$shrunken$atleast_1)),
                                            cv.costs.res = NULL, costs = C_set, kern.func.list = kernels,
                                            kern.hash = kernels.hash.table, kern.ind.list = k.gr, params = parameters, folds = 5,
                                            gr.list = list(gr), n = 1, unit.norm = T, cores = 4)
save(tr.gr.lasso, file = "environment/methylation/tr_gr_lasso.RData")

##BEMKL
tr.bemkl <- list()
tr.bemkl[["shrunken"]] <- list()
tr.bemkl[["varSelRF"]] <- list()

tr.bemkl[["varSelRF"]][["atleast_1"]] <- train.bemkl.trial(train.data.list = list(t(logit2(meth.req.data[,train.common.index]))), 
                                            test.data.list = list(t(logit2(meth.req.data[,test.common.index]))), 
                                            stages.train = comb.stage[train.common.index], stages.test = comb.stage[test.common.index], 
                                            features.list = list(list("f" = meth.fea.list$varSelRF$atleast_1)),
                                            cv.costs.res = NULL, gamma.prior = list(c(1,1), c(1e-10, 1e+10)), kern.func.list = kernels,
                                            kern.hash = kernels.hash.table, kern.ind.list = k.gr, params = params, folds = 5, 
                                            gr.list = list(gr), n = 1, unit.norm = T, cores = 4)
tr.bemkl[["shrunken"]][["atleast_1"]] <- train.bemkl.trial(train.data.list = list(t(logit2(meth.req.data[,train.common.index]))), 
                                                test.data.list = list(t(logit2(meth.req.data[,test.common.index]))),
                                                stages.train = comb.stage[train.common.index], stages.test = comb.stage[test.common.index],
                                                features.list = list(list("f" = meth.fea.list$shrunken$atleast_1)),
                                                cv.costs.res = NULL, gamma.prior = list(c(1,1), c(1e-10, 1e+10)), kern.func.list = kernels,
                                                kern.hash = kernels.hash.table, kern.ind.list = k.gr, params = params, folds = 5,
                                                gr.list = list(gr), n = 1, unit.norm = T, cores = 4)
save(tr.bemkl, file = "environment/methylation/tr_bemkl.RData")
