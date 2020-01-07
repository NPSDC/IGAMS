source("main/feature_sel/get_features.R")
source("main/build_models.R")
source("integration/load_basic.R")
source("integration/kernel_init.R")

load("environment/rnaseq/vst_rna_meth_common.RData")
load("environment/rnaseq/rna_fea.RData")
load("environment/methylation/comb_stages.RData") ##Stage information for the samples
load("environment/methylation/train_common_ind.RData") ##Training index within 250 samples
load("environment/methylation/test_common_ind.RData") ##Test index within the 250 samples


##Training methylation model
pap.rna.train.model <- get.train.model(tr.data = vst.rna.req[train.common.index,], fea.list = rna.fea, 
                                        stages.train = comb.stage[train.common.index], C_set = C_set, n.trees = n.trees, cores = 4)
save(pap.rna.train.model, file = "environment/rnaseq/pap_rna_train_model.RData")

##Cross Validation
pap.rna.cv.model <- get.cv.model(tr.data = vst.rna.req[train.common.index,], 
                                  fea.list = rna.fea, folds = 5,
                                  stages.train = comb.stage[train.common.index],
                                  tr.model = pap.rna.train.model, cores = 4)
save(pap.rna.cv.model, file = "environment/rnaseq/pap_rna_cv_model.RData")

##Predicting test stages
pap.rna.test.pred <- get.test.pred(tr.data = vst.rna.req[train.common.index,], 
                                    te.data = vst.rna.req[test.common.index,],
                                    fea.list = rna.fea, stages.tr = comb.stage[train.common.index],
                                    tr.model = pap.rna.train.model, cv.model = pap.rna.cv.model, cores = 4)
save(pap.rna.test.pred, file = "environment/rnaseq/pap_rna_test_pred.RData")

pap.rna.test.res <- get.test.res(stages.test = comb.stage[test.common.index], test.pred = pap.rna.test.pred)
save(pap.rna.test.res, file = "environment/rnaseq/pap_rna_test_res.RData")

##MKL
gr <- build.groups(stages = comb.stage[train.common.index], num.group = 5, strat = T) 

##Group Lasso
tr.rna.gr.lasso <- list()
tr.rna.gr.lasso[["varSelRF"]] <- list()

tr.rna.gr.lasso[["varSelRF"]][["atleast_1"]] <- train.group.lasso.trial(train.data.list = list(vst.rna.req[train.common.index,]), 
                                                                    test.data.list = list(vst.rna.req[test.common.index,]), 
                                                                    stages.train = comb.stage[train.common.index], stages.test = comb.stage[test.common.index], 
                                                                    features.list = list(list("f" = rna.fea$varSelRF$atleast_1)),
                                                                    cv.costs.res = NULL, costs = C_set, kern.func.list = kernels, 
                                                                    kern.hash = kernels.hash.table, kern.ind.list = k.gr, params = parameters, folds = 5, 
                                                                    gr.list = list(gr), n = 1, unit.norm = T, cores = 4)

save(tr.rna.gr.lasso, file = "environment/rnaseq/tr_rna_gr_lasso.RData")

##BEMKL
tr.rna.bemkl <- list()
tr.rna.bemkl[["varSelRF"]] <- list()

tr.rna.bemkl[["varSelRF"]][["atleast_1"]] <- train.bemkl.trial(train.data.list = list(vst.rna.req[train.common.index,]), 
                                                           test.data.list = list(vst.rna.req[test.common.index,]), 
                                                           stages.train = comb.stage[train.common.index], stages.test = comb.stage[test.common.index], 
                                                           features.list = list(list("f" = rna.fea$varSelRF$atleast_1)),
                                                           cv.costs.res = NULL, gamma.prior = list(c(1,1), c(1e-10, 1e+10)), kern.func.list = kernels,
                                                           kern.hash = kernels.hash.table, kern.ind.list = k.gr, params = params, folds = 5, 
                                                           gr.list = list(gr), n = 1, unit.norm = T, cores = 4)
save(tr.rna.bemkl, file = "environment/rnaseq/tr_rna_bemkl.RData")

