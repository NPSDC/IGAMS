##Setting the kernel paramaters

##RBF kernels
#rbf.kernels.params <- c(5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1, 1,10, 100, 1000, 1e4) #12 sigma
rbf.kernels.params <- c(5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1, 1,2,5,8, 10,16,32,50,64,75,100,128,256,512,1000, 1e4) #12 sigma
names(rbf.kernels.params) <- paste0('rbf_sigma_', as.character(rbf.kernels.params))

##Polynomial kernels
poly.kernels.params <- c(1,2,3,4) # 4 degrees
names(poly.kernels.params) <- paste0('poly_degree_', as.character(poly.kernels.params))

##Kernels we plan on using
kernels <- list()
kernels[names(rbf.kernels.params)] <- lapply(rbf.kernels.params, function(si){
                create.kernel.function(kernel.name = 'rbf', kernel.param = c('sigma'=si))
  })
kernels[names(poly.kernels.params)] <- lapply(poly.kernels.params, function(degree){
  create.kernel.function(kernel.name = 'polynomial', kernel.param = c('degree' = degree))
})
#kernels[['gower']] <- create.kernel.function(kernel.name = 'gower')
#kernels[['jaccard']] <- create.kernel.function(kernel.name = 'jaccard')
remove(rbf.kernels.params)
remove(poly.kernels.params)
##Maintaining hash table
kernels.hash.table <- seq(length(kernels))
names(kernels.hash.table) <- names(kernels)

#List maintaining indexes of kernels to be applied on each parition of training set
#kernels.ind.list <- list()
#kernels.ind.list[['all_genes']] <- seq(length(kernels))

C_set <- c(0.0001, 0.001, 0.01, 0.1, 1, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5)
n.trees <- c(100,250,500,1000,1500,2000,2500)
parameters <- list(epsilon = 1e-3, iteration_count = 200, C = 1)

#k.gr <- list(c(1:10, 13:16), c(2:10, 13:15), c(1:16), c(2:16), c(3:16), c(4:16), c(2:11, 13:15), 
#              c(2:12, 13:15))

k.gr <- list(c(1:24), c(2:24), c(4:24), c(1:18, 21:24), c(2:18, 21:23), c(2:19, 21:23), c(2:20, 21:23))

##Bemkl
params <- list()
params$alpha_lambda <- 1
params$beta_lambda <- 1

#set the hyperparameters of gamma prior used for bias
params$alpha_gamma <- 1
params$beta_gamma <- 1

#set the hyperparameters of gamma prior used for kernel weights
params$alpha_omega <- 1
params$beta_omega <- 1

params$iteration <- 200

#set the margin parameter
params$margin <- 1

#determine whether you want to store the lower bound values
params$progress <- 0

#set the seed for random number generator used to initalize random variables
params$seed <- 1606

#set the standard deviation of intermediate representations
params$sigma_g <- 0.1
