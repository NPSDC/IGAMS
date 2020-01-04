group_lasso_multiple_kernel_classification_train <- function(Km, y, parameters) {
  P <- dim(Km)[3]
  eta <- rep(1 / P, P)
  #print(eta)
  Keta <- calculate_Keta(Km, eta)
  #print("ss")
  model <- solve_classification_svm(Keta, y, parameters$C, parameters$epsilon)
  
  #model <- ksvm(as.kernelMatrix(Keta), y = y, C = parameters$C, epsilon = parameters$epsilon)
  objectives <- model$objective
  #objectives <- model@obj
  #print(c(length(objectives), model$objective))
  #print(c(length(objectives), model@obj))
  k <- 1

  while(1) {
    #print("ss")
    start_objective <- model$objective
    #start_objective <- model@obj
    start_eta <- eta
    #print(k)
    k = k+1
    for (m in 1:P) {
     # alpha <- rep(0, length(model@fitted[,1]))
    #  alpha[model@alphaindex] <- model@alpha
      eta[m] <- eta[m] * sqrt(t(model$alpha) %*% Km[,,m] %*% model$alpha)
    #  eta[m] <- eta[m] * sqrt(t(alpha) %*% Km[,,m] %*% alpha)
    }
    eta <- eta / sum(eta)
    eta[eta < parameters$epsilon] <- 0
    eta <- eta / sum(eta)
    Keta <- calculate_Keta(Km, eta)
    #model <- ksvm(as.kernelMatrix(Keta), y = y, C = parameters$C, epsilon = parameters$epsilon)
    modelnew <- tryCatch(
      {
        solve_classification_svm(Keta, y, parameters$C, parameters$epsilon) ##remember this would never return NULL, put model in front of trycatch
      }, error = function(e){ return(NULL)})

    if(is.null(modelnew))
    {
      print("svm cant be solved further")
      break
    }
    model <- modelnew
    #print(paste(k,' k'))
    #objectives <- c(objectives, model@obj)
    objectives <- c(objectives, model$objective)
    #print(c(length(objectives), model$objective))
    if (length(objectives) == parameters$iteration_count) {
      break
    }
  }
  if(is.null(model))
    return(NULL)
  state <- list(alpha = model$alpha, b = model$b, eta = eta, objectives = objectives, parameters = parameters)
}

group_lasso_multiple_kernel_classification_train_new <- function(Km, y, parameters) {
  P <- dim(Km)[3]
  eta <- rep(1 / P, P)
  Keta <- calculate_Keta(Km, eta)
  model <- ksvm(as.kernelMatrix(Keta), y = y, scaled = F, C = parameters$C, epsilon = parameters$epsilon)
  objectives <- model@obj
  
  alpha <- rep(0, length(model@fitted[,1]))
  alpha[model@alphaindex] <- model@alpha
  
  print(c(length(objectives), model@obj))
  k <- 1
  while(1) {
    start_objective <- model@obj
    start_eta <- eta
    #print(k)
    k = k+1

    for (m in 1:P) {
      eta[m] <- eta[m] * sqrt(t(alpha) %*% Km[,,m] %*% alpha)
    }
    #print(length(alpha))
    eta <- eta / sum(eta)
    eta[eta < parameters$epsilon] <- 0
    eta <- eta / sum(eta)
    Keta <- calculate_Keta(Km, eta)
    #print(eta)
    model <- ksvm(as.kernelMatrix(Keta), y = y, scaled = F, C = parameters$C, epsilon = parameters$epsilon)
    alpha[model@alphaindex] <- model@alpha
    objectives <- c(objectives, model@obj)
    #print(c(length(objectives), model$objective))
    if (length(objectives) == parameters$iteration_count) {
      break
    }
  }
  state <- list(alpha = alpha, b = model@b, eta = eta, objectives = objectives, parameters = parameters)
}
