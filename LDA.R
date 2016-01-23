# V <- 10              # number of words in the vocabulary
# M <- 4               # number of documents in the corpus
# k <- 5               # number of topics
# w <- matrix(0,M,V)   # matrix of counts
# 

library(gtools)

source("E_step.R")
source("M_step.R")
source("helpers.R")


LDA = function(W, n_topics = 3, max_iter = 100, max_iter_E_step = 10, convergence_threshold=1e-6){
  k = n_topics
  V = ncol(W)
  M = nrow(W)
  alpha <- rep(1, k)
  beta <- rdirichlet(k, rep(1, V))
  gamma <- matrix((alpha + V/k), M, k) #lapply(1:M, function(y) alpha + V/k)
  phi <- array(0, dim = c(V, k, M))
  
  likelihood = rep(NA, max_iter)
  for(i in 1:max_iter){
    # E step
    obj = E_step(gamma, phi, alpha, beta, W, max_iter_E_step, convergence_threshold)
    likelihood[i] = obj$likelihood
    phi = obj$phi
    gamma = obj$gamma
    # M step
    beta = update_beta(phi, W, k)
    alpha = update_alpha(alpha, gamma, M)
    
    if(check_convergence(likelihood, i, convergence_threshold)) break
  }
  return(list(likelihood = likelihood[1:i], phi = obj$phi, gamma = obj$gamma, alpha = alpha, beta=beta))
}

