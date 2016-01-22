V <- 10              # number of words in the vocabulary
M <- 4               # number of documents in the corpus
k <- 5               # number of topics
w <- matrix(0,M,V)   # matrix of counts


library(gtools)


source("E_step.R")
source("M_step.r")
source("helpers.R")


LDA = function(W, n_topics = 3, max_iter = 100){
  k = 5
  V = ncol(W)
  M = nrow(W)
  alpha <- rep(10, k)
  beta <- rdirichlet(k, rep(1,V))
  gamma <- matrix((alpha + V/k), M, k) #lapply(1:M, function(y) alpha + V/k)
  phi <- array(1/k, dim = c(V,k,M))
  
  likelihoods = rep(NA, max_iter)
  for(i in 1:max_iter){
    # E step
    obj = E_step(gamma, phi, alpha, beta, W, 10, 1e-3)
    likelihoods[i] = obj$likelihood
    phi = obj$phi
    gamma = obj$gamma
    # M step
    beta = update_beta(phi, W, k)
    alpha = update_alpha(alpha[1], gamma, M)
    
    if(check_convergence(likelihoods, i, convergence_threshold)) break
  }
  return(list(likelihoods = likelihoods[1:i], phi = obj$phi, gamma = obj$gamma, alpha = alpha, beta=beta))
}



LDA(W)
