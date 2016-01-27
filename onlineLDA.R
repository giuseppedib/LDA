# Stochastic Variational EM

library(foreach)
library(iterators)
library(doParallel)

source("M_step.R")
source("E_step.R")
source("helpers.R")

# compute perplexity
calculate_perplexity = function(test_W, alpha, beta, max_iter=30, convergence_threshold=1e-6){
  n_topics = length(alpha)
  M = dim(test_W)[1]
  V = dim(test_W)[2]
  
  gamma <- matrix((alpha + V/n_topics), M, n_topics)
  phi <- array(0, dim = c(V, n_topics, M))
  
  res_obj = E_step(gamma, phi, alpha, beta, test_W, max_iter, convergence_threshold)
  
  loglikelihood = res_obj$likelihood
  perplexity = exp(- loglikelihood / sum(test_W))
  return(perplexity)
}



onlineLDA <- function(w, testw, n_topics, n_iter = 100, conv_threshold = 1e-3, MAX_ITER_var_EM = 1000) {
  k = n_topics
  V = ncol(w)
  M = nrow(w)
  eta <- 1
  alpha <- rep(1, k)
  gamma <- matrix((alpha + V/k), M, k)
  phi <- array(0, dim = c(V, k, M))
  lambda <- matrix(rgamma(k*V, 100, 100), nrow=k, ncol=V)
  
  M <- nrow(gamma)
  conv_likelihood <- 100
  perplexity = rep(NA, n_iter %/% 10)
  for(i in 1:n_iter){
    rho <- (64 + i)**(-0.75)
    doc <- sample(c(1:M), 1)
    
    update_obj = online_update_single_doc(alpha, phi[,,doc], gamma[doc,], w[doc, ], lambda, conv_threshold, MAX_ITER_var_EM)
    phi[, , doc] = update_obj$Phi
    gamma[doc, ] = update_obj$gamma
    
    mean_lambda <- eta + M * t(phi[, , doc])
    lambda <- (1 - rho) * lambda + rho * mean_lambda
    
    # update alpha and eta
    newalpha <- update_alpha(alpha, gamma[doc, , drop=FALSE], 1)
    alpha <- (1 - rho) * alpha + rho * newalpha
    cat("Finished iter", i, "\n")
    cat("alpha", alpha[1], "\n\n")
    
    beta <- lambda / rowSums(lambda)
    if(i %% 100 == 0){
      j = i %/% 100
      perplexity[j] = calculate_perplexity(testw, alpha, beta)
    }
  }
  return(list(alpha=alpha, eta=eta, gamma = gamma, lambda = lambda, beta=beta, perplexity = perplexity))
}

online_update_single_doc = function(alpha, Phi, gamma, w_doc, lambda, conv_threshold, MAX_ITER_var_EM){
  iter2 <- 0
  conv_value <- 100
  while ( conv_value > conv_threshold && iter2 < MAX_ITER_var_EM) {
    old_phi <- Phi
    old_gamma <- gamma
    tmp_phi <- t(exp( digamma(gamma)- digamma(sum(gamma)) + digamma(lambda) - digamma(rowSums(lambda)) ))
    Phi <- tmp_phi / rowSums(tmp_phi) * w_doc
    gamma <- alpha + colSums(Phi)
    conv_value <- mean(abs(old_gamma - gamma))
    iter2 <- iter2 + 1
  }
  return(list(Phi = Phi, gamma=gamma))
}

onlineLDA_minibatch <- function(w, testw, n_topics, batch_size = 8, n_iter = 100, conv_threshold = 1e-3, MAX_ITER_var_EM = 1000) {
  k = n_topics
  V = ncol(w)
  M = nrow(w)
  eta <- 1
  alpha <- rep(1, k)
  gamma <- matrix((alpha + V/k), M, k)
  phi <- array(0, dim = c(V, k, M))
  lambda <- matrix(rgamma(k*V, 100, 100), nrow=k, ncol=V)
  
  `%op%` <- if (getDoParRegistered()) `%dopar%` else `%do%`
  
  M <- nrow(gamma)
  conv_likelihood <- 100
  perplexity = rep(NA, n_iter %/% 100)
  for(i in 1:n_iter){
    rho <- (10 + i)**(-0.75)
    doc <- sample(c(1:M), batch_size)
    
    input_list = create_iter_list(gamma[doc, ], phi[, , doc], w[doc, ])
    res = foreach(input = iter(input_list), .export = c("online_update_single_doc")) %op% {
      online_update_single_doc(alpha, input$phi, input$gamma, input$W, lambda, conv_threshold, MAX_ITER_var_EM)
    }
    
    # update_obj = online_update_single_doc(alpha, phi[,,doc], gamma[doc,], w[doc, ], lambda, conv_threshold, MAX_ITER_var_EM)
    new_lambda = matrix(eta, nrow(lambda), ncol(lambda))
    for(k in 1:length(doc)){
      j = doc[k]
      phi[, , j] = res[[k]]$Phi
      gamma[j, ] = res[[k]]$gamma
      new_lambda = new_lambda + t(phi[, , j])
    }
    
    # mean_lambda <- eta + M * t(phi[, , doc])
    lambda <- (1 - rho) * lambda + rho * (M / batch_size)* new_lambda
    
    # update alpha and eta
    newalpha <- update_alpha(alpha, gamma[doc, , drop=FALSE], batch_size)
    alpha <- (1 - rho) * alpha + rho * newalpha
    cat("Finished iter", i, "\n")
    cat("alpha", alpha[1], "\n\n")
    
    beta <- lambda / rowSums(lambda)
    if(i %% 100 == 0){
      j = i %/% 100
      perplexity[j] = calculate_perplexity(testw, alpha, beta)
    }
  }
  return(list(alpha=alpha, eta=eta, gamma = gamma, lambda = lambda, beta=beta, perplexity = perplexity))
}

online_update_single_doc = function(alpha, Phi, gamma, w_doc, lambda, conv_threshold, MAX_ITER_var_EM){
  iter2 <- 0
  conv_value <- 100
  
  while ( conv_value > conv_threshold && iter2 < MAX_ITER_var_EM) {
    old_phi <- Phi
    old_gamma <- gamma
    tmp_phi <- t(exp( digamma(gamma)- digamma(sum(gamma)) + digamma(lambda) - digamma(rowSums(lambda)) ))
    Phi <- tmp_phi / rowSums(tmp_phi) * w_doc
    gamma <- alpha + colSums(Phi)
    conv_value <- mean(abs(old_gamma - gamma))
    iter2 <- iter2 + 1
  }
  return(list(Phi = Phi, gamma=gamma))
}


# compute_likelihood = function(gamma, Phi, alpha, lambda, eta){
#   E_log_theta = digamma(gamma) - digamma(sum(gamma))
#   
#   E_log_beta = digamma(lambda) - digamma(rowSums(lambda))
#   
#   V = ncol(lambda)
#   out = sum(Phi * (t(E_log_theta + E_log_beta) - ifelse(Phi>0, log(Phi), 0))) + 
#     -lgamma(sum(gamma)) + sum(lgamma(gamma)) + sum((alpha - gamma)*E_log_theta) +
#     sum(-lgamma(rowSums(lambda)) + rowSums((eta - lambda)*E_log_beta + lgamma(lambda))) +
#     lgamma(sum(alpha)) - sum(lgamma(alpha)) + (lgamma(V*eta) - V*lgamma(eta))
#   return(out)
# }
