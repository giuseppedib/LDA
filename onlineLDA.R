# Stochastic Variational EM

source("M_step.R")
source("E_step.R")

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
    doc <- sample(c(1:M),1)
    iter2 <- 0
    conv_value <- 100
    while ( conv_value > conv_threshold && iter2 < MAX_ITER_var_EM) {
      old_phi <- phi[,,doc]
      old_gamma <- gamma[doc,]
      tmp_phi <- t(exp( digamma(gamma[doc,])- digamma(sum(gamma[doc, ])) + digamma(lambda) - digamma(rowSums(lambda)) ))
      phi[ , ,doc] <- tmp_phi / rowSums(tmp_phi) * w[doc, ]
      gamma[doc,] <- alpha + colSums(phi[,,doc])
      conv_value <- mean(abs(old_gamma - gamma[doc,])) # + mean(abs(old_phi - phi[,,doc]))
      iter2 <- iter2 + 1
    }
    cat("Finished iter", i, "This took", iter2, "steps.\n")
    mean_lambda <- eta + M * t(phi[, , doc])
    lambda <- (1 - rho) * lambda + rho * mean_lambda
    
    # update alpha and eta
    newalpha <- update_alpha(alpha, gamma[doc, , drop=FALSE], 1)
    alpha <- (1 - rho) * alpha + rho * newalpha #- rho * df(alpha[1], gamma, k, M) / d2f(alpha[1], gamma, k, M)
    # eta <- eta - rho * df(eta, lambda, V, k) / d2f(eta, lambda, V, k)
    cat("alpha", alpha[1] "\n\n")
    
    beta <- lambda / rowSums(lambda)
    if(i %% 10 == 0){
      j = i %/% 10
      perplexity[j] = calculate_perplexity(testw, alpha, beta)
    }
  }
  return(list(alpha=alpha, eta=eta, gamma = gamma, lambda = lambda, beta=beta, perplexity = perplexity))
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
