# Stochastic Variational EM

source("M_step.R")

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


onlineLDA_EM <- function(w, testw, n_iter = 100, conv_threshold = 1e-3, MAX_ITER_var_EM = 1000) {
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
  for(i in 1:n_iter){
    rho <- (2 + i)**(-0.75)
    doc <- sample(c(1:M),1)
    gamma[doc,] <- alpha + ncol(w)/length(alpha)
    iter2 <- 0
    conv_value <- 100
    # old_likelihood <- compute_likelihood(gamma[doc, ], phi[,,doc], alpha, lambda, eta)
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
    alpha <- update_alpha(alpha, gamma, M)
    eta <- update_eta(eta, lambda, M)
    
    beta <- lambda / rowSums(lambda)
    perplexity = rep(NA, n_iter %/% 100)
    if(i %% 100 == 0){
      j = i %/% 100
      perplexity[j] = calculate_perplexity(test_W, alpha, beta)
      cat("Iter", i, "Test set perplexity", perplexity[j], "\n")
    }

    #likelihood <- compute_likelihood(gamma[doc, ], phi[,,doc], alpha, lambda, eta)
    #conv_likelihood <- abs((old_likelihood - likelihood) / old_likelihood)
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
