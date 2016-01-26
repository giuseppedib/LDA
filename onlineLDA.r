compute_likelihood = function(gamma, Phi, alpha, beta){
  temp = digamma(gamma) - digamma(sum(gamma))
  # formula from appendix A3
  out = lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha-1)*temp) + 
    sum(t(as.matrix(Phi)) * temp) +     
    sum(Phi * t(ifelse(beta > 0, log(beta), 0))) -   # this bit can be written as   sum( (w[d,]*Phi) * t(log(beta)) )
    ( lgamma(sum(gamma)) - sum(lgamma(gamma)) + sum((gamma-1)*temp) ) - 
    sum(ifelse(Phi>0, Phi * log(Phi), 0))
  return(out)
}




# E step for a specific document
E_step_single_doc = function(gamma, Phi, alpha, beta, W_doc, max_iter, convergence_threshold){
  likelihood = rep(NA, max_iter)
  for(i in 1:max_iter){
    # loop over vocabulary, but only those words that have count > 0
    subset_words = which(W_doc > 0)
    for(n in subset_words){
      Phi[n, ] = beta[, n] * exp(digamma(gamma))
    }
    # normalise Phi
    Phi = Phi / rowSums(Phi)
    Phi[is.nan(Phi)] = 0
    Phi = Phi * W_doc
    gamma = (alpha + colSums(Phi))
    likelihood[i] = compute_likelihood(gamma, Phi, alpha, beta)
    if(check_convergence(likelihood, i, convergence_threshold)) break
  }
  # print(likelihood[1:i])
  return(list(gamma = gamma, Phi = Phi, likelihood = likelihood[i]))
}


# E_step = function(gamma, Phi, alpha, beta, W, max_iter, convergence_threshold){
#   likelihoods = rep(NA, nrow(W))
#   for(i in 1:nrow(W)){
#     obj = E_step_single_doc(gamma[i, ], Phi[, , i], alpha, beta, W[i, ], max_iter, convergence_threshold)
#     likelihoods[i] = obj$likelihood
#     Phi[, , i] = obj$Phi
#     gamma[i, ] = obj$gamma
#   }
#   return(list(gamma = gamma, phi=Phi, likelihood = sum(likelihoods)))
# }

E_step = function(gamma, Phi, alpha, beta, W, max_iter, convergence_threshold){
  likelihoods = rep(NA, nrow(W))
  # Use package foreach to parallelise E step
  `%op%` <- if (getDoParRegistered()) `%dopar%` else `%do%`
  input_list = create_iter_list(gamma, Phi, W)
  res = foreach(input = iter(input_list), .export = c("E_step_single_doc", "compute_likelihood", "check_convergence")) %op% {
    E_step_single_doc(input$gamma, input$Phi, alpha, beta, input$W, max_iter, convergence_threshold)
  }
  for(i in 1:nrow(W)){
    likelihoods[i] = res[[i]]$likelihood
    Phi[, , i] = res[[i]]$Phi
    gamma[i, ] = res[[i]]$gamma
  }
  return(list(gamma = gamma, phi=Phi, likelihood = sum(likelihoods)))
}




                                        # Stochastic Variational EM

mean_lambda <- rep(100,k)
lambda <- matrix(100,k,V)    # Each row is a Dirichlet parameter of lenght V

compute_likelihood = function(gamma, Phi, alpha, lambda, eta){
  E_log_theta = digamma(gamma) - digamma(sum(gamma))
  
  E_log_beta = digamma(lambda) - digamma(rowSums(lambda))
  
  V = ncol(lambda)
  out = sum(Phi * (t(E_log_theta + E_log_beta) - ifelse(Phi>0, log(Phi), 0))) + 
    -lgamma(sum(gamma)) + sum(lgamma(gamma)) + sum((alpha - gamma)*E_log_theta) +
    sum(-lgamma(rowSums(lambda)) + rowSums((eta - lambda)*E_log_beta + lgamma(lambda))) +
    lgamma(sum(alpha)) - sum(lgamma(alpha)) + (lgamma(V*eta) - V*lgamma(eta))
  return(out)
}

onlineLDA <- function(w, alpha, gamma, phi, lambda, eta, max_iter = 100, conv_threshold = 1e-3, MAX_ITER_var_EM = 1000) {
  M <- nrow(gamma)
  conv_likelihood <- 100
  for(i in 1:max_iter){
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
    cat("Finished iter", i, "took", iter2, "steps.\n")
    mean_lambda <- eta + M * t(phi[,,doc])
    lambda <- (1 - rho) * lambda + rho * mean_lambda

    #likelihood <- compute_likelihood(gamma[doc, ], phi[,,doc], alpha, lambda, eta)
    #conv_likelihood <- abs((old_likelihood - likelihood) / old_likelihood)
  }
  return(list(alpha=alpha, gamma = gamma, lambda=lambda))
}
