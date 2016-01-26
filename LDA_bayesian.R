LDA = function(W, n_topics = 3, max_iter = 100, max_iter_E_step = 30, convergence_threshold=1e-6){
  k = 7
  V = ncol(W)
  M = nrow(W)
  alpha <- rep(1, k)
  # beta <- rdirichlet(k, rep(1, V))
  eta <- 1
  lambda <- matrix(rgamma(k*V, 100, 100), nrow=k, ncol=V)
  gamma <- matrix((alpha + V/k), M, k)
  phi <- array(0, dim = c(V, k, M))
  
  likelihood = rep(NA, max_iter)
  for(i in 1:max_iter){
    # E step
    obj = E_step(gamma, phi, alpha, lambda, eta, W, max_iter_E_step, convergence_threshold)
    
    likelihood[i] = obj$likelihood
    phi = obj$phi
    gamma = obj$gamma
    # M step
    # update lambda
    for(j in 1:k){
      lambda[j, ] = eta + colSums(t(phi[, j, ]) * W)
    }
    alpha = update_alpha(alpha, gamma, M)
    # eta = update_eta(eta, lambda, k)
    
    cat(sprintf("Iteration %d of EM completed. Likelihood: %1.3f \n", i, likelihood[i]))
    #cat("Max beta", apply(beta, 1, max), "\n")
    if(check_convergence(likelihood, i, convergence_threshold)) break
  }
  return(list(likelihood = likelihood[1:i], phi = obj$phi, gamma = obj$gamma, alpha = alpha, eta=eta, lambda=lambda))
}

compute_likelihood = function(gamma, Phi, alpha, lambda, eta, temp, E_log_theta, E_log_beta, M){
  V = ncol(lambda)
  out = temp + 
    -lgamma(sum(gamma)) + sum(lgamma(gamma)) + sum((alpha - gamma)*E_log_theta) +
    sum(-lgamma(rowSums(lambda)) + rowSums((eta - lambda)*E_log_beta + lgamma(lambda))) / M +
    lgamma(sum(alpha)) - sum(lgamma(alpha)) + (lgamma(V*eta) - V*lgamma(eta)) / M
  return(out)
}


check_convergence = function(likelihood, i, convergence_threshold){
  if(i == 1){
    return(FALSE)
  }
  else{
    relative_change = (likelihood[i-1] - likelihood[i]) / likelihood[i-1]
    return(ifelse(relative_change < convergence_threshold, TRUE, FALSE))
  }
}

# E step for a specific document
E_step_single_doc = function(gamma, Phi, alpha, lambda, eta, W_doc, M, max_iter, convergence_threshold){
  likelihood = rep(NA, max_iter)
  for(i in 1:max_iter){
    # loop over vocabulary, but only those words that have count > 0
    subset_words = which(W_doc > 0)
    E_log_theta = digamma(gamma) - digamma(sum(gamma))
    temp_likelihood_contribution = matrix(0, ncol(lambda), nrow(lambda))
    E_log_beta = matrix(0, nrow(lambda), ncol(lambda))
    for(n in subset_words){
      E_log_beta[, n] = digamma(lambda[, n]) - digamma(rowSums(lambda))
      Phi[n, ] = exp(E_log_theta + E_log_beta[, n])
      temp_likelihood_contribution[n, ] = E_log_theta + E_log_beta[, n]
    }
    # normalise Phi
    Phi = Phi / rowSums(Phi)
    Phi[is.nan(Phi)] = 0
    Phi = Phi * W_doc
    gamma = (alpha + colSums(Phi))
    temp_likelihood_contribution = sum(Phi*(temp_likelihood_contribution - ifelse(Phi>0, log(Phi), 0)))
    likelihood[i] = compute_likelihood(gamma, Phi, alpha, lambda, eta, temp_likelihood_contribution, E_log_theta, E_log_beta, M)
    if(check_convergence(likelihood, i, convergence_threshold)) break
  }
  # print(likelihood[1:i])
  return(list(gamma = gamma, Phi = Phi, likelihood = likelihood[i]))
}


E_step = function(gamma, Phi, alpha, lambda, eta, W, max_iter, convergence_threshold){
  likelihoods = rep(NA, nrow(W))
  for(i in 1:nrow(W)){
    obj = E_step_single_doc(gamma[i, ], Phi[, , i], alpha, lambda, eta, W[i, ], nrow(W), max_iter, convergence_threshold)
    likelihoods[i] = obj$likelihood
    Phi[, , i] = obj$Phi
    gamma[i, ] = obj$gamma
  }
  return(list(gamma = gamma, phi=Phi, likelihood = sum(likelihoods)))
}
# 
# E_step = function(gamma, Phi, alpha, beta, W, max_iter, convergence_threshold){
#   likelihoods = rep(NA, nrow(W))
#   # Use package foreach to parallelise E step
#   `%op%` <- if (getDoParRegistered()) `%dopar%` else `%do%`
#   input_list = create_iter_list(gamma, Phi, W)
#   res = foreach(input = iter(input_list), .export = c("E_step_single_doc", "compute_likelihood", "check_convergence")) %op% {
#     E_step_single_doc(input$gamma, input$Phi, alpha, beta, input$W, max_iter, convergence_threshold)
#   }
#   for(i in 1:nrow(W)){
#     likelihoods[i] = res[[i]]$likelihood
#     Phi[, , i] = res[[i]]$Phi
#     gamma[i, ] = res[[i]]$gamma
#   }
#   return(list(gamma = gamma, phi=Phi, likelihood = sum(likelihoods)))
# }



f = function(a, gamma, K, M){
  ss = sum(colSums(digamma(gamma) - digamma(rowSums(gamma))))
  return(M * (lgamma(K * a) - K * lgamma(a)) + (a - 1) * ss)
}

df = function(a, gamma, K, M){
  ss = sum(colSums(digamma(gamma) - digamma(rowSums(gamma))))
  out = M * (K * digamma(K * a) - K * digamma(a)) + ss
  return(out)
}

d2f = function(a, gamma, K, M){
  out = M * (K * K * trigamma(K * a) - K * trigamma(a))
  return(out)
}

update_alpha = function(alpha, gamma, M){
  MAXITER <- 1000
  EPSILON <- 1e-6
  step = 0
  conv = 10
  
  init_a = 200
  log_a = log(init_a)
  cat("Optimising alpha with 1-dim Newton-Rhapson.\n")
  while (step < MAXITER && conv > EPSILON) {
    a = exp(log_a)
    fvalue = f(a, gamma, length(alpha), M)
    dfvalue = df(a, gamma, length(alpha), M)[1]
    d2fvalue = d2f(a, gamma, length(alpha), M)[1]
    log_a = log_a - dfvalue/(d2fvalue * a + dfvalue);
    conv = abs(dfvalue)
    cat("\t", a, "\n")
    step = step + 1
  } 
  cat("Alpha optimisation completed.\n")
  out = rep(exp(log_a), length(alpha))
  return(out)
}

update_eta = function(eta, lambda, M){
  MAXITER <- 1000
  EPSILON <- 1e-6
  step = 0
  conv = 10
  
  init_eta = 100
  log_a = log(init_eta)
  cat("Optimising eta with 1-dim Newton-Rhapson.\n")
  while (step < MAXITER && conv > EPSILON) {
    a = exp(log_a)
    fvalue = f(a, lambda, length(eta), M)
    dfvalue = df(a, lambda, length(eta), M)
    d2fvalue = d2f(a, lambda, length(eta), M)
    log_a = log_a - dfvalue/(d2fvalue * a + dfvalue);
    conv = abs(dfvalue)
    cat("\t", a, "\n")
    step = step + 1
  } 
  cat("Alpha optimisation completed.\n")
  out = exp(log_a)
  return(out)
}
