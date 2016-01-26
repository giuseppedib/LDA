# V <- 10              # number of words in the vocabulary
# M <- 4               # number of documents in the corpus
# k <- 5               # number of topics
# w <- matrix(0,M,V)   # matrix of counts
# 
# (alpha <- rep(1,k) )
# (beta <- rdirichlet(k, rep(1,V)) )
# (phi <- array(1/V, dim = c(V,k,M)) )        
# (gamma <- matrix( alpha + V/k, M, k) )

# update beta
update_beta = function(phi, w, k){
  beta = matrix(NA, k, ncol(w))
  for (i in seq(1:k))
    beta[i,] <- colSums(t(phi[,i,]) * (w>0))
  beta <- beta / rowSums(beta)
  return(beta)
}

# update alpha
Dlalpha <- function(alpha, gamma, M) {
  return (M * (digamma(sum(alpha)) - digamma(alpha)) + colSums(digamma(gamma) - digamma(rowSums(gamma))))
}


update_alpha_vec = function(alpha, gamma, M){
  MAXITER <- 100
  EPSILON <- 0.00001
  step = 0
  conv = 10
  log_alpha = log(alpha)
  while (step < MAXITER && conv > EPSILON) {
    Dalpha <- Dlalpha(alpha, gamma, M)
    gradf <- Dalpha*alpha
    H <- M*(trigamma(sum(alpha)) - diag(trigamma(alpha)))
    grad2f <- H * alpha + diag(gradf)
    #hvalue <- h(alpha, M)
    #c <- sum( Dalpha / hvalue ) / (1 / (M*trigamma(sum(alpha))) + sum(1 / hvalue) )
    #HinvD <- (Dalpha - c) / hvalue
    HinvD <- solve(H, gradf)
    log_alpha <- log_alpha - HinvD
    alpha <- exp(log_alpha)
    alpha <- alpha - HinvD
    
    conv <- sum(HinvD * HinvD)
    step <- step + 1
    cat(alpha,"\n")
  } 
  cat(alpha, "\n")
  return(alpha)
}


#####################################
### Follows the 1D alpha optimisation
#####################################

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
  
  init_a = 100
  log_a = log(init_a)
  cat("Optimising alpha with 1-dim Newton-Rhapson.\n")
  while (step < MAXITER && conv > EPSILON) {
    a = exp(log_a)
    fvalue = f(a, gamma, length(alpha), M)
    dfvalue = df(a, gamma, length(alpha), M)[1]
    d2fvalue = d2f(a, gamma, length(alpha), M)[1]
    log_a = log_a - dfvalue/(d2fvalue * a + dfvalue)
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
  cat("Eta optimisation completed.\n")
  out = exp(log_a)
  return(out)
}
