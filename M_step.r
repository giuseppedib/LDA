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
    beta[i,] <- colSums(t(phi[,i,]) * w)
  beta <- beta / rowSums(beta)
  return(beta)
}

# update alpha
Dlalpha <- function(alpha, gamma) {
  return (M * (digamma(sum(alpha)) - digamma(alpha)) + colSums(digamma(gamma) - digamma(rowSums(gamma))))
}

Mtrialpha <- function(alpha) {
  return (M * trigamma(alpha)) 
}

update_alpha = function(alpha, gamma){
  MAXITER <- 1000
  EPSILON <- 0.00001
  step = 0
  conv = 10
  while (step < MAXITER && conv > EPSILON) {
    Dalpha <- Dlalpha(alpha, gamma)
    Malpha <- Mtrialpha(alpha)
    c <- sum( Dalpha / Malpha ) / (1 / trigamma(sum(alpha)) + sum(1 / Malpha) )
    HinvD <- (Dalpha - c) / Malpha
    alpha <- alpha - HinvD
    
    conv <- sum(HinvD * HinvD)
    step <- step + 1
    cat(alpha ,"\n")
  } 
  return(alpha)
}
