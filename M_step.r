V <- 10              # number of words in the vocabulary
M <- 4               # number of documents in the corpus
k <- 5               # number of topics
w <- matrix(0,M,V)   # matrix of counts

MAXITER <- 1000
EPSILON <- 0.00001

library("gtools")

(alpha <- rep(1,k) )
(beta <- rdirichlet(k, rep(1,V)) )
(phi <- array(1/V, dim = c(V,k,M)) )        
(gamma <- matrix( alpha + V/k, M, k) )

                                        # M - step
# update beta
for (i in seq(1:k))  
  beta[i,] <- colSums(t(phi[,i,]) * w)
beta <- beta / rowSums(beta)

# update alpha
Dlalpha <- function(alpha, gamma) {
  return (M * (digamma(sum(alpha)) - digamma(alpha)) + colSums(digamma(gamma) - digamma(rowSums(gamma))))
}

Mtrialpha <- function(alpha) {
  return (M * trigamma(alpha)) 
}

step <- 0
conv <- 10
while (step < MAXITER && conv > EPSILON) {
  Dalpha <- Dlalpha(alpha, gamma)
  Malpha <- Mtrialpha(alpha)
  c <- sum( Dalpha / Malpha ) / (1 / trigamma(sum(alpha)) + sum(1 / Malpha) )
  HinvD <- (Dalpha - c) / Malpha
  alpha <- alpha - HinvD

  conv <- sum(HinvD * HinvD)
  step <- step + 1
  cat(conv ,"\n")
} 
alpha
step
conv


                                        # log_likelihood lower bound

#loglike lowbound <- log(gamma(sum(alpha))) -
#  sum((alpha - 1) * (digamma(gamma) - digamma(sum(gamma))) - log(gamma(alpha))) +
#  sum( (digamma(gamma) - digamma(sum(gamma))) * 
