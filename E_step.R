
# E step for a specific document
E_step = function(gamma, Phi, alpha, beta, W_doc, n_iter){
  for(i in 1:n_iter){
    # loop over vocabulary, but only those words that have count > 0
    subset_words = which(W_doc > 0)
    for(n in subset_words){
      Phi[n, ] = beta[, n] * exp(digamma(gamma)) * W_doc[n]
    }
    # normalise Phi
    Phi = Phi / rowSums(Phi)
    Phi[is.nan(Phi)] = 0
    gamma = alpha + colSums(Phi)
  }
  return(list(gamma = gamma, Phi = Phi))
}
