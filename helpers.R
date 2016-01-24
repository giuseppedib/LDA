check_convergence = function(likelihood, i, convergence_threshold){
  if(i == 1){
    return(FALSE)
  }
  else{
    relative_change = (likelihood[i-1] - likelihood[i]) / likelihood[i-1]
    return(ifelse(relative_change < convergence_threshold, TRUE, FALSE))
  }
}

create_iter_list = function(gamma, Phi, W){
  iter_list = list()
  for(i in 1:nrow(W)){
    iter_list[[i]] = list(gamma = gamma[i, ], Phi = Phi[, , i], W = W[i, ])
  }
  return(iter_list)
}
