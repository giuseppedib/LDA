check_convergence = function(likelihood, i, convergence_threshold){
  if(i == 1){
    return(FALSE)
  }
  else{
    relative_change = (likelihood[i-1] - likelihood[i]) / likelihood[i-1]
    return(ifelse(relative_change < convergence_threshold, TRUE, FALSE))
  }
}
