compute_denomtable_bybox <- function(x_arr, mu_arr, sigma){
  denom_res <- array(dim = dim(x_arr)[-1])
  for(i in seq.int(dim(denom_res)[1])){
    for(j in seq.int(dim(denom_res)[2])){
      for(k in seq.int(dim(denom_res)[3])){
        denom_res[i, j, k] <- eval_proposal(x_arr[, i, j, k], mu_arr[, i, k], sigma)
      }
    }
  }
  log(denom_res)
}

compute_denomtable_byrow <- function(x_arr, mu_arr, sigma){
  denom_res <- array(dim = dim(x_arr)[-1])
  for(i in seq.int(dim(denom_res)[1])){
    for(j in seq.int(dim(denom_res)[2])){
      for(k in seq.int(dim(denom_res)[3])){
        denom_res[i, j, k] <- mean(apply(mu_arr[, i,], 2, function(mu) eval_proposal(x_arr[, i, j, k], mu, sigma)))
      }
    }
  }
  log(denom_res)
}

compute_denomtable_bytable <- function(x_arr, mu_arr, sigma){
  denom_res <- array(dim = dim(x_arr)[-1])
  for(i in seq.int(dim(denom_res)[1])){
    for(j in seq.int(dim(denom_res)[2])){
      for(k in seq.int(dim(denom_res)[3])){
        denom_res[i, j, k] <- mean(apply(apply(mu_arr, 2, identity), 1, function(mu) eval_proposal(x_arr[, i, j, k], mu, sigma)))
      }
    }
  }
  log(denom_res)
}
