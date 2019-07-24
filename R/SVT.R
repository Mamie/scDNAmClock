# Candes, E.J., Sing-Long, C.A., and Trzasko, J.D.,
# ``Unbiased Risk Estimates for Singular Value Thresholding''
sure_svt <- function (lambda, sigma, X, s = NULL, is_real = T, svThreshold = 1e-8) {
  M = nrow(X)
  N = ncol(X)
  if (is.null(s)) {
    s <- svd(X)$d
    s <- s[1:min(M,N)]
  }
  R <- sure_div_svt(lambda, s, is_real, M, N, svThreshold)
  R  <- -M * N * sigma^2 + sum(pmin(lambda^2, s^2)) + 2*sigma^2*R
  if (!is_real) {
    R  <- R - M * N * sigma^2
  }
  return(R)
}


sure_div_svt <- function(lambda, s, is_real, M, N, svThreshold) {
  # check multiplicities of singular values in a robust manner
  # browser()
  z <- s[-1]
  s <- matrix(c(s[1], 1), nrow = 1)
  Is <- 1
  while(sum(!is.na(z)) > 0) {
    idx <- which(abs(z - s[Is, 1]) < svThreshold)
    if(length(idx) == 0) {
      s <- rbind(s, matrix(c(z[Is], 1), nrow = 1))
      z[Is] <- NA
      Is <- Is + 1
    }
    z[idx] <- NA
    s[Is, 2] <- s[Is, 2] + length(idx)
    #cat(Is, "\n")
  }
  
  # warns the user about using SURE with a non-simple, not-full rank matrix
  if (any(s[, 2] < svThreshold)) 
    stop('   +   [SURE_SVT] Warning: argument might be rank-deficient.\n')
  if( any( s[, 2] > 1 ))
    stop('   +   [SURE_SVT] Warning: argument might have repeated singular values.\n')
  # find singular values above the threshold
  idx_p <- s[, 1] > lambda
  
  if( is_real )
    x  <- div_svt_real(lambda, s, idx_p, M, N)
  else
    x  <- div_svt_complex(lambda, s, idx_p, M, N)
  return(x)
}


div_svt_real <- function(lambda, s, idx_p, M, N) {
  x <- 0
  if(any(idx_p)) {
    x   =   x + sum(0.5*s[idx_p, 2] * (s[idx_p, 2] + 1) )
    x   =   x + sum( (abs(M-N)*s[idx_p, 2] + 0.5*s[idx_p, 2] * (s[idx_p, 2] - 1)) * (pmax(0, s[idx_p, 1] - lambda)/s[idx_p, 1]))
  }
  D = matrix(0, nrow(s), nrow(s))
  for (Ik in 1:nrow(s)) {
    D[, Ik] <- s[Ik, 2] * s[, 2] * s[, 1] * pmax(0, (s[, 1] - lambda)/(s[, 1]^2 - s[Ik, 1]^2))
  }
  D[ is.na(D) | is.infinite(D) | abs(D) > 1e6 ]     =   0
  x <- x + 2*sum(D)
  return(x)
}


div_svt_complex <- function(lambda, s, idx_p, M, N) {
  x <- 0
  if(any(idx_p)) {
    x <- x + sum( s[idx_p, 2]^2 )
    x <- x + sum((2*abs(M-N) + 1 + s[idx_p, 2]*(s[idx_p, 2] - 1))*(pmax(0, s[idx_p, 1] - lambda)/s[idx_p, 1]))
  }
  
  D = matrix(0, nrow(s), nrow(s))
  for (Ik in 1:nrow(s)) {
    D[, Ik] <- s[Ik, 2]*s[, 2]*s[, 1]*pmax(0, s[, 1] - lambda)/(s[, 1]^2 - s[Ik, 1]^2)
  }
  D[ is.na(D) | is.infinite(D) | abs(D) > 1e6] <- 0
  x 	=   x + 4 * sum(D)
  return(x)
}