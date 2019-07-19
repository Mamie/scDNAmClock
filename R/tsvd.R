#' Denoise using truncated SVD
#' 
#' @inheritParams rsvd::rsvd
#' @return A matrix of denoised data
#' @export
TSVD_denoise <- function(A, k = 100) {
  decomposed <- rsvd::rsvd(A, k =  k, q = 2)
  reconstructed <- decomposed$u %*% diag(decomposed$d) %*% t(decomposed$v)
  colnames(reconstructed) <- colnames(A)
  rownames(reconstructed) <- rownames(A)
  return(reconstructed)
}

#' Decide the rank of the truncated SVD based on Linderman 2018 bioRxiv
#' 
#' Adpated from https://github.com/KlugerLab/ALRA/blob/master/alra.R
#' Heuristic for choosing rank k for the low rank approximation based on
#' statistics of the spacings between consecutive singular values. Finds
#' the smallest singular value $\sigma_i$ such that $\sigma_i - \sigma_{i-1}$
#' is significantly different than spacings in the tail of the singular values.
#'
#' @inheritParams rsvd::rsvd
#' @param n Number of resampling iteration
#' @param pval_thresh The threshold for ``significance''
#' @param noise_start Index for which all smaller singular values are considered noise
#' @return A list with three items: Chosen k, P values of each possible k, 
#' singular values of the matrix A
#' @export
choose_k <- function (A, K=100, pval_thresh=1E-10, noise_start=80, q=2) {
  if (K > min(dim(A))) {
    stop("For an m by n matrix, K must be smaller than the min(m,n).\n")
  }
  if (noise_start >K-5) {
    stop("There need to be at least 5 singular values considered noise.\n")
  }
  noise_svals <- noise_start:K
  rsvd_out <- rsvd::rsvd(A, K, q=q)
  diffs <- diff(rsvd_out$d)
  pvals <- pnorm(diffs, mean(diffs[noise_svals-1]), sd(diffs[noise_svals-1]))
  k <- max(which(pvals < pval_thresh))
  return (list(k = k, pvals = pvals,d = rsvd_out$d))
}