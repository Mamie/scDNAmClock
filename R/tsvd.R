#' Denoise using truncated SVD
#' 
#' @inheritParams rsvd::rsvd
#' @return A matrix of denoised data
#' @export
SVT_denoise <- function(lambda, svd) {
  d <- soft_threshold(svd$d, lambda)
  reconstructed <- svd$u %*% diag(d) %*% t(svd$v)
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

plot_SURE <- function(lambda, SURE) {
  df <- data.frame(lambda = lambda, SURE = SURE)
  min_lambda <- df$lambda[which.min(df$SURE)]
  ggplot(data = df) +
    geom_point(aes(x = lambda, y = SURE), size = 0.2) +
    geom_vline(aes(xintercept = min_lambda), size = 0.3) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(bquote(lambda[min] == .(round(min_lambda, 3))))  +
    xlab(bquote(lambda))
}
