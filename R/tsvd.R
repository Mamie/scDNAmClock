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

#' Decide the rank of the truncated SVD 
#' @inheritParams rsvd::rsvd
#' @param n Number of resampling iteration
decide_k <- function(A, k = 100, n = 10, ...) {
  decomposed <- rsvd::rsvd(A, k =  k, q = 2)
  print(ggplot(data = data.frame(x = 1:k, singular_values = decomposed$d)) +
    geom_point(aes(x = x, y = singular_values), color = "steelblue") +
    scale_y_log10() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank()) +
    ylab("sigular values"))
  return(decomposed$d)
}