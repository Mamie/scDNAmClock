#' Get coefficients from glmnet
#' 
#' @param x A glmnet object
#' @param s lambda value for glmnet
#' @param coef_names The names of the coefficients
#' @import dplyr
tidy_coef <- function(x, s, coef_names = NULL){
  coefs <- coef(x, s = s) %>%
    matrix() %>%   
    data.frame() %>%  
    tibble::rownames_to_column() %>%  
    setNames(c("term","estimate")) 
  if (!is.null(coef_names)) {
    coefs$term <- coef_names
  }
  coefs %>%
    filter(estimate != 0)
}