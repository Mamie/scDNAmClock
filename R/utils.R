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

#This is a wrapper function for colorRampPalette. It allows for the
#definition of the number of intermediate colors between the main colors.
#Using this option one can stretch out colors that should predominate
#the palette spectrum. Additional arguments of colorRampPalette can also
#be added regarding the type and bias of the subsequent interpolation.
# from https://menugget.blogspot.com/2011/11/define-color-steps-for-colorramppalette.html#more
color.palette <- function(steps, n.steps.between=NULL, ...){
  
  if(is.null(n.steps.between)) n.steps.between <- rep(0, (length(steps)-1))
  if(length(n.steps.between) != length(steps)-1) stop("Must have one less n.steps.between value than steps")
  
  fill.steps <- cumsum(rep(1, length(steps))+c(0,n.steps.between))
  RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
  RGB[,fill.steps] <- col2rgb(steps)
  
  for(i in which(n.steps.between>0)){
    col.start=RGB[,fill.steps[i]]
    col.end=RGB[,fill.steps[i+1]]
    for(j in seq(3)){
      vals <- seq(col.start[j], col.end[j], length.out=n.steps.between[i]+2)[2:(2+n.steps.between[i]-1)]  
      RGB[j,(fill.steps[i]+1):(fill.steps[i+1]-1)] <- vals
    }
  }
  
  new.steps <- rgb(RGB[1,], RGB[2,], RGB[3,], maxColorValue = 255)
  pal <- colorRampPalette(new.steps, ...)
  return(pal)
}