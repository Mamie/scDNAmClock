#' @import ggplot2 
setMethod("plot",
           signature = signature("methCall", "character"),
           function(x, y) {
            ggplot(data = data.frame(coverage = x@data$met_reads + x@data$unmet_reads)) +
              geom_histogram(aes(x = coverage), bins = 10, fill = "steelblue") +
              ggtitle(y) +
              theme_classic() +
              scale_x_log10() +
              scale_y_continuous(labels = fancy_scientific) +
              theme(plot.title = element_text(hjust = 0.5))
          })

#' Plot coverage distribution for a list of methCall objects
#' @param meth_list A list of methCall objects
#' @return A list of ggplot objects
#' @export
plot_coverage <- function(meth_list) {
  titles <- names(meth_list)
  p_list <- list()
  for (i in seq_along(meth_list)) {
    p <- plot(meth_list[[i]], titles[i])
    p_list[[titles[i]]] <- p
  }
  return(p_list)
}

# from https://groups.google.com/forum/#!topic/ggplot2/a_xhMoQyxZ4
fancy_scientific <- function(l) { 
  l <- format(l, scientific = TRUE) 
  l <- gsub("e+00", "", l, fixed = T)
  l <- gsub("+", "", l, fixed = T)
  # quote the part before the exponent to keep all the digits 
  l <- gsub("^(.*)e", "'\\1'e", l) 
  # turn the 'e+' into plotmath format 
  l <- gsub("e", "%*%10^", l) 
  # return this as an expression 
  parse(text=l) 
} 