# setMethod("coverage_plot",
#           signature = "SimpleList",
#           function(object) {
#             for (i in seq_along(object)) {
#               
#             }
#           })

#' @import ggplot2 
.coverage_plot <- function(methCall_obj, id) {
  ggplot(data = data.frame(coverage = methCall_obj$met_reads + methCall_obj$unmet_reads)) +
    geom_histogram(aes(x = coverage), bins = 10, fill = "steelblue") +
    ggtitle(id) +
    theme_classic() +
    scale_x_log10() +
    scale_y_continuous(labels = fancy_scientific)
}

# from https://groups.google.com/forum/#!topic/ggplot2/a_xhMoQyxZ4
fancy_scientific <- function(l) { 
  browser()
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