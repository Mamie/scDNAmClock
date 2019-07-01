#' Number of rows
#' 
#' @param x A methCall object
#' @rdname n_row-methods
#' @export
setGeneric("n_row", function(x) standardGeneric("n_row"))

#' @rdname n_row-methods
setMethod("n_row",
          signature = signature("methCall"),
          function(x) {
            return(nrow(x@data))
          })

#' Filter sites
#' 
#' @param x A methCall object
#' @param f A predicate function that return a logical vector for whether to keep 
#' the sites
#' @rdname filter_sites-methods
#' @export
setGeneric("filter_sites", function(x, f) standardGeneric("filter_sites"))

#' @rdname filter_sites-methods
setMethod("filter_sites",
          signature = signature("methCall", "function"),
          function(x, f) {
            return(new("methCall", data = x@data[f(x),]))
          })

#' Convert methCall object to data frame
#' 
#' @param x A methCall object
#' @rdname as_df-methods
#' @export
setGeneric("as_df", function(x, name) standardGeneric("as_df"))

#' @rdname as_df-methods
setMethod("as_df",
          signature = signature("methCall"),
          function(x) {
            if(n_row(x) == 0) return()
            df <- as.data.frame(x@data)
            return(df)
          })

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

#' Convert methCall object to data frame
#' 
#' @param x A methCall object
#' @rdname as_GRanges-methods
#' @export
setGeneric("as_GRanges", function(x) standardGeneric("as_GRanges"))

#' @rdname as_GRanges-methods
setMethod("as_GRanges",
          signature = signature("methCall"),
          function(x) {
            if(n_row(x) == 0) return()
            else {
              granges <- GenomicRanges::GRanges(seqnames = Rle(x@data$chr),
                                 ranges = IRanges::IRanges(start = as.integer(x@data$position), 
                                                           end = as.integer(x@data$position)),
                                 strand = Rle(x@data$strand))
              return(granges)
            }
          })

#' Liftover methCall object coordinate
#' 
#' @param x A methCall object
#' @param y The path to the chain file
#' @rdname map_coord-methods
#' @export
setGeneric("map_coord", function(x, y) standardGeneric("map_coord"))

#' @rdname map_coord-methods
setMethod("map_coord",
          signature = signature(x = "methCall", 
                                y = "character"),
          function(x, y) {
            grObject <- as_GRanges(x)
            chainObject <- rtracklayer::import.chain(y)
            results <- as.data.frame(rtracklayer::liftOver(grObject, chainObject))
            x@data$chr <- Rle(results$seqnames)
            x@data$position <- Rle(results$start)
            x@data$strand <- Rle(results$strand)
            x
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
