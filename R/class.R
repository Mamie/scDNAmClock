setClass("methCall",
         representation(data = "DataFrame"),
         validity = function(object) {
           if (sum(c("chr", "position", "strand", "met_reads", "unmet_reads") %in% colnames(object@data)) < 5) {
             stop("The columns must have chr, position, strand, met_reads, unmet_reads.")
           } else {
             return(TRUE)
           }
         })
