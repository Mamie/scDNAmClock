#' Extract the relevant columns from the files
#' 
#' @param files The path to tsv files for methylation status
#' @param chr_idx The index of column for chromosome
#' @param pos_idx The index of column for position
#' @param met_idx The index of column for the number of methylated reads
#' @param unmet_idx The index of column for the number of unmethylated reads
#' @param strand_idx The index of column containing strand if available
#' @param id The id to each file
#' @param deduplicate Whether to sum the reads for rows with the same genomic location
#' @importFrom S4Vectors Rle DataFrame SimpleList
#' @import dplyr
#' @export
read_meth <- function(files, chr_idx, pos_idx, met_idx, unmet_idx, strand_idx = NULL,
                      id, deduplicate = T) {
  meth_list <- SimpleList()
  
  p <- dplyr::progress_estimated(length(files))
  for (i in seq_along(files)) {
    data <- read.delim(files[i], sep = "\t", header = F, stringsAsFactors = F)
    if (!is.null(strand_idx)) {
      strand <- data[, strand_idx]
    } else {
      strand <- rep("*", nrow(data))
    }
    
    if (deduplicate)
      data <- data.frame(chr = data[, chr_idx],
                       position = data[, pos_idx],
                       strand = strand,
                       met_reads = data[, met_idx],
                       unmet_reads = data[, unmet_idx]) %>%
      group_by(chr, position, strand) %>%
      summarize_all(sum)
    
    meth <- new("methCall", 
                data = with(data, 
                            DataFrame(
                              chr = Rle(chr),
                              position = Rle(position),
                              strand = Rle(strand),
                              met_reads = Rle(met_reads),
                              unmet_reads = Rle(unmet_reads))))
                 
    
    meth_list[[id[i]]] <- meth
    p$pause(0.1)$tick()$print()
  }
  
  return(meth_list)
}

#' Join a list of methCall object by site
#' 
#' @param meth_list A named list of methCall object (from read_meth)
#' @param all.x Whether to include all rows in x
#' @param all.y Whether to include all rows in y
#' @return A DataFrame that is a inner join of all the methCall
#' @export
join_meth_list <- function(meth_list, all.x = FALSE, all.y = FALSE) {
  browser()
  joined <- x[[1]]@data
  colnames(joined)[4:5] <- paste0(colnames(joined)[4:5], ".", names(x)[1])
  for (i in seq_along(x)[-1]) {
    df <- x[[i]]@data
    colnames(df)[4:5] <- paste0(colnames(df)[4:5], ".", names(x)[i])
    joined <- SparkR::merge(joined, df, by = c("chr", "position", "strand"),
                            all.x = FALSE, all.y = FALSE)
  }
  return(joined)
}