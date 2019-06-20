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
  browser()
  meth_list <- SimpleList()
  
  p <- dplyr::progress_estimated(length(files))
  for (i in seq_along(files)) {
    data <- data.table::fread(files[i], sep = "\t", header = F, stringsAsFactors = F)
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
