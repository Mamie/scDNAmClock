#' Extract the relevant columns from the files
#' 
#' @param files The path to tsv files for methylation status
#' @param chr_idx The index of column for chromosome
#' @param pos_idx The index of column for position
#' @param met_idx The index of column for the number of methylated reads
#' @param id The id to each file
#' @param coverage_idx The index of column containing coverage; optional 
#' @param unmet_idx The index of column for the number of unmethylated reads
#' @param strand_idx The index of column containing strand; optional 
#' @param deduplicate Whether to sum the reads for rows with the same genomic location
#' @param header Whether to there is a header
#' @param sep delimiter for methylation status file
#' @importFrom S4Vectors Rle DataFrame SimpleList
#' @import dplyr
#' @export
read_meth <- function(files, chr_idx, pos_idx, met_idx, id, coverage_idx = NULL, strand_idx = NULL,
                      unmet_idx = NULL, deduplicate = T, header = T, sep = '\t') {
  meth_list <- SimpleList()
  
  p <- dplyr::progress_estimated(length(files))
  for (i in seq_along(files)) {
    data <- read.delim(files[i], sep = sep, header = header, stringsAsFactors = F)
    if (!is.null(strand_idx)) {
      strand <- data[, strand_idx]
    } else {
      strand <- rep("*", nrow(data))
    }
    
    
    if (!is.null(coverage_idx)) coverage <- data[, coverage_idx]
    if (!is.null(unmet_idx)) coverage <- data[, met_idx] + data[, unmet_idx]
    data <- data.frame(chr = gsub("chr", "", data[, chr_idx]),
                      position = data[, pos_idx],
                      strand = strand,
                      met_reads = data[, met_idx],
                      coverage = coverage)

    if (deduplicate)
      data <- dplyr::distinct(data) %>%
        group_by(chr, position, strand) %>%
        summarize_all(sum)
    
    meth <- new("methCall", 
                data = with(data, 
                            DataFrame(
                              chr = Rle(chr),
                              position = Rle(position),
                              strand = Rle(strand),
                              met_reads = Rle(met_reads),
                              coverage = Rle(coverage))))
    
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
  joined <- meth_list[[1]]@data
  colnames(joined)[4:5] <- paste0(colnames(joined)[4:5], ".", names(meth_list)[1])
  for (i in seq_along(meth_list)[-1]) {
    df <- meth_list[[i]]@data
    colnames(df)[4:5] <- paste0(colnames(df)[4:5], ".", names(meth_list)[i])
    joined <- SparkR::merge(joined, df, by = c("chr", "position", "strand"),
                            all.x = FALSE, all.y = FALSE)
  }
  return(joined)
}
