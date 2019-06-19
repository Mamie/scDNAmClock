#' Extract the relevant columns from the files
#' 
#' @param dir The directory containing tsv.gz files for methylation status
#' @param col_idx The index of column for chromosome, start, methylation level
#' @param out_dir The output directory
#' @import dplyr
#' @export
extract_methylation <- function(dir, col_idx, out_dir) {
  files <- list.files(dir, pattern = "gz")
  p <- dplyr::progress_estimated(length(files))
  for (file in files) {
    data <- read.delim(file.path(dir, file), sep = "\t", header = F)
    data <- data[, col_idx]
    readr::write_tsv(data, path = file.path(out_dir, file), col_names = F)
    p$pause(0.1)$tick()$print()
  }
}

# find_intersecting_regions <- function(dir, col_idx) {
#   files <- list.files(dir, pattern = "tsv.gz")
#   for (file in files) {
#     data <- readr::read_tsv(file.path(dir, file), col_names = F)
#     data <- data[, col_idx]
#     names <- 
#   }
# }