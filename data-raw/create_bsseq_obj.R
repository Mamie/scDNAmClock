# download GSE68642
in_dir <- "~/Downloads/GSE68642_coverage"
cov_files <- list.files(in_dir, pattern = "cov.gz", full.names = T)
id <- purrr::map_chr(cov_files, ~strsplit(.x, split = "[/.]")[[1]][6])

test <- read_meth(cov_files[1:2],
          chr_idx = 1,
          pos_idx = 2,
          met_idx = 5, 
          unmet_idx = 6,
          strand_idx = NULL,
          id = id[1:2]) # 129.4 vs 51.6 vs 49.2

.coverage_plot(test[[1]]@data, names(test)[1])
test <- merge_meth()
