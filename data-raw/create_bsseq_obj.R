#' ---
#' title: "Mouse ESC"
#' author: "Mamie Wang"
#' date: "June 20 2019"
#' ---

#' We are testing the utility functions to preprocess the single cell bisulfite
#' sequencing data. We used scM&T-seq dataset from [Angermueller 2016 Nature Methods](https://www.nature.com/articles/nmeth.3728) as an example.
#' The dataset (GEO accession: [GSE68642](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68642)) provides 81 mouse embryonic stem cells that passed the quality filter, 
#' 16 of which was cultured in '2i' media and 65 in serum. The paper mentioned 
#' that in serum, ESCs will be metastable and heterogenous, switching between  
#' transcriptional states, whereas 2i ESC will have hypomethylation. 

in_dir <- "/gpfs/ysm/project/mw957/data/public/scMT-seq"
cov_files <- list.files(in_dir, pattern = "cov.gz", full.names = T)
id <- purrr::map_chr(cov_files, ~strsplit(.x, split = "[/.]")[[1]][9])

data <- read_meth(
  files = cov_files,
  chr_idx = 1,
  pos_idx = 2,
  met_idx = 5, 
  unmet_idx = 6,
  strand_idx = NULL,
  id = id,
  deduplicate = T)
saveRDS(data, file = "/gpfs/ysm/project/mw957/data/processed/scMT-seq/deduplicated_meth_call.rds")

data <- readRDS("/gpfs/ysm/project/mw957/data/processed/Angermueller2016/deduplicated_meth_call.rds")
joined <- join_meth_list(data, all.x = FALSE, all.y = FALSE) 
saveRDS(joined, "/gpfs/ysm/project/mw957/data/processed/Angermueller2016/inner_joined_meth.rds")

p_list <- plot_coverage(data)
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = "v")
ggsave(p, filename = "~/Desktop/coverage.png", width = 10, height = 3)

