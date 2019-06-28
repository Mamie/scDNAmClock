# preprocess individual methylation call files
horvath <- list.files("/gpfs/ysm/project/mw957/data/public/Horvath2018/Liver",
                      pattern = "txt.gz", full.names = T) # mm10 
# chr1 Cg 3014928 3 3 1.000 0.444
# chr1 Cg 3037802 13 13 1.000 0.135
# chromosome, dinucleotide, position, methylated counts, total counts, methylation frequency, 95% confidence interval for methylation frequency

# just needed, position, methylated counts, total counts
id <- purrr::map_chr(horvath, ~strsplit(.x, split = "_|/")[[1]][10])

data <- read_meth(
  files = horvath,
  chr_idx = 1,
  pos_idx = 3,
  met_idx = 4, 
  coverage_idx = 5,
  id = id,
  deduplicate = F,
  header = F
)