# preprocess individual methylation call files
# chr1 Cg 3014928 3 3 1.000 0.444
# chr1 Cg 3037802 13 13 1.000 0.135
# chromosome, dinucleotide, position, methylated counts, total counts, methylation frequency, 95% confidence interval for methylation frequency

# just needed, position, methylated counts, total counts
horvath <- list.files("/gpfs/ysm/project/mw957/data/public/Horvath2018/Liver",
                      pattern = "txt.gz", full.names = T) # mm10 

id <- purrr::map_chr(horvath, ~strsplit(.x, split = "_|/")[[1]][10])
data <- read_meth(
  files = horvath,
  chr_idx = 1,
  pos_idx = 3,
  met_idx = 4, 
  coverage_idx = 5,
  id = id,
  deduplicate = F,
  header = F,
  sep = " "
)

saveRDS(data, file = "/gpfs/ysm/project/mw957/data/processed/mouse_liver/Horvath2018.rds")

# Hahn2017 dataset
# 4       3050345 3050345 100     1       0
# 4       3050346 3050346 33.3333333333333        2       4
# BS-seq: The Bismark CpG coverage report is tab-delimited, uses 1-based genomic coordinates for every covered cytosine position in the experiment and is in the following format: <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count non-methylated>
hahn <- list.files("/gpfs/ysm/project/mw957/data/public/Hahn2017/",
                   pattern = "txt.gz", full.names = T)
id <- purrr::map_chr(hahn, ~strsplit(.x, split = "_|/")[[1]][10])
data <- read_meth(
  files = hahn,
  chr_idx = 1,
  pos_idx = 2,
  met_idx = 5, 
  unmet_idx = 6, 
  id = id,
  deduplicate = F,
  header = F
)
saveRDS(data, file = "/gpfs/ysm/project/mw957/data/processed/mouse_liver/Hahn2017.rds")

# Cole2017 dataset
# chr1	3000573	1	0
# chr1	3000574	1	0
# chr <tab> position <tab> methylated reads <tab> unmethylated reads)
cole <- list.files("/gpfs/ysm/project/mw957/data/public/Cole2017",
                   pattern = "txt.gz", full.names = T)
id <- purrr::map_chr(cole, ~strsplit(.x, split = "_|/")[[1]][9])
data <- read_meth(
  files = cole,
  chr_idx = 1,
  pos_idx = 2,
  met_idx = 3, 
  unmet_idx = 4, 
  id = id,
  deduplicate = F,
  header = F
)

data <- purrr::map(data, ~map_coord(.x, chain_file = "/gpfs/ysm/project/mw957/data/public/liftOver_chain/mm9ToMm10.over.chain"))

saveRDS(data, file = "/gpfs/ysm/project/mw957/data/processed/mouse_liver/Cole2017.rds")

