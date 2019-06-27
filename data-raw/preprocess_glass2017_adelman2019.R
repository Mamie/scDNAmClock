library(ggplot2)
library(dplyr)
# preprocess Glass2017 and Adelman2019 ERRBS dataset
sample_info <- readr::read_csv("data-raw/glass2017_adelman2019_sample_info.csv")
table(sample_info$sex, sample_info$cell_type)


p <- sample_info %>%
  mutate(age = as.numeric(age)) %>%
ggplot(data = .) +
  geom_histogram(aes(x = age, fill = sex)) +
  facet_wrap(~cell_type, ncol = 1) +
  theme_bw() +
  theme(panel.background = element_blank(), 
        strip.background = element_blank(),
        panel.grid = element_blank())
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/scRNAseq_RRBS/figs/eRRBS_sample.png",
       width = 5, height = 4)

adelman2019 <- list.files("/gpfs/ysm/project/mw957/data/public/Adelman2019",
                          pattern = "myCpG.txt.gz")

# format 
# chrBase	chr	base	strand	coverage	freqC	freqT
# chr1.10631	chr1	10631	F	112	95.536	4.464
# chr1.10542	chr1	10542	F	80	100.000	0.000
#
# C represents methylated cytosine
# T represents unmethylated cytosine

id <- purrr::map_chr(adelman2019, ~strsplit(.x, split = "_|/")[[1]][3])

data <- read_meth(
  files = cov_files,
  chr_idx = 2,
  pos_idx = 3,
  met_idx = 6, 
  unmet_idx = 7,
  strand_idx = 4,
  coverage_idx = 5,
  id = id,
  deduplicate = F,
  header = T
)
saveRDS(data, file = "/gpfs/ysm/project/mw957/data/processed/glass2017_adelman2019/adelman2019_meth_call.rds")


glass2017 <- list.files("/gpfs/ysm/project/mw957/data/public/Glass2017", 
                        pattern = "mincov10.txt.gz")

# format
# chrBase	chr	base	strand	coverage	freqC	freqT
# chr1.10497	chr1	10497	F	501	97.01	2.99
# chr1.10563	chr1	10563	F	11	100.00	0.00

id <- purrr::map_chr(glass2017, ~strsplit(.x, split = "_|/")[[1]][3])

data <- read_meth(
  files = cov_files,
  chr_idx = 2,
  pos_idx = 3,
  met_idx = 6, 
  unmet_idx = 7,
  strand_idx = 4,
  coverage_idx = 5,
  id = id,
  deduplicate = F,
  header = T
)
saveRDS(data, file = "/gpfs/ysm/project/mw957/data/processed/glass2017_adelman2019/glass2017_meth_call.rds")