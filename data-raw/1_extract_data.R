library(dplyr)

# select samples that are technical replicates from series family file
samples <- readr::read_tsv("~/Downloads/GSE55763.txt", col_names = F)
samples <- t(samples)
samples <- samples[-1,]
colnames(samples) <- c("name", "geo_accession", "dataset", "gender", "age", "desc")
samples <- as.data.frame(samples)
samples$name <- gsub("Peripheral blood, ", "", samples$name)
samples$dataset <- gsub("dataset: ", "", samples$dataset)
samples$gender <- gsub("gender: ", "", samples$gender)
samples$age <- as.numeric(gsub("age: ", "", samples$age))

technical_rep_samples <- samples %>%
  filter(grepl("technical replication study", dataset))

technical_rep_samples$group <- as.numeric(purrr::map_chr(technical_rep_samples$desc, 
                                              ~stringr::str_match(.x, "group ([0-9])")[2]))
technical_rep_samples$sample <- as.numeric(purrr::map_chr(technical_rep_samples$desc,
                                               ~stringr::str_match(.x, "sample ([0-9]+)")[2]))
readr::write_csv(technical_rep_samples[, c(1, 2, 4, 5, 7, 8)],
                path = "data-raw/technical_rep_sample_lists.csv")
  
# data from GSE55763 (the normalized beta values)
samples <- readr::read_csv("/gpfs/ysm/project/mw957/repos/scDNAmClock/data-raw/technical_rep_sample_lists.csv")

# Levine CpG clock
cpg_clocks <- readr::read_csv("/gpfs/ysm/project/mw957/repos/scDNAmClock/data-raw/cpgAllonetime.cpgInfo.csv")
Levine_clock  <- cpg_clocks %>%
  filter(`Levine-DNAm-PhenoAge` == "Levine") # 513 clock genes


tech_rep_data <- readr::read_tsv("/gpfs/ysm/project/mw957/data/public/Lehne2015/GSE55763_normalized_betas.txt.gz", col_names = T) 
tech_rep_data <- tech_rep_data %>%
  filter(ID_REF %in% Levine_clock$prob) %>%
  select_(.dots = c("ID_REF", samples$name))

readr::write_tsv(tech_rep_data, file = "/gpfs/ysm/project/mw957/data/processed/Lehne2015/tech_rep_data.tsv")

