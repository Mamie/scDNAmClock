library(dplyr)

# evaluate the missingness from the dataset
inner_joined <- readRDS("~/Downloads/mm9_mm10_inner_join.rds")
dim(inner_join) # 12689082      307
# 152 samples before filtering

# compute the methylation percentage
as.data.frame(inner_joined) %>%
  tidyr::gather(category, value, -c(chr, position, strand)) %>%
  tidyr::separate(category, c("type", "geo_accession"), sep = "[.]") %>%
  tidyr::spread(type, value) %>%
  mutate(`% methylated` = met_reads/coverage) %>%
  select(chr, position, geo_accession, `% methylated`) %>%
  tidyr::spread(geo_accession, `% methylated`)
  