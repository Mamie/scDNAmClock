load("data-raw/Caclulation_DNAmPhenoAge.RData")
pheno_age_dat <- list(CpG = CpGs, intercept = Intercept, weight = Weights$Weight)
devtools::use_data(pheno_age_dat, internal = TRUE)




