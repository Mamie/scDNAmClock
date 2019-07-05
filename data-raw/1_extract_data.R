# Levine CpG clock
cpg_clocks <- readr::read_csv("~/Dropbox/600 Presentations/Yale projects/scMethylseq/figs/cpgAllonetime.cpgInfo.csv")
Levine_clock  <- cpg_clocks %>%
  filter(`Levine-DNAm-PhenoAge` == "Levine") # 513 clock genes
Levine_clock$prob

# samples that are technical replicates


# data from GSE55763 (the normalized beta values)
readr::read_tsv("/gpfs/ysm/project/mw957/data/public/Lehne2015/GSE55763_normalized_betas.txt.gz", col_names = T) %>%
  filter(ID_REF %in% Levine_clock$prob) %>%
  select()