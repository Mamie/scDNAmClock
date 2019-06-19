# example from scmeth
library(Melissa)
library()

# download GSE68642
in_dir = "~/Downloads/GSE68642_coverage"
extract_methylation(in_dir, c(1, 2, 4), out_dir = file.path(in_dir, "binarised"))
