# 9 meta analysis of hazard ratio
# https://stats.stackexchange.com/questions/158985/meta-analysis-of-trials-given-hazard-ratios

library(metafor)

orig <- data.frame(logOR = 1, SE = 1)
SVT <- data.frame(logOR = 1, SE = 1)

orig_res <- rma(yi = logOR, sei = SE, data = orig, method = "FE")
SVT_res <- rma(yi = logOR, sei = SE, data = SVT, method = "FE")

# save the results in a csv file