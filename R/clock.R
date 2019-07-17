#' Compute DNAm PhenoAge
#' @param data A matrix with rows as CpGs and columns as samples (rownames should
#' be probe id)
#' @return PhenoAge prediction
#' @export
PhenoAge <- function(data) {
  common <- intersect(rownames(data), pheno_age_dat$CpG)
  cat("Number of CpGs", length(common), "of 513\n")
  matched <- match(pheno_age_dat$CpG[pheno_age_dat$CpG %in% rownames(data)], 
                   rownames(data))
  Ax <- data[matched,] * pheno_age_dat$weight[pheno_age_dat$CpG %in% rownames(data)[matched]]
  DNAm_PhenoAge <- colSums(Ax) + pheno_age_dat$intercept
  return(list(y = DNAm_PhenoAge, Ax = Ax, b = pheno_age_dat$intercept))
}