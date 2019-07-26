library(glmnet)
library(scDNAmClock)
library(survival)
library(tidyverse)
library(tidymodels)

# retraining the PhenoAge model with TSVD values
load("data-raw/NewPhenoAge_SVD.RData")

# examine the original dataset
set.seed(120)
# cross validation to select the best lambda
CV = cv.glmnet(datMeth_InCHIANTI,Pheno_InCHIANTI$PredAge, nfolds=10, alpha=0.5, family="gaussian")
CV_coef <- scDNAmClock:::tidy_coef(CV, s = CV$lambda.min, coef_names = c("intercept", colnames(datMeth_InCHIANTI))) 

fit = glmnet(datMeth_InCHIANTI, Pheno_InCHIANTI$PredAge, family="gaussian", alpha=0.5, nlambda=100)

# original performance 
DNAmPhenoAge_Train <- as.numeric(predict(fit, datMeth_InCHIANTI, type="response", s=CV$lambda.min))
# training error
metrics(data.frame(truth = Pheno_InCHIANTI$PredAge, estimate = DNAmPhenoAge_Train),
        truth, estimate) %>%
  .[,-2] %>%
  kableExtra::kable(., digits = 2) %>%
  kableExtra::kable_styling()# r squared, root mean squared, mean absolute error
# .metric .estimator .estimate
# <chr>   <chr>          <dbl>
# 1 rmse    standard       3.07 
# 2 rsq     standard       0.976
# 3 mae     standard       2.34

DNAmPhenoAge_Test <- as.numeric(predict(fit, datMeth_WHI, type="response", s=CV$lambda.min))
head(DNAmPhenoAge_Test)
head(Pheno_WHI$DNAmPhenoAge)

# correlation of age with predicted PhenoAge
WGCNA::verboseScatterplot(Pheno_WHI$agewhi,DNAmPhenoAge_Test,xlab="Age",ylab="New DNAmPhenoAge",main="Validation Set")
abline(0,1,col="red") # correlation of predicted phenoage with actual age in test set: 0.66 
WGCNA::verboseScatterplot(Pheno_WHI$agewhi,Pheno_WHI$DNAmPhenoAge,xlab="Age",ylab="Original DNAmPhenoAge",main="Validation Set") # ? not the original DNAmPhenoAge
abline(0,1,col="red")

# perform SVD 
set.seed(120)
k <- 500
datMeth_InCHIANTI_k <- choose_k(datMeth_InCHIANTI, K = k, pval_thresh = 1e-10, noise_start = k * 0.8)
InCHIANTI_SVD <- rsvd::rsvd(datMeth_InCHIANTI, k = min(dim(datMeth_InCHIANTI)))

datMeth_WHI_k <- choose_k(datMeth_WHI, K = k, pval_thresh = 1e-10, noise_start = k * 0.8)
datMeth_WHI_SVD <- TSVD_denoise(datMeth_WHI, k = datMeth_WHI_k$k)

# examine the difference between the denoised dataset vs the original dataset
levine_clocks <- scDNAmClock:::pheno_age_dat$CpG


# arrange the patients by predicted PhenoAge
SVD_normed <- apply(datMeth_WHI_SVD[,levine_clocks], 2, scale)
col_dend <- hclust(dist(t(SVD_normed)), method = "ward.D2") %>%
  as.dendrogram() %>%
  ladderize
row_dend <- hclust(dist(SVD_normed), method = "ward.D2") %>%
  as.dendrogram() %>%
  ladderize

png(file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/svd_WHI_heatmap.png", width = 600, height = 500)
gplots::heatmap.2(SVD_normed,
                  trace = "none", 
                  dendrogram = "column",
                  Colv = col_dend, 
                  Rowv = row_dend, labRow = F, labCol = F, 
                  col = colorRampPalette(c("#0047BB", "white", "#D1350F"))(10), 
                  breaks = seq(-5, 5, 1), 
                  margins = c(1, 1))
dev.off()


normed <- apply(datMeth_WHI[,levine_clocks], 2, scale)

png(file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/orig_WHI_heatmap.png", width = 600, height = 500)
gplots::heatmap.2(normed, 
                  trace = "none", 
                  Colv = col_dend, 
                  Rowv = row_dend, labRow = F, labCol = F, 
                  col = colorRampPalette(c("#0047BB", "white", "#D1350F"))(10), breaks = seq(-5, 5, 1),  
                  margins = c(1, 1))
dev.off()


set.seed(120)
CV_SVD = cv.glmnet(datMeth_InCHIANTI_SVD, Pheno_InCHIANTI$PredAge, nfolds=10, alpha=0.5, family="gaussian")
fit_SVD = glmnet(datMeth_InCHIANTI_SVD,Pheno_InCHIANTI$PredAge, family="gaussian", alpha=0.5, nlambda=100)

CV_coef_SVD <- scDNAmClock::tidy_coef(CV_SVD, CV_SVD$lambda.min, c("intercept", colnames(datMeth_InCHIANTI_SVD)))
plot(CV_SVD, main="CV Elastic Net")

#This Predicts the Age in the Training data
DNAmPhenoAge_SVD_Train <- as.numeric(predict(fit_SVD, datMeth_InCHIANTI_SVD,type="response",s=CV_SVD$lambda.min))
DNAmPhenoAge_SVD_Test <- as.numeric(predict(fit_SVD, datMeth_WHI_SVD,type="response",s=CV_SVD$lambda.min))

metrics(data.frame(truth = Pheno_InCHIANTI$PredAge, estimate = DNAmPhenoAge_SVD_Train),
        truth, estimate) %>%
  .[,-2] %>%
  kableExtra::kable(., digits = 2) %>%
  kableExtra::kable_styling()

#Plot the predicted Age versus the DNAmPhenoAge
WGCNA::verboseScatterplot(Pheno_WHI$agewhi,DNAmPhenoAge_SVD_Test,xlab="Age",ylab="New DNAmPhenoAge",main="Validation Set")
abline(0,1,col="red")
#Plot New vs. Original DNAmPhenoAge
WGCNA::verboseScatterplot(DNAmPhenoAge_SVD_Test, Pheno_WHI$DNAmPhenoAge,xlab="New DNAmPhenoAge",ylab="Original DNAmPhenoAge",main="Training Set")
abline(0,1,col="red")

test_perf <- summary(coxph(Surv(Pheno_WHI$ENDFOLLOWALLDY, Pheno_WHI$DEATHALL) ~ DNAmPhenoAge_SVD_Test+Pheno_WHI$agewhi))
test_perf$coefficients %>% 
  kableExtra::kable(., digits = 3) %>%
  kableExtra::kable_styling()

original_perf <- summary(coxph(Surv(Pheno_WHI$ENDFOLLOWALLDY, Pheno_WHI$DEATHALL) ~ Pheno_WHI$DNAmPhenoAge+Pheno_WHI$agewhi))
original_perf$coefficients %>% 
  kableExtra::kable(., digits = 3) %>%
  kableExtra::kable_styling()

length(intersect(CV_coef_SVD$term, scDNAmClock:::pheno_age_dat$CpG)) 

# hyperparameter search (k vs parameters in the elastic net model)
# performance of the data on the original model
DNAmPhenoAge_Train_2 <- as.numeric(predict(fit, datMeth_InCHIANTI_SVD, type="response", s=CV$lambda.min))
DNAmPhenoAge_Test <- as.numeric(predict(fit, datMeth_WHI_SVD, type="response", s=CV$lambda.min))
metrics(data.frame(truth = Pheno_InCHIANTI$PredAge, estimate = DNAmPhenoAge_Train_2),
        truth, estimate) %>%
  .[,-2] %>%
  kableExtra::kable(., digits = 2) %>%
  kableExtra::kable_styling()
test_perf_2 <- summary(coxph(Surv(Pheno_WHI$ENDFOLLOWALLDY, Pheno_WHI$DEATHALL) ~ DNAmPhenoAge_Test+Pheno_WHI$agewhi))
test_perf_2$coefficients %>% 
  kableExtra::kable(., digits = 3) %>%
  kableExtra::kable_styling()
