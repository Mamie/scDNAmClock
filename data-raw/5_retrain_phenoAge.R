library(glmnet)
library(scDNAmClock)
library(survival)

# retraining the PhenoAge model with TSVD values
load("~/Downloads/NewPhenoAge_SVD.RData")

# examine the original dataset
set.seed(120)
CV = cv.glmnet(datMeth_InCHIANTI,Pheno_InCHIANTI$PredAge, nfolds=10, alpha=0.5, family="gaussian")
fit = glmnet(datMeth_InCHIANTI, Pheno_InCHIANTI$PredAge, family="gaussian", alpha=0.5, nlambda=100)
CV_coef <- scDNAmClock:::tidy_coef(CV, s = CV$lambda.min, coef_names = c("intercept", colnames(datMeth_InCHIANTI))) %>%
  filter(estimate != 0)
plot(CV, main="CV Elastic Net")

DNAmPhenoAge_Train <- as.numeric(predict(fit, datMeth_InCHIANTI, type="response", s=CV$lambda.min))
DNAmPhenoAge_Test <- as.numeric(predict(fit, datMeth_WHI, type="response", s=CV$lambda.min))
# correlation of age with predicted PhenoAge
WGCNA::verboseScatterplot(Pheno_WHI$agewhi,DNAmPhenoAge_Test,xlab="Age",ylab="New DNAmPhenoAge",main="Validation Set")
abline(0,1,col="red") 
WGCNA::verboseScatterplot(Pheno_WHI$agewhi,Pheno_WHI$DNAmPhenoAge,xlab="Age",ylab="Original DNAmPhenoAge",main="Validation Set")
abline(0,1,col="red")

# perform SVD 
set.seed(120)
k <- round(0.5 * min(dim(datMeth_InCHIANTI)), 0)
datMeth_InCHIANTI_k <- choose_k(datMeth_InCHIANTI, K = k, pval_thresh = 1e-3, noise_start = k * 0.8)
k <- round(0.5 * min(dim(datMeth_WHI)), 0)
datMeth_WHI_k <- choose_k(datMeth_WHI, K = k, pval_thresh = 1e-3, noise_start = k * 0.8)
datMeth_InCHIANTI_SVD <- TSVD_denoise(datMeth_InCHIANTI, datMeth_InCHIANTI_k$k)
datMeth_WHI_SVD <- TSVD_denoise(datMeth_WHI, k = datMeth_WHI_k$k)
set.seed(120)
CV_SVD = cv.glmnet(datMeth_InCHIANTI_SVD, Pheno_InCHIANTI$PredAge, nfolds=10, alpha=0.5, family="gaussian")
fit_SVD = glmnet(datMeth_InCHIANTI_SVD,Pheno_InCHIANTI$PredAge, family="gaussian", alpha=0.5, nlambda=100)

CV_coef_SVD <- scDNAmClock::tidy_coef(CV_SVD, CV_SVD$lambda.min, c("intercept", colnames(datMeth_InCHIANTI_SVD))) %>%
  filter(estimate != 0)
#This is plotting the error as a function of the lambda value. It will also show how many CpGs are selected as a function of lambda
plot(CV_SVD, main="CV Elastic Net")

#This Predicts the Age in the Training data
DNAmPhenoAge_SVD_Train=as.numeric(predict(fit_SVD, datMeth_InCHIANTI_SVD,type="response",s=CV_SVD$lambda.min))
DNAmPhenoAge_SVD_Test=as.numeric(predict(fit_SVD, datMeth_WHI_SVD,type="response",s=CV_SVD$lambda.min))

#Plot the predicted Age versus the DNAmPhenoAge
WGCNA::verboseScatterplot(Pheno_WHI$agewhi,DNAmPhenoAge_SVD_Test,xlab="Age",ylab="New DNAmPhenoAge",main="Validation Set")
abline(0,1,col="red")

#Plot New vs. Original DNAmPhenoAge
WGCNA::verboseScatterplot(DNAmPhenoAge_SVD_Test, Pheno_WHI$DNAmPhenoAge,xlab="New DNAmPhenoAge",ylab="Original DNAmPhenoAge",main="Training Set")
abline(0,1,col="red")

coxph(Surv(Pheno_WHI$ENDFOLLOWALLDY, Pheno_WHI$DEATHALL) ~ DNAmPhenoAge_SVD_Test+Pheno_WHI$agewhi)
coxph(Surv(Pheno_WHI$ENDFOLLOWALLDY, Pheno_WHI$DEATHALL) ~ Pheno_WHI$DNAmPhenoAge+Pheno_WHI$agewhi)

length(intersect(CV_coef_SVD$term, scDNAmClock:::pheno_age_dat$CpG)) # just 89 sites (49 are in the original model)

