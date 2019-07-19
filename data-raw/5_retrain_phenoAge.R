library(glmnet)
library(scDNAmClock)
library(survival)

# retraining the PhenoAge model with TSVD values
load("data-raw/NewPhenoAge_SVD.RData")

# examine the original dataset
set.seed(120)
CV= cv.glmnet(datMeth_InCHIANTI,Pheno_InCHIANTI$PredAge, nfolds=10, alpha=0.5, family="gaussian")
sum(CV_coef$estimate != 0)
lambda.Min = CV$lambda.min
pos <- match(colnames(datMeth_InCHIANTI)[which(CV_coef$estimate != 0) - 1], scDNAmClock:::pheno_age_dat$CpG)

fit = glmnet(datMeth_InCHIANTI, Pheno_InCHIANTI$PredAge, family="gaussian", alpha=0.5, nlambda=100)
CV_coef <- scDNAmClock:::tidy_coef(CV, s = CV$lambda.min, coef_names = c("intercept", colnames(datMeth_InCHIANTI))) %>%
  filter(estimate != 0)
ref <- data.frame(term = scDNAmClock:::pheno_age_dat$CpG, estimate = scDNAmClock:::pheno_age_dat$weight) %>%
  arrange(term) # same results as the reference
plot(CV, main="CV Elastic Net")
DNAmPhenoAge_Train=as.numeric(predict(fit,datMeth_InCHIANTI,type="response",s=lambda.Min))
DNAmPhenoAge_Test=as.numeric(predict(fit,datMeth_WHI,type="response",s=lambda.Min))
WGCNA::verboseScatterplot(Pheno_WHI$agewhi,DNAmPhenoAge_Test,xlab="Age",ylab="New DNAmPhenoAge",main="Validation Set")
abline(0,1,col="red")
WGCNA::verboseScatterplot(Pheno_WHI$agewhi,Pheno_WHI$DNAmPhenoAge,xlab="Age",ylab="Original DNAmPhenoAge",main="Validation Set")
abline(0,1,col="red")

coxph(Surv(Pheno_WHI$ENDFOLLOWALLDY, Pheno_WHI$DEATHALL) ~ DNAmPhenoAge_Test+Pheno_WHI$agewhi)
coxph(Surv(Pheno_WHI$ENDFOLLOWALLDY, Pheno_WHI$DEATHALL) ~ Pheno_WHI$DNAmPhenoAge+Pheno_WHI$agewhi)


# perform SVD 
js.svd <- jackstraw::jackstraw_rpca(datMeth_InCHIANTI)
datMeth_InCHIANTI_SVD <- TSVD_denoise(datMeth_InCHIANTI, k = 50)
datMeth_WHI_SVD <- TSVD_denoise(datMeth_WHI, k = 50)

CV_SVD = cv.glmnet(datMeth_InCHIANTI_SVD,Pheno_InCHIANTI$PredAge, nfolds=10, alpha=0.5, family="gaussian")
fit = glmnet(datMeth_InCHIANTI_SVD,Pheno_InCHIANTI$PredAge, family="gaussian", alpha=0.5, nlambda=100)
lambda.Min = CV$lambda.min

#This is plotting the error as a function of the lambda value. It will also show how many CpGs are selected as a function of lambda
plot(CV, main="CV Elastic Net")

#This Predicts the Age in the Training data
DNAmPhenoAge_Train=as.numeric(predict(fit,datMeth_InCHIANTI_SVD,type="response",s=lambda.Min))
DNAmPhenoAge_Test=as.numeric(predict(fit,datMeth_WHI_SVD,type="response",s=lambda.Min))

#Plot the predicted Age versus the DNAmPhenoAge
WGCNA::verboseScatterplot(Pheno_WHI$agewhi,DNAmPhenoAge_Test,xlab="Age",ylab="New DNAmPhenoAge",main="Validation Set")
abline(0,1,col="red")
WGCNA::verboseScatterplot(Pheno_WHI$agewhi,Pheno_WHI$DNAmPhenoAge,xlab="Age",ylab="Original DNAmPhenoAge",main="Validation Set")
abline(0,1,col="red")

#Plot New vs. Original DNAmPhenoAge
WGCNA::verboseScatterplot(DNAmPhenoAge_Test,Pheno_WHI$DNAmPhenoAge,xlab="New DNAmPhenoAge",ylab="Original DNAmPhenoAge",main="Training Set")
abline(0,1,col="red")