# 8 test set performance
library(glmnet)
library(survival)

model_650 <- readRDS("data-raw/hyperparameters/650.rds") # code from 20190723
load("data-raw/NewPhenoAge_SVD.RData")

set.seed(120)
CV_SVD = cv.glmnet(model_650, Pheno_InCHIANTI$PredAge, nfolds=10, alpha=0.5, family="gaussian")
fit_SVD = glmnet(model_650,Pheno_InCHIANTI$PredAge, family="gaussian", alpha=0.5, nlambda=100)
CV_coef_SVD <- scDNAmClock:::tidy_coef(CV_SVD, CV_SVD$lambda.min, c("intercept", colnames(model_650)))
dim(CV_coef_SVD) # 494 probes
metrics(data.frame(truth = Pheno_InCHIANTI$PredAge, estimate = DNAmPhenoAge_Train),
        truth, estimate) %>%
  .[,-2] %>%
  kableExtra::kable(., digits = 2) %>%
  kableExtra::kable_styling()

# compute asus as surrogate of MSE to truth to choose k 
# required input - estimator and variance of error from the original dataset ()
# how to estimate variance of error (smallest variance in the probes) 

# test set PhenoAge (difference from the original matrix)
datMeth_WHI_SVD <- TSVD_denoise(datMeth_WHI, k = 650)
DNAmPhenoAge_Test <- as.numeric(predict(fit_SVD, datMeth_WHI_SVD, type="response", s=CV_SVD$lambda.min))
test_perf <- summary(coxph(Surv(Pheno_WHI$ENDFOLLOWALLDY, Pheno_WHI$DEATHALL) ~ DNAmPhenoAge_Test+Pheno_WHI$agewhi))
test_perf$coefficients %>% 
  kableExtra::kable(., digits = 3) %>%
  kableExtra::kable_styling()

