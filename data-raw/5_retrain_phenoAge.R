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
InCHIANTI_SVD <- rsvd::rsvd(datMeth_InCHIANTI, k = min(dim(datMeth_InCHIANTI)))
set.seed(120)
WHI_SVD <- rsvd::rsvd(datMeth_WHII, k = min(dim(datMeth_WHI)))
InCHIANTI_SVD <- readRDS("data-raw/InCHIANTI_SVD.rds")
WHI_SVD <- readRDS("data-raw/WHI_SVD.rds")

# select the lambda shrinkage parameter using the SURE
tau <- 0.013
lambda <- seq(0, tau * 200, length.out = 100)
InCHIANTI_SURE <- c()
WHI_SURE <- c()
for (l in lambda) {
  InCHIANTI_SURE <- c(InCHIANTI_SURE, sure_svt(l, tau, datMeth_InCHIANTI, s = InCHIANTI_SVD$d, is_real = T, svThreshold = 1e-8))
}

for (l in lambda) {
  WHI_SURE <- c(WHI_SURE, sure_svt(l, tau, datMeth_WHI, s = WHI_SVD$d, is_real = T, svThreshold = 1e-8))
}

p1 <- plot_SURE(lambda, InCHIANTI_SURE) + ylab("InCHIANTI SURE")
p2 <- plot_SURE(lambda, WHI_SURE) + ylab("WHI SURE")
p <- cowplot::plot_grid(p1, p2, nrow = 2, align = "h")
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/InCHIANTI_WHI_lambda.png", width = 2, height = 4)
InCHINATI_SVT <- SVT_denoise(datMeth_InCHIANTI, 
                             lambda = lambda[which.min(InCHIANTI_SURE)],
                             svd = InCHIANTI_SVD)
PhenoAge_InCHIANTI <- scDNAmClock::PhenoAge(t(datMeth_InCHIANTI))$y
PhenoAge_InCHIANTI_SVT <- scDNAmClock::PhenoAge(t(InCHINATI_SVT))$y

metrics(data.frame(truth = Pheno_InCHIANTI$PredAge, estimate = PhenoAge_InCHIANTI),
        truth, estimate) %>%
  .[,-2] %>%
  kableExtra::kable(., digits = 2) %>%
  kableExtra::kable_styling()
metrics(data.frame(truth = Pheno_InCHIANTI$PredAge, estimate = PhenoAge_InCHIANTI_SVT),
        truth, estimate) %>%
  .[,-2] %>%
  kableExtra::kable(., digits = 2) %>%
  kableExtra::kable_styling()


p <- ggplot(data = data.frame(orig = PhenoAge_InCHIANTI, SVT = PhenoAge_InCHIANTI_SVT)) +
  geom_point(aes(x = orig, y = SVT), size = 0.1, alpha = 0.2) +
  theme_bw() +
  geom_abline() +
  xlab("PhenoAge") +
  ylab("PhenoAge SVT") +
  theme(panel.grid = element_blank())
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/InCHIANTI_SVT_orig.png", width = 3, height = 3)


WHI_SVT <- SVT_denoise(datMeth_WHI, 
                       lambda = lambda[which.min(WHI_SURE)],
                       svd = WHI_SVD)
PhenoAge_WHI <- scDNAmClock::PhenoAge(t(datMeth_WHI))$y
PhenoAge_WHI_SVT <- scDNAmClock::PhenoAge(t(WHI_SVT))$y


p <- ggplot(data = data.frame(orig = PhenoAge_WHI, SVT = PhenoAge_WHI_SVT)) +
  geom_point(aes(x = orig, y = SVT), size = 0.1, alpha = 0.2) +
  theme_bw() +
  geom_abline() +
  xlab("PhenoAge") +
  ylab("PhenoAge SVT") +
  theme(panel.grid = element_blank())
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/WHI_SVT_orig.png", width = 3, height = 3)

WHI_cox <- coxph(Surv(ENDFOLLOWALLDY, DEATHALL) ~ PhenoAge_WHI + agewhi, data = Pheno_WHI)
WHI_cox_res <- summary(WHI_cox)

WHI_cox_SVT <- coxph(Surv(ENDFOLLOWALLDY, DEATHALL) ~ PhenoAge_WHI_SVT + agewhi, data = Pheno_WHI)
WHI_cox_SVT_res <- summary(WHI_cox_SVT)
coefs <- data.frame(original = WHI_cox_res$conf.int[1,c(1,3,4)],
           SVT = WHI_cox_SVT_res$conf.int[1,c(1,3,4)])
coefs$var <- c("hazard ratio", "lwr", "upr")
p <- coefs %>%
  tidyr::gather(model, value, -var) %>%
  tidyr::spread(var, value) %>%
  ggplot(data = .) +
  geom_point(aes(x = model, y = `hazard ratio`)) +
  geom_errorbar(aes(x = model, ymin = lwr, ymax = upr), width = 0.2) +
  coord_flip() +
  theme_bw() +
  theme(axis.line = element_blank(), panel.border = element_blank(),
        axis.title.y = element_blank()) +
  geom_hline(aes(yintercept = 1), linetype = "dashed")
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/WHI_cox_SVT.png", width = 4, height = 3)

  
  
  
  
set.seed(120)
CV_SVD = cv.glmnet(InCHINATI_SVT, Pheno_InCHIANTI$PredAge, nfolds=10, alpha=0, family="gaussian")
fit_SVD = glmnet(InCHINATI_SVT,Pheno_InCHIANTI$PredAge, family="gaussian", alpha=0, nlambda=100)
# CV_coef_SVD <- scDNAmClock::tidy_coef(CV_SVD, CV_SVD$lambda.min, c("intercept", colnames(datMeth_InCHIANTI)))
# nrow(CV_coef_SVD) # 389 

#This Predicts the Age in the Training data
DNAmPhenoAge_SVD_Train <- as.numeric(predict(fit_SVD, InCHINATI_SVT,type="response",s=CV_SVD$lambda.min))

metrics(data.frame(truth = Pheno_InCHIANTI$PredAge, estimate = DNAmPhenoAge_SVD_Train),
        truth, estimate) %>%
  .[,-2] %>%
  kableExtra::kable(., digits = 2) %>%
  kableExtra::kable_styling()


DNAmPhenoAge_SVD_Test <- as.numeric(predict(fit_SVD, WHI_SVT, type="response",s=CV_SVD$lambda.min))

test_perf <- summary(coxph(Surv(Pheno_WHI$ENDFOLLOWALLDY, Pheno_WHI$DEATHALL) ~ DNAmPhenoAge_SVD_Test+Pheno_WHI$agewhi))
test_perf$coefficients %>%
  kableExtra::kable(., digits = 3) %>%
  kableExtra::kable_styling()

#Plot the predicted Age versus the DNAmPhenoAge
WGCNA::verboseScatterplot(Pheno_WHI$agewhi,DNAmPhenoAge_SVD_Test,xlab="Age",ylab="New DNAmPhenoAge",main="Validation Set")
abline(0,1,col="red")
#Plot New vs. Original DNAmPhenoAge
WGCNA::verboseScatterplot(DNAmPhenoAge_SVD_Test, Pheno_WHI$DNAmPhenoAge,xlab="New DNAmPhenoAge",ylab="Original DNAmPhenoAge",main="Training Set")
abline(0,1,col="red")



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
