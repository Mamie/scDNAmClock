library(glmnet)
library(scDNAmClock)
library(survival)
library(tidyverse) 
library(tidymodels)
library(hqreg)

# retraining the PhenoAge model with TSVD values
load("data-raw/NewPhenoAge_SVD.RData")

# examine the original dataset
set.seed(120)
# cross validation to select the best lambda
CV = cv.glmnet(datMeth_InCHIANTI, Pheno_InCHIANTI$PredAge, nfolds=10, alpha=0.5, family="gaussian")
plot(CV)
CV_coef <- scDNAmClock:::tidy_coef(CV, s = CV$lambda.min, coef_names = c("intercept", colnames(datMeth_InCHIANTI))) 

# original performance 
DNAmPhenoAge_Train <- as.numeric(predict(CV, datMeth_InCHIANTI, type="response", s=CV$lambda.min))
# training error
metrics(data.frame(truth = Pheno_InCHIANTI$PredAge, estimate = DNAmPhenoAge_Train),
        truth, estimate) %>%
  .[,-2] %>%
  kableExtra::kable(., digits = 2) %>%
  kableExtra::kable_styling() # r squared, root mean squared, mean absolute error
# .metric .estimator .estimate
# <chr>   <chr>          <dbl>
# 1 rmse    standard       3.07 
# 2 rsq     standard       0.976
# 3 mae     standard       2.34

DNAmPhenoAge_Test <- as.numeric(predict(fit, datMeth_WHI, type="response", s=CV$lambda.min))
WHI_cox <- coxph(Surv(ENDFOLLOWALLDY, DEATHALL) ~ DNAmPhenoAge_Test + agewhi, data = Pheno_WHI)
summary(WHI_cox)

library(hqreg)
CV_robust <- cv.hqreg(datMeth_InCHIANTI, Pheno_InCHIANTI$PredAge, method = "huber", nfolds=10, alpha=0.5, seed = 120, lambda.min = 0.01)
plot(CV_robust)
CV_coef <- coef(CV_robust, lambda = "lambda.1se") %>%
  as.matrix() %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%  
  setNames(c("term","estimate")) %>%
  filter(estimate != 0)

DNAmPhenoAge_Train <- as.numeric(predict(CV_robust, datMeth_InCHIANTI, lambda = "lambda.1se"))
# training error
metrics(data.frame(truth = Pheno_InCHIANTI$PredAge, estimate = DNAmPhenoAge_Train),
        truth, estimate) %>%
  .[,-2] %>%
  kableExtra::kable(., digits = 2) %>%
  kableExtra::kable_styling()

DNAmPhenoAge_Test <- as.numeric(predict(CV_robust, datMeth_WHI, lambda = "lambda.min"))

WHI_cox <- coxph(Surv(ENDFOLLOWALLDY, DEATHALL) ~ DNAmPhenoAge_Test + agewhi, data = Pheno_WHI)
summary(WHI_cox)




# perform SVD (transform data)

# InCHIANTI_M <- apply(datMeth_InCHIANTI, 2, function(x) asinh(log2(x/(1-x))))
# InCHIANTI_M <- t(apply(InCHIANTI_M, 1, function(x) truncate_outlier(x, a = 3)))
# 
# WHI_M <- apply(datMeth_WHI, 2, function(x) asinh(log2(x/(1-x))))
# WHI_M <- t(apply(WHI_M, 1, function(x) truncate_outlier(x, a = 3)))
# 
# # the following code was ran on HPC
# set.seed(120)
# InCHIANTI_SVD <- rsvd::rsvd(InCHIANTI_M, k = min(dim(InCHIANTI_M)))
# set.seed(120)
# WHI_SVD <- rsvd::rsvd(WHI_M, k = min(dim(WHI_M)))
InCHIANTI_M_asinh <- readRDS("data-raw/InCHIANTI_M_asinh.rds")
InCHIANTI_SVD <- readRDS("data-raw/InCHIANTI_M_SVD.rds")
WHI_M_asinh <- readRDS("data-raw/WHI_M_asinh.rds")
WHI_SVD <- readRDS("data-raw/WHI_M_SVD.rds")

# find a way to estimate the noise level

# select the lambda shrinkage parameter using the SURE
tau <- 0.104 # use the tau from the replicate analysis
n_tau <- length(tau)
n_lambda <- 50
lambda_max <- 500
lambda <- matrix(NA, nrow = n_tau, ncol = n_lambda)
for (i in seq_along(tau)) {
  lambda[i,] <- seq(0, tau[i] * lambda_max, length.out = n_lambda)
}

InCHIANTI_SURE <- WHI_SURE <- matrix(NA, nrow = n_tau, ncol = n_lambda)


for (i in seq(n_tau)) {
  for (j in seq(n_lambda)) {
    InCHIANTI_SURE[i, j] <- sure_svt(lambda[i,j], tau[i], InCHIANTI_M_asinh, s = InCHIANTI_SVD$d, is_real = T, svThreshold = 1e-8)
  }
}

for (i in seq(n_tau)) {
  for (j in seq(n_lambda)) {
    WHI_SURE[i, j] <- sure_svt(lambda[i,j], tau[i], WHI_M_asinh, s = WHI_SVD$d, is_real = T, svThreshold = 1e-8)
  }
}

i = 1
p1 <- plot_SURE(lambda[i,], InCHIANTI_SURE[i,]) + ylab("InCHIANTI SURE") +
  scale_y_continuous(labels = scDNAmClock:::fancy_scientific)
p2 <- plot_SURE(lambda[i,], WHI_SURE[i,]) + ylab("WHI SURE") +
  scale_y_continuous(labels = scDNAmClock:::fancy_scientific)
p <- cowplot::plot_grid(p1, p2, nrow = 2, align = "h")
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/InCHIANTI_WHI_lambda.png", width = 2, height = 4)
InCHIANTI_SVT <- SVT_denoise(InCHIANTI_M_asinh, 
                             lambda = lambda[i, which.min(InCHIANTI_SURE[i,])],
                             svd = InCHIANTI_SVD)
InCHIANTI_beta <- apply(InCHIANTI_SVT, 2, function(x) 2^sinh(x)/(1 + 2^sinh(x)))

CV_robust <- cv.hqreg(InCHIANTI_beta, Pheno_InCHIANTI$PredAge, method = "huber", nfolds=10, alpha=0.5, 
                      lambda.min = 0.01, seed = 120)
plot(CV_robust)
CV_coef <- coef(CV_robust, lambda = "lambda.1se") %>%
  as.matrix() %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%  
  setNames(c("term","estimate")) %>%
  filter(estimate != 0)

DNAmPhenoAge_Train <- as.numeric(predict(CV_robust, InCHIANTI_beta, lambda = "lambda.1se"))
# training error
metrics(data.frame(truth = Pheno_InCHIANTI$PredAge, estimate = DNAmPhenoAge_Train),
        truth, estimate) %>%
  .[,-2] %>%
  kableExtra::kable(., digits = 2) %>%
  kableExtra::kable_styling()

# PhenoAge_InCHIANTI <- scDNAmClock::PhenoAge(t(datMeth_InCHIANTI))$y
# PhenoAge_InCHIANTI_SVT <- scDNAmClock::PhenoAge(t(InCHIANTI_beta))$y
# 
# metrics(data.frame(truth = Pheno_InCHIANTI$PredAge, estimate = PhenoAge_InCHIANTI),
#         truth, estimate) %>%
#   .[,-2] %>%
#   kableExtra::kable(., digits = 2) %>%
#   kableExtra::kable_styling()
# metrics(data.frame(truth = Pheno_InCHIANTI$PredAge, estimate = PhenoAge_InCHIANTI_SVT),
#         truth, estimate) %>%
#   .[,-2] %>%
#   kableExtra::kable(., digits = 2) %>%
#   kableExtra::kable_styling()


# p <- ggplot(data = data.frame(orig = PhenoAge_InCHIANTI, SVT = PhenoAge_InCHIANTI_SVT)) +
#   geom_point(aes(x = orig, y = SVT), size = 0.1, alpha = 0.2) +
#   theme_bw() +
#   geom_abline() +
#   xlab("PhenoAge") +
#   ylab("PhenoAge SVT") +
#   theme(panel.grid = element_blank())
# ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/InCHIANTI_SVT_orig.png", width = 3, height = 3)


WHI_SVT <- SVT_denoise(WHI_M_asinh, 
                       lambda = lambda[i,which.min(WHI_SURE[i,])],
                       svd = WHI_SVD)
WHI_beta <- apply(WHI_SVT, 2, function(x) 2^sinh(x)/(1 + 2^sinh(x)))

PhenoAge_WHI_SVT <- as.numeric(predict(CV_robust, WHI_beta, lambda = "lambda.1se"))
PhenoAge_WHI <- scDNAmClock::PhenoAge(t(datMeth_WHI))$y
PhenoAge_WHI_SVT <- scDNAmClock::PhenoAge(t(WHI_beta))$y


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
