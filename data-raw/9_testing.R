# 8 test set performance
library(glmnet)
library(survival)

model_650 <- readRDS("data-raw/hyperparameters/650.rds") # code from 20190723
load("data-raw/NewPhenoAge_SVD.RData")

dim(datMeth_InCHIANTI)
summary(apply(datMeth_InCHIANTI, 2, sd))

nl <- 50
taus <- c(1.4e-3, 1.5e-3, 1.6e-3, 1.7e-3)
lambdas <- matrix(NA, length(taus), nl)   # thresholds
SURE <- matrix(NA, length(taus), nl)
rownames(SURE) <- taus
colnames(SURE) <- 1:nl
lambda_max <- 50

sd <- apply(datMeth_InCHIANTI, 2, sd)
quantile(sd, probs = c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2))  %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

InCHIANTI_SVD <- rsvd::rsvd(datMeth_InCHIANTI, k = min(dim(datMeth_InCHIANTI)))
saveRDS(InCHIANTI_SVD, file = "data-raw/InCHIANTI_SVD.rds")

for (i in seq_along(taus)) {
  cat("tau:", taus[i], "\n")
  lambdas[i,] <- seq(0, lambda_max * taus[i], length.out = nl) 
  for (j in 1:nl) {
    cat("lambda:", lambdas[i, j], "\n")
    SURE[i, j] <- sure_svt(lambdas[i, j], taus[i], datMeth_InCHIANTI, s = InCHIANTI_SVD$d, is_real = T, svThreshold = 1e-8) 
  }
}

SURE_data <- as.data.frame(SURE) 
SURE_data$tau <- taus
colnames(SURE_data) <- c(1:nl, "tau")

SURE_data <- SURE_data %>%
  tidyr::gather(lambda, sure, -tau) %>%
  mutate(lambda = as.numeric(lambda)) %>%
  arrange(tau, lambda) %>%
  mutate(lambda_val = as.numeric(t(lambdas)))
ggplot(data = SURE_data %>%
         filter(tau == 0.017)) +
  geom_point(aes(x = lambda_val, y = sure, group = tau)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  #facet_wrap(~tau, scale = "free") +
  ylab("SURE") +
  xlab("lambda") +
  scale_x_continuous(labels = scDNAmClock:::fancy_scientific) 
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/choice_lambda.png", width = 6, height = 3) 

SURE_data %>%
  group_by(tau) %>% filter(sure == min(sure))
# get an estimate for noise level


shrinked <- 

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

# test sure_svt