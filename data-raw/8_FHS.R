# test on FHS dataset
load("data-raw/Phenoc1c2-FHS-mortality.Rdata")

View(Phenoc1.FHS.mortality) # 3936
View(Phenoc2.FHS.mortality) # 213

FHSDNAmc2_M <- readRDS('data-raw/FHSDNAmc2_M_svd.rds')

n_lambda <- 50
lambda_max <- 1000
tau <- 0.13
n_tau <- length(tau)
lambda <- matrix(NA, nrow = n_tau, ncol = n_lambda)
for (i in seq_along(tau)) {
  lambda[i,] <- seq(0, tau[i] * lambda_max, length.out = n_lambda)
}

c2_SURE <- matrix(NA, nrow = n_tau, ncol = n_lambda)

for (i in seq(n_tau)) {
  for (j in seq(n_lambda)) {
    c2_SURE[i, j] <- sure_svt(lambda[i,j], tau[i], SVD = FHSDNAmc2_M, is_real = T, svThreshold = 1e-8)
  }
}

i <- 1
p1 <- plot_SURE(lambda[i,], c2_SURE[i,]) + ylab("FHS c2 SURE")
ggsave(p1, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/FHSc2_lamba.png", width = 2, height = 2)

reconstructed_c2 <- SVT_denoise(lambda = lambda[which.min(c2_SURE)], svd = FHSDNAmc2_M) # 443206    213
reconstructed_c2 <- apply(reconstructed_c2, 2, function(x) 2^sinh(x)/(1 + 2^sinh(x)))

# load the original data (M values)
original_c2 <- readRDS("data-raw/FHSDNAmc2_M.rds")
sample_id <- colnames(original_c2)[-1]
probe_id <- original_c2$probe.id
original_c2 <- apply(original_c2[,-1], 2, function(x) 2^sinh(x)/(1 + 2^sinh(x)))

# PhenoAge modeling
rownames(reconstructed_c2) <- probe_id
colnames(reconstructed_c2) <- sample_id
reconstructed_c2_phenoage <- scDNAmClock::PhenoAge(reconstructed_c2)

rownames(original_c2) <- probe_id
original_c2_phenoage <- scDNAmClock::PhenoAge(original_c2)

# comparison of the moribidity 
head(Phenoc2.FHS.mortality)

summary_stats <- data.frame(
  shareid = as.numeric(names(original_c2_phenoage$y)),
  original = original_c2_phenoage$y, 
  SVT = reconstructed_c2_phenoage$y) %>%
  left_join(Phenoc2.FHS.mortality, by = "shareid")


p <- ggplot(data = summary_stats) +
  geom_point(aes(x = original, y = SVT), size = 0.1, alpha = 0.2) +
  theme_bw() +
  geom_abline() +
  xlab("PhenoAge") +
  ylab("PhenoAge SVT") +
  theme(panel.grid = element_blank())
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/FHSc2_SVT_orig.png", width = 3, height = 3)


original_haz <- coxph(Surv(TimetoDeath2, Deathyesno) ~ original + age.upd, data = summary_stats)
SVT_haz <- coxph(Surv(TimetoDeath2, Deathyesno) ~ SVT + age.upd, data = summary_stats)

coefs <- data.frame(original = summary(original_haz)$conf.int[1,c(1,3,4)],
                    SVT = summary(SVT_haz)$conf.int[1,c(1,3,4)])
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
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/FHSc2_cox_SVT.png", width = 4, height = 3)
