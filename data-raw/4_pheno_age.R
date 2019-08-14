# evaluate the age prediction using the phenoage model
library(scDNAmClock)
library(dplyr)
library(ggplot2)

samples <- readr::read_csv("data-raw/technical_rep_sample_lists.csv")
M_data <- readRDS("data-raw/Lehne2015_M_Levine_1000_random.rds")
all_data <- M_data
all_data[,-1] <- apply(all_data[,-1], 2, function(x) 2^(x)/(2^x + 1))
data <- all_data[, c("ID_REF", samples$name)]
colnames(data)[1] <- "probe"

mapped <- data %>%
  tidyr::gather(pid, beta, -probe) %>%
  left_join(samples[, -2], by = c("pid" = "name")) %>%
  select(-pid) %>%
  tidyr::spread(group, beta)
colnames(mapped)[5:6] <- c("replicate 1", "replicate 2")

probes_data <- split(mapped, mapped$probe)
probe_ICC <- purrr::map(probes_data, 
                        ~irr::icc(.x[, c(5, 6)], 
                                  model = "oneway",
                                  type = "agreement",
                                  unit = "single"))
probe_ICC <- data.frame(
  probe = names(probes_data),
  ICC = purrr::map_dbl(probe_ICC, ~.x$value)
)
probe_ICC$ICC[probe_ICC$ICC < 0] <- 0

# phenoage prediction
data_mat <- as.matrix(data[,-1])
rownames(data_mat) <- data$probe
pheno_age <- scDNAmClock::PhenoAge(data_mat)

plot_data <- data.frame(name = names(pheno_age$y), 
                        pheno_age = unname(pheno_age$y)) %>%
  left_join(samples, by = "name") %>%
  select(-c(name, geo_accession)) %>%
  tidyr::spread(group, pheno_age) %>%
  mutate(diff = `1` - `2`, abs_diff = abs(`1` - `2`))

cols <- c("DNAm PhenoAge" = "steelblue", 
          "Actual age" = "orange")

p <- ggplot(data = plot_data) +
  geom_abline(alpha = 0.6, linetype = "dashed") +
  geom_abline(intercept = 2, alpha = 0.3, linetype = "dashed") +
  geom_abline(intercept = -2, alpha = 0.3, linetype = "dashed") +
  geom_point(aes(x = `1`, y = `2`, color = "DNAm PhenoAge"), size = 3) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(name = "", values = cols) +
  theme(legend.position = c(0.7, 0.2)) +
  xlab("Replicate 1") +
  ylab("Replicate 2") 
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/phenoage2.pdf", width = 6, height = 5)  
  

p <- ggplot(data = plot_data) +
  geom_histogram(aes(x = abs_diff), fill = "steelblue", bins = 10, alpha = 0.8) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("Absolute PhenoAge replicate difference") +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) 
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/phenoage_hist.pdf", width = 4, height = 2)  
res <- t(as.matrix(summary(plot_data$abs_diff)))
kableExtra::kable(res, digits = 2) %>%
  kableExtra::kable_styling()

# examine variability contribution
summary_data <- as.data.frame(pheno_age$Ax) %>% 
  mutate(probe = rownames(pheno_age$Ax)) %>%
  tidyr::gather(name, value, -probe) %>%
  left_join(samples, by = "name") %>%
  select(-c(name, geo_accession, gender, age)) %>%
  tidyr::spread(group, value) %>%
  mutate(abs_deviation = abs(`1` - `2`)) %>%
  group_by(probe) %>%
  summarize(abs_deviation = mean(abs_deviation)) %>%
  left_join(data.frame(probe = scDNAmClock:::pheno_age_dat$CpG,
                       weight = scDNAmClock:::pheno_age_dat$weight), by = "probe") %>%
  left_join(probe_ICC, by = "probe")

p <- ggplot(data = summary_data, 
            aes(x = weight, y = abs_deviation, color = ICC)) +
  geom_point(alpha = 0.9) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.9, 0.2)) +
  ylab("weight * mean replicate difference") +
  scale_color_viridis_c() +
  xlab("weight")
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/probe_variation.pdf", width = 5, height = 4.5)

# annotation
annot_probes <- summary_data %>%
  filter(abs_deviation > 0.4)

p <- p + geom_text(data = annot_probes, aes(label = probe), position = position_jitter(width = 0.01, height = 0.01))
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/probe_variation.pdf", width = 7, height = 4.5)

probe_ICC$category <- "excellent (> 0.9)"
probe_ICC$category[probe_ICC$ICC < 0.9] <- "good (0.75 - 0.9)"
probe_ICC$category[probe_ICC$ICC < 0.75] <- "moderate (0.5 - 0.75)"
probe_ICC$category[probe_ICC$ICC < 0.5] <- "poor (< 0.5)"

summary_data <- left_join(summary_data, probe_ICC) %>%
  mutate(category = factor(category, 
                           levels = c("excellent (> 0.9)",
                                      "good (0.75 - 0.9)",
                                      "moderate (0.5 - 0.75)",
                                      "poor (< 0.5)")))

p <- ggplot(data = summary_data,
            aes(x = category, y = abs_deviation, color = ICC)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.9) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.1, 0.8)) +
  ylab("weight * mean replicate difference") +
  scale_color_viridis_c()
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/probe_variation_box.pdf", width = 7, height = 4.5)

p <- ggplot(data = summary_data,
            aes(x = category, y = weight, color = ICC)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.9) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.1, 0.8)) +
  ylab("weight") +
  scale_color_viridis_c()
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/probe_weight_box.pdf", width = 7, height = 4.5)

summary_data$abs_dev_div_weight <- summary_data$abs_deviation/abs(summary_data$weight)
summary_data %>%
  arrange(desc(abs_dev_div_weight))

# detailed source of variation 
anova_data <- as.data.frame(pheno_age$Ax) %>%
  mutate(probe = rownames(pheno_age$Ax)) %>%
  tidyr::gather(name, Ax, -probe) %>%
  dplyr::left_join(select(samples, c(name, group, sample))) %>%
  select(-name) %>% 
  tidyr::spread(group, Ax) %>%
  mutate(`rep1_rep2` = abs(`1` - `2`))

p1 <- ggplot(data = anova_data) +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed") +
  geom_point(aes(x = reorder(probe, rep1_rep2), y = rep1_rep2), size = 0.1, alpha = 0.5,
             color = "steelblue") +
  theme_classic() +
  coord_flip() +
  theme(panel.grid = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  xlab("CpG sites") +
  ylab("Weighted absolute PhenoAge replicate difference")  +
  scale_y_continuous(limits = c(0, 2.5))


ordering <- levels(with(anova_data, reorder(probe, rep1_rep2)))
# plot the weights
p2 <- ggplot(data = data.frame(probe = scDNAmClock:::pheno_age_dat$CpG, weight = abs(scDNAmClock:::pheno_age_dat$weight)) %>%
           left_join(probe_ICC, by = "probe") %>% 
             mutate(probe = factor(probe, levels = ordering)) ) +
  geom_point(aes(x = probe, y = weight, color = ICC), size = 1, alpha = 0.5) +
  theme_classic() +
  coord_flip() +
  theme(panel.grid = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(0.9, 0.2)) +
  xlab("CpG sites") +
  ylab("Absolute weight") +
  scale_color_viridis_c()
p <- cowplot::plot_grid(p1, p2, nrow = 1, align = "h", rel_widths = c(1, 1))
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/probe_contribution.png", width = 7.5, height = 5)

# preprcessing for the full dataset
data <- readRDS("/gpfs/ysm/project/mw957/data/processed/Lehne2015/Lehne2015_beta.rds") # 473864
# filter for probes with sd > 0.0042 (based on minimum sd of the Levine PhenoAge probe)
data_sd <- apply(data[,-1], 1, function(x) sd(x, na.rm = T))
data <- data[data_sd > 0.0042,] # remove probes with really low standard deviation
# 472222

probe_missing <- apply(data[,-1], 1, function(x) sum(is.na(x)))
data <- data[probe_missing < 0.5 * (ncol(data) - 1),] # remove probe with > 50% missing 
# 472009

#patient_missing <- apply(data[,-1], 2, function(x) sum(is.na(x))) # examine patient missingness OK

# knn imputation 
data[,2:ncol(data)] <- impute::impute.knn(data[,2:ncol(data)])$data
saveRDS(data, file = "/gpfs/ysm/project/mw957/data/processed/Lehne2015/Lehne2015_beta_imputed.rds")

# extract the Levine clock data
cpg_clocks <- readr::read_csv("/gpfs/ysm/project/mw957/repos/scDNAmClock/data-raw/cpgAllonetime.cpgInfo.csv")
Levine_clock  <- cpg_clocks %>%
  filter(`Levine-DNAm-PhenoAge` == "Levine") 
Levine_data <- data %>%
  filter(ID_REF %in% Levine_clock$prob)
saveRDS(Levine_data, file = "/gpfs/ysm/project/mw957/data/processed/Lehne2015/Lehne2015_beta_Levine.rds")

# apply variance stabilizing transformation
M_data_scaled <- asinh(as.matrix(M_data[,-1]))
replicate_1_names <- filter(samples, group == 1) %>%
  .$name
replicate_2_names <- filter(samples, group == 2) %>%
  .$name
other_names <- setdiff(colnames(M_data)[-1], samples$name)

# rep_1_data <- as.matrix(M_data_scaled[,c(replicate_1_names, other_names)])
# rownames(rep_1_data) <- M_data$ID_REF
# rep_2_data <- as.matrix(M_data_scaled[,c(replicate_2_names, other_names)])
# rownames(rep_2_data) <- M_data$ID_REF
# 
# a <- 3
# 
# rep_1_data <- t(apply(rep_1_data, 1, function(x) truncate_outlier(x, a = a)))
# rep_2_data <- t(apply(rep_2_data, 1, function(x) truncate_outlier(x, a = a)))

rep_data <- as.matrix(M_data_scaled)
rownames(rep_data) <- M_data$ID_REF

set.seed(1)
# replicate_1_svd <- rsvd::rsvd(rep_1_data, k = min(dim(rep_1_data)), q = 2)
# replicate_2_svd <- rsvd::rsvd(rep_2_data, k = min(dim(rep_2_data)), q = 2)
rep_svd <- rsvd::rsvd(rep_data, k = min(dim(rep_data)), q = 2)

# estimate the noise standard deviation 
noises <- rep_data[,replicate_1_names] - rep_data[,replicate_2_names]
dim(noises) <- NULL
#library(ggpubr)
#p <- ggqqplot(noises)
#ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/qqplot.png", width = 5, height = 4)
# expression <- (rep_1_data[,replicate_1_names] + rep_2_data[,replicate_2_names])/2
# dim(expression) <- NULL
# p <- ggplot(data = data.frame(expression = expression, noise = noises)) +
#   geom_point(aes(x = expression, y = noise), color = "steelblue", alpha = 0.4, size = 0.1) +
#   theme_classic() +
#   theme(panel.grid = element_blank()) +
#   xlab(expression("arcinsh("~bar(hat(M)[ij])~")")) +
#   ylab("arcinsh("~hat(M)[ij1]~")-arcinsh("~hat(M)[ij2]~")")
# ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/heteroscedasticity.png", width = 5, height = 4)

library(MASS)
t_fit <- fitdistr(unlist(noises), "t", start = list(m = mean(unlist(noises)), s = sd(unlist(noises)), df = 39))
normal_fit <- fitdistr(unlist(noises), "normal") # normal_fit

 p <- ggplot(data = data.frame(noise = unlist(noises)), aes(x = noise)) +
   geom_histogram(bins = 120, aes(y=..density..), fill = "steelblue") +
   theme_classic() +
   theme(panel.grid = element_blank()) +
  #stat_function(fun = metRology::dt.scaled, args = list(mean = normal_fit$estimate[1], sd = normal_fit$estimate[2], df = normal_fit$estimate[3]), n = 500, linetype = "dashed", size = 0.3) +
   stat_function(fun = dnorm, args = list(mean = 0, sd = normal_fit$estimate[2])) +
   xlab("arcsinh("~hat(M)[ij1]~")-arcsinh("~hat(M)[ij2]~")")
# ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/density.png", width = 5, height = 4)

n_lambda <- 50
lambda_max <- 100
tau <- normal_fit$estimate[2]
n_tau <- length(tau)
lambda <- matrix(NA, nrow = n_tau, ncol = n_lambda)
for (i in seq_along(tau)) {
  lambda[i,] <- seq(0, tau[i] * lambda_max, length.out = n_lambda)
}

# rep1_SURE <- rep2_SURE <- matrix(NA, nrow = n_tau, ncol = n_lambda)
# 
# for (i in seq(n_tau)) {
#   for (j in seq(n_lambda)) {
#     rep1_SURE[i, j] <- sure_svt(lambda[i,j], tau[i], rep_1_data, s = replicate_1_svd$d, is_real = T, svThreshold = 1e-8)
#   }
# }
# 
# for (i in seq(n_tau)) {
#   for (j in seq(n_lambda)) {
#     rep2_SURE[i, j] <- sure_svt(lambda[i,j], tau[i], rep_2_data, s = replicate_2_svd$d, is_real = T, svThreshold = 1e-8)
#   }
# }

rep_SURE <- matrix(NA, nrow = n_tau, ncol = n_lambda)

for (i in seq(n_tau)) {
  for (j in seq(n_lambda)) {
    rep_SURE[i, j] <- sure_svt(lambda[i,j], tau[i], rep_data, s = rep_svd$d, is_real = T, svThreshold = 1e-8)
  }
}

p <- plot_SURE(lambda[i,], rep_SURE[i,]) + ylab("SURE") 

# p1 <- plot_SURE(lambda[i,], rep1_SURE[i,]) + ylab("Replicate 1 SURE")
# p2 <- plot_SURE(lambda[i,], rep2_SURE[i,]) + ylab("Replicate 2 SURE")
# p <- cowplot::plot_grid(p1, p2, nrow = 1, align = "h")
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/replicate_lambda_2.png", width = 3, height = 2)

# reconstructed_1 <- SVT_denoise(rep_1_data, lambda = lambda[which.min(rep1_SURE)], svd = replicate_1_svd)
# reconstructed_2 <- SVT_denoise(rep_2_data, lambda = lambda[which.min(rep2_SURE)], svd = replicate_2_svd)
# truncated_reconstructed <- as.data.frame(cbind(reconstructed_1[, replicate_1_names], 
#                                                reconstructed_2[, replicate_2_names]))
reconstructed <- SVT_denoise(rep_data, lambda = lambda[which.min(rep_SURE)], svd = rep_svd)
#reconstructed <- apply(reconstructed, 2, function(x) truncate_outlier(x, a = 1.5))
truncated_reconstructed <- as.data.frame(cbind(reconstructed[, replicate_1_names], reconstructed[, replicate_2_names]))
truncated_reconstructed$probe <- all_data$ID_REF
truncated_reconstructed_mat <- apply(truncated_reconstructed[,-ncol(truncated_reconstructed)], 2, function(x) 2^sinh(x)/(1 + 2^sinh(x)))
rownames(truncated_reconstructed_mat) <- all_data$ID_REF
pheno_age_2 <- scDNAmClock::PhenoAge(truncated_reconstructed_mat)
plot_data_2 <- data.frame(name = names(pheno_age_2$y), 
                        pheno_age = unname(pheno_age_2$y)) %>%
  left_join(samples, by = "name") %>%
  dplyr::select(-c(name, geo_accession)) %>%
  tidyr::spread(group, pheno_age) %>%
  mutate(diff = `1` - `2`, abs_diff = abs(`1` - `2`)) 

combined_data <- cbind(plot_data_2[, 4:6], plot_data[, 4:6])
colnames(combined_data)[1:3] <- c("SVD_1", "SVD_2", "SVD_diff")

cols <- c("SVT" = "steelblue", 
          "Original" = "gray")

p <- ggplot(data = combined_data) +
  geom_abline(alpha = 0.6, linetype = "dashed") +
  geom_abline(intercept = 2, alpha = 0.3, linetype = "dashed") +
  geom_abline(intercept = -2, alpha = 0.3, linetype = "dashed") +
  geom_point(aes(x = `1`, y = `2`, color = "Original"), size = 1) +
  geom_point(aes(x = SVD_1, y = SVD_2, color = "SVT"), size = 2) +
  geom_segment(aes(x = `1`, y = `2`, xend = SVD_1, yend = SVD_2), 
               arrow = arrow(length = unit(0.01, "npc")), alpha = 0.3) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(name = "", values = cols) +
  theme(legend.position = c(0.7, 0.2)) +
  xlab(expression("Replicate 1"~beta)) +
  ylab(expression("Replicate 2"~beta)) 
ggsave(p, file = paste0("~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/phenoage_svd_truncated_2.pdf"), width = 4, height = 3)  

cols <- c("SVT" = "steelblue", 
          "Original" = "gray")

p <- ggplot(data = combined_data) +
  geom_histogram(aes(x = abs(diff), fill = "Original"), bins = 10, alpha = 0.8) +
  geom_histogram(aes(x = abs(SVD_diff), fill = "SVT"), bins = 10, alpha = 0.8) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8, 0.8)) +
  xlab("Absolute PhenoAge replicate difference") +
  scale_fill_manual(name = "", values = cols)
  
ggsave(p, file = paste0("~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/phenoage_hist_2.pdf"), width = 4, height = 3)
dat <- combined_data %>%
  tibble::rownames_to_column() %>%
  select(rowname, SVD_diff, diff) %>%
  tidyr::gather(method, value, -rowname) %>%
  mutate(value = abs(value))
dat$method[dat$method == "SVD_diff"] <- "SVD"
dat$method[dat$method == "diff"] <- "original"

p <- ggplot(data = dat, aes(x = method, y = value, color = method)) +
  geom_hline(aes(yintercept = 2), linetype = 'dashed') +
  geom_boxplot(width = 0.2, alpha = 0.2) +
  geom_point(aes(group = rowname), size = 0.5) +
  geom_line(aes(group = rowname), size = 0.2) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.8, 0.8)) +
  ylab("Absolute replicate PhenoAge difference") +
  ggtitle("Paired t-test: t = -4.73, df = 35, p = 3.6e-5") +
  scale_color_manual(values = c("gray", "steelblue"))
ggsave(p, file = paste0("~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/phenoage_boxplot_2.pdf"), width = 4, height = 3)

t.test(dat$value[1:36], dat$value[37:72], paired = T)

sum(dat$value[1:36] > 2)

# examine variability contribution
summary_data_2 <- as.data.frame(pheno_age_2$Ax) %>%
  mutate(probe = rownames(pheno_age_2$Ax)) %>%
  tidyr::gather(name, value, -probe) %>%
  left_join(samples, by = "name") %>%
  dplyr::select(-c(name, geo_accession, gender, age)) %>%
  tidyr::spread(group, value) %>%
  mutate(abs_deviation = abs(`1` - `2`)) %>%
  group_by(probe) %>%
  summarize(abs_deviation = mean(abs_deviation)) %>%
  left_join(data.frame(probe = scDNAmClock:::pheno_age_dat$CpG,
                       weight = scDNAmClock:::pheno_age_dat$weight), by = "probe") %>%
  left_join(probe_ICC, by = "probe")

summary_data_combined <- summary_data %>%
  mutate(truncated_abs_deviation = summary_data_2$abs_deviation)

cols <- c("SVT" = "circle",
          "Original" = "square")
p <- ggplot(data = summary_data_combined) +
  geom_point(aes(x = weight, y = abs_deviation, shape = "Original"), color = "gray", alpha = 0.6) +
  geom_point(aes(x = weight, y = truncated_abs_deviation, color = ICC, shape = "SVT"), alpha = 0.9) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ylab("weight * mean replicate difference") +
  scale_color_viridis_c() +
  xlab("weight")  +
  scale_shape_manual(name = "", values = cols)
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/probe_variation_truncated.pdf", width = 6, height = 4)

p <- ggplot(data = summary_data_combined,
            aes(x = category)) +
  geom_boxplot(aes(color = ICC, y = truncated_abs_deviation)) +
  geom_jitter(aes(color = ICC, y = truncated_abs_deviation), alpha = 0.9) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.1, 0.8)) +
  ylab("weight * mean replicate difference") +
  scale_color_viridis_c() 

p <- ggplot(data = summary_data_combined) +
  geom_abline(size = 0.2) +
  geom_point(aes(x = abs_deviation, y = truncated_abs_deviation, color = ICC), size = 0.3, alpha = 0.8) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.7)) +
  facet_wrap(~category, ncol = 2) +
  scale_color_viridis_c() +
  xlab("Absolute PhenoAge difference (original)") +
  ylab("Absolute PhenoAge difference (SVT)")
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/probe_variation_box_2.pdf", width = 4, height = 3)


# examine the contribution of probes to deviation 
anova_data_2 <- as.data.frame(pheno_age_2$Ax) %>%
  mutate(probe = rownames(pheno_age_2$Ax)) %>%
  tidyr::gather(name, Ax, -probe) %>%
  dplyr::left_join(dplyr::select(samples, c(name, group, sample))) %>%
  dplyr::select(-name) %>% 
  tidyr::spread(group, Ax) %>%
  mutate(`rep1_rep2` = abs(`1` - `2`))  %>%
  mutate(probe = factor(probe, levels = ordering))

beta_noise <- as.matrix(all_data[,replicate_1_names] - all_data[,replicate_2_names] )
beta_noise_after <- truncated_reconstructed_mat[,replicate_1_names] - truncated_reconstructed_mat[, replicate_2_names]
dim(beta_noise_after) <- NULL
dim(beta_noise) <- NULL
p_1 <- ggplot(data = data.frame(noise = beta_noise, noise_after = beta_noise_after)) +
  geom_abline() +
  geom_hline(aes(yintercept = 0)) +
  geom_point(aes(x = noise, y = noise_after), color = "steelblue", alpha = 0.3, size = 0.2) +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  xlab(hat(beta)[ij1]~"-"~hat(beta)[ij2]~"(original)") +
  ylab(tilde(beta)[ij1]~"-"~tilde(beta)[ij2]~"(SVT)") +
  scale_x_continuous(limits = c(-0.3, 0.3)) +
  scale_y_continuous(limits = c(-0.3, 0.3))

p1_2 <- ggplot(data = anova_data_2) +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed") +
  geom_point(aes(x = probe, y = rep1_rep2), size = 0.1, alpha = 0.5,
             color = "steelblue") +
  theme_classic() +
  coord_flip() +
  theme(panel.grid = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  xlab("CpG sites") +
  ylab("Weighted absolute PhenoAge replicate difference") +
  scale_y_continuous(limits = c(0, 2.5))


p <- cowplot::plot_grid(p_1, p1_2, nrow = 1, align = "h", rel_widths = c(1, 1))
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/probe_contribution_2.png", width = 7.5, height = 5)

# patients who have PhenoAge difference greater than 2 years
patient_diff <- data.frame(y = pheno_age_2$y) %>%
  mutate(name = colnames(pheno_age_2$Ax)) %>%
  dplyr::left_join(dplyr::select(samples, c(name, group, sample))) %>%
  dplyr::select(-name) %>% 
  tidyr::spread(group, y) %>%
  mutate(`rep1_rep2` = abs(`1` - `2`)) %>%
  arrange(desc(rep1_rep2))# 15 patients


all_probes <-  anova_data_2 %>%
  dplyr::select(-c(`1`, `2`)) %>%
  group_by(probe) %>%
  mutate(mean = mean(rep1_rep2)) %>%
  mutate(sample = factor(sample, levels = patient_diff$sample)) %>%
  tidyr::spread(sample, rep1_rep2) %>%
  arrange(desc(mean)) 
  
all_probes_mat <- as.matrix(all_probes[,-c(1,2)])

ha <- ComplexHeatmap::HeatmapAnnotation(`PhenoAge difference` = patient_diff$rep1_rep2)
ht1 <- ComplexHeatmap::Heatmap(all_probes_mat, name = "original", top_annotation = ha, show_row_names = F, cluster_rows = F, cluster_columns = F)
png("~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/probe_contribution_3.png", width = 500, height = 500)
ht1
dev.off()



# Examine the outliers
outliers <- as.data.frame(M_data_scaled[,replicate_1_names] - M_data_scaled[,replicate_2_names]) %>%
  mutate(probe = scDNAmClock:::pheno_age_dat$CpG) %>%
  tidyr::gather(name, noise, -probe) %>%
  mutate(gt.05 = noise > 0.5 | noise < -0.5) %>%
  group_by(probe) %>%
  summarize(outlier_counts = sum(gt.05))
p <- ggplot(data = as.data.frame(outliers)) +
  geom_histogram(aes(x = outlier_counts)) +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  xlab("outlier counts")
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/outlier_probes.png", width = 5, height = 4)

# Can these outliers be corrected by truncation? 
# plot of noise and plot of actual values
dat <- data.frame(
  rep_1 = unlist(M_data_scaled[which(scDNAmClock:::pheno_age_dat$CpG %in% "cg17133388"), replicate_1_names]),
  rep_2 = unlist(M_data_scaled[which(scDNAmClock:::pheno_age_dat$CpG %in% "cg17133388"), replicate_2_names])
) %>%
  mutate(noise = rep_1 - rep_2) %>%
  mutate(outlier = abs(rep_1 - median(rep_1)) > 1.5 * IQR(rep_1) | abs(rep_2 - median(rep_2)) > 1.5 * IQR(rep_2) ) 
dat$trunc_rep_1 <- truncate_outlier(dat$rep_1)
dat$trunc_rep_2 <- truncate_outlier(dat$rep_2)

# shrink outliers
NA_outlier <- function(x, a = 1.5) {
  iqr <- IQR(x)
  median <- median(x)
  x[x > median + a * iqr] <- NA#median + a * iqr
  x[x < median - a * iqr] <- NA#median - a * iqr
  x
}

truncate_outlier <- function(x, a = 1.5) {
  iqr <- IQR(x)
  median <- median(x)
  x[x > median + a * iqr] <- median + a * iqr
  x[x < median - a * iqr] <- median - a * iqr
  x
}

p1 <- ggplot(data = dat, aes(x = rep_1, y = rep_2), color = "gray") +
  geom_abline(size = 0.1) +
  geom_abline(aes(slope = 1, intercept = 0.5), size = 0.1) +
  geom_abline(aes(slope = 1, intercept = -0.5), size = 0.1) +
  geom_segment(aes(xend = trunc_rep_1, yend = trunc_rep_2), size = 0.2) +
  geom_point() +
  geom_point(aes(x = trunc_rep_1, y = trunc_rep_2), color = "steelblue") +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  xlab("outlier counts") +
  geom_rug(col=rgb(.5,0,0,alpha=.2))

ggsave(p1, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/outlier_probes_ex.png", width = 5, height = 4)



# replicate data visualization
original_dat <- all_data[,c(replicate_1_names,replicate_2_names)]
SVT_dat <- truncated_reconstructed_mat[,c(replicate_1_names,replicate_2_names)]
labels <- c(seq(36), seq(36))
SVT_normed <- t(apply(SVT_dat, 1, scale))
normed <- t(apply(original_dat, 1, scale))

ha <- ComplexHeatmap::HeatmapAnnotation(sample = as.character(labels),
                                        show_legend = c(sample = F))
ht1 <- ComplexHeatmap::Heatmap(normed, name = "original", col = colorRampPalette(c("#0047BB", "white", "#D1350F"))(10), top_annotation = ha, show_row_names = F)
ht2 <- ComplexHeatmap::Heatmap(SVT_normed, name = "SVT", col = colorRampPalette(c("#0047BB", "white", "#D1350F"))(10), top_annotation = ha, show_row_names = F)

png("~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/replicate_heatmap.png", width = 600, height = 500)
ht1 + ht2
dev.off()




# applying capping on the dataset
cap <- function(data, map_list) {
  purrr::imap(data, function(x, y) {
      trunc_list <- map_list[[y]]
      x[x < trunc_list[1]] <- trunc_list[1]
      x[x > trunc_list[2]] <- trunc_list[2]
      x
    }) %>%
    dplyr::bind_rows()
}

quantile_capping <- readRDS("data-raw/quantile.rds")
rep_1_data_capped <- cap(as.data.frame(t(rep_1_data)), quantile_capping)
rep_2_data_capped <- cap(as.data.frame(t(rep_2_data)), quantile_capping)
rep_1_data_capped <- t(as.matrix(rep_1_data_capped))
colnames(rep_1_data_capped) <- colnames(rep_1_data)
rep_2_data_capped <- t(as.matrix(rep_2_data_capped))
colnames(rep_2_data_capped) <- colnames(rep_2_data)


i = 2

par(mfrow = c(1,2))
hist(unlist(rep_1_data[i,]), main = "original", xlab = "original beta")
hist(unlist(rep_1_data_capped[,i]), main = "truncated", xlab = "truncated beta")



set.seed(100)
rep_1_svd_3 <- rsvd::rsvd(rep_1_data_capped, k = min(dim(rep_1_data_capped)), q = 2)
rep_2_svd_3 <- rsvd::rsvd(rep_2_data_capped, k = min(dim(rep_2_data_capped)), q = 2)



n_lambda <- 50
lambda_max <- 60
tau <- sd(rep_1_data_capped[,replicate_1_names] - rep_2_data_capped[,replicate_2_names])/sqrt(2)
n_tau <- length(tau)
lambda <- matrix(NA, nrow = n_tau, ncol = n_lambda)
for (i in seq_along(tau)) {
  lambda[i,] <- seq(0, tau[i] * n_lambda, length.out = n_lambda)
}

rep1_SURE <- rep2_SURE <- matrix(NA, nrow = n_tau, ncol = n_lambda)

for (i in seq(n_tau)) {
  for (j in seq(n_lambda)) {
    rep1_SURE[i, j] <- sure_svt(lambda[i,j], tau[i], rep_1_data_capped, s = rep_1_svd_3$d, is_real = T, svThreshold = 1e-8)
  }
}

for (i in seq(n_tau)) {
  for (j in seq(n_lambda)) {
    rep2_SURE[i, j] <- sure_svt(lambda[i,j], tau[i], rep_2_data_capped, s = rep_2_svd_3$d, is_real = T, svThreshold = 1e-8)
  }
}

i <- 1
p1 <- plot_SURE(lambda[i,], rep1_SURE[i,]) + ylab("Replicate 1 SURE")
p2 <- plot_SURE(lambda[i,], rep2_SURE[i,]) + ylab("Replicate 2 SURE")
p <- cowplot::plot_grid(p1, p2, nrow = 1, align = "h")

reconstructed_1 <- SVT_denoise(rep_1_data_capped, lambda = lambda[which.min(rep1_SURE)], svd = rep_1_svd_3)
reconstructed_2 <- SVT_denoise(rep_2_data_capped, lambda = lambda[which.min(rep2_SURE)], svd = rep_2_svd_3)

# truncated_reconstructed <- as.data.frame(cbind(reconstructed_1[, replicate_1_names], 
#                                                reconstructed_2[, replicate_2_names]))
# truncated_reconstructed_mat <- as.matrix(truncated_reconstructed)
# rownames(truncated_reconstructed_mat) <- all_data$ID_REF
# truncated_reconstructed$probe <- all_data$ID_REF

truncated_reconstructed <- as.data.frame(cbind(rep_1_data_capped[, replicate_1_names], 
                                               rep_2_data_capped[, replicate_2_names]))
truncated_reconstructed_mat <- as.matrix(truncated_reconstructed)
rownames(truncated_reconstructed_mat) <- all_data$ID_REF
truncated_reconstructed$probe <- all_data$ID_REF

pheno_age_3 <- scDNAmClock::PhenoAge(truncated_reconstructed_mat)

plot_data_3 <- data.frame(name = names(pheno_age_3$y), 
                          pheno_age = unname(pheno_age_3$y)) %>%
  left_join(samples, by = "name") %>%
  select(-c(name, geo_accession)) %>%
  tidyr::spread(group, pheno_age) %>%
  mutate(diff = `1` - `2`, abs_diff = abs(`1` - `2`)) 

combined_data <- cbind(plot_data_3[, 4:6], plot_data[, 4:6])
colnames(combined_data)[1:3] <- c("SVD_1", "SVD_2", "SVD_diff")

cols <- c("SVT" = "steelblue", 
          "Original" = "gray")

p <- ggplot(data = combined_data) +
  geom_abline(alpha = 0.6, linetype = "dashed") +
  geom_abline(intercept = 5, alpha = 0.3, linetype = "dashed") +
  geom_abline(intercept = -5, alpha = 0.3, linetype = "dashed") +
  
  geom_point(aes(x = `1`, y = `2`, color = "Original"), size = 1) +
  geom_point(aes(x = SVD_1, y = SVD_2, color = "SVT"), size = 2) +
  geom_segment(aes(x = `1`, y = `2`, xend = SVD_1, yend = SVD_2), 
               arrow = arrow(length = unit(0.01, "npc")), alpha = 0.3) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(name = "", values = cols) +
  theme(legend.position = c(0.7, 0.2)) +
  xlab("Replicate 1") +
  ylab("Replicate 2") 
p
#ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/phenoage_svd_truncated_2.pdf", width = 4, height = 3)  

cols <- c("SVT" = "steelblue", 
          "Original" = "gray")

p <- ggplot(data = combined_data) +
  geom_histogram(aes(x = abs(diff), fill = "Original"), bins = 10, alpha = 0.8) +
  geom_histogram(aes(x = abs(SVD_diff), fill = "SVT"), bins = 10, alpha = 0.8) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8, 0.8)) +
  xlab("Absolute PhenoAge replicate difference") +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  scale_fill_manual(name = "", values = cols)

ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/phenoage_hist_2.pdf", width = 4, height = 3)  

