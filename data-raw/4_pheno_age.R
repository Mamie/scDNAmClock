# evaluate the age prediction using the phenoage model
library(scDNAmClock)
library(dplyr)
library(ggplot2)

samples <- readr::read_csv("data-raw/technical_rep_sample_lists.csv")
all_data <- readRDS("data-raw/Lehne2015_beta_Levine.rds")
data <- all_data[, c("ID_REF", samples$name)]
colnames(data)[1] <- "probe"

# Let's fill out the matrix by nearest neighbor imputation
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
  geom_abline(intercept = 5, alpha = 0.3, linetype = "dashed") +
  geom_abline(intercept = -5, alpha = 0.3, linetype = "dashed") +
  #geom_segment(aes(x = `1`, y = `2`, xend = age, yend = age), size = 0.06) +
  geom_point(aes(x = `1`, y = `2`, color = "DNAm PhenoAge"), size = 3) +
  #geom_point(aes(x = age, y = age, color = "Actual age"), size = 2) +
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

# check effect of doing SVD on using more patient samples 
replicate_1_names <- filter(samples, group == 1) %>%
  .$name
replicate_2_names <- filter(samples, group == 2) %>%
  .$name
other_names <- setdiff(colnames(all_data)[-1], samples$name)
k = round(nrow(all_data) * 0.3, 0)
set.seed(200)
replicate_1_k <- choose_k(all_data[,c(replicate_1_names, other_names)], K = k, pval_thresh = 1e-10, noise_start = k * 0.8)
replicate_2_k <- choose_k(all_data[,c(replicate_2_names, other_names)], K = k, pval_thresh = 1e-10, noise_start = k * 0.8)

replicate_1_svd <- rsvd::rsvd(all_data[,c(replicate_1_names, other_names)], k = replicate_1_k$k, q = 2)
replicate_2_svd <- rsvd::rsvd(all_data[,c(replicate_2_names, other_names)], k = replicate_2_k$k, q = 2)

reconstructed_1 <- replicate_1_svd$u %*% diag(replicate_1_svd$d) %*% t(replicate_1_svd$v)
reconstructed_2 <- replicate_2_svd$u %*% diag(replicate_2_svd$d) %*% t(replicate_2_svd$v)

colnames(reconstructed_1) <- c(replicate_1_names, other_names)
colnames(reconstructed_2) <- c(replicate_2_names, other_names)

truncated_reconstructed <- as.data.frame(cbind(reconstructed_1[, replicate_1_names], 
                                 reconstructed_2[, replicate_2_names]))
truncated_reconstructed_mat <- as.matrix(truncated_reconstructed)
rownames(truncated_reconstructed_mat) <- all_data$ID_REF
truncated_reconstructed$probe <- all_data$ID_REF

pheno_age_2 <- scDNAmClock::PhenoAge(truncated_reconstructed_mat)

plot_data_2 <- data.frame(name = names(pheno_age_2$y), 
                        pheno_age = unname(pheno_age_2$y)) %>%
  left_join(samples, by = "name") %>%
  select(-c(name, geo_accession)) %>%
  tidyr::spread(group, pheno_age) %>%
  mutate(diff = `1` - `2`, abs_diff = abs(`1` - `2`)) 

combined_data <- cbind(plot_data_2[, 4:6], plot_data[, 4:6])
colnames(combined_data)[1:3] <- c("SVD_1", "SVD_2", "SVD_diff")

cols <- c("SVD denoised" = "steelblue", 
          "Original" = "gray")

p <- ggplot(data = combined_data) +
  geom_abline(alpha = 0.6, linetype = "dashed") +
  geom_abline(intercept = 5, alpha = 0.3, linetype = "dashed") +
  geom_abline(intercept = -5, alpha = 0.3, linetype = "dashed") +
  
  geom_point(aes(x = `1`, y = `2`, color = "Original"), size = 1) +
  geom_point(aes(x = SVD_1, y = SVD_2, color = "SVD denoised"), size = 2) +
  geom_segment(aes(x = `1`, y = `2`, xend = SVD_1, yend = SVD_2), 
               arrow = arrow(length = unit(0.01, "npc")), alpha = 0.3) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(name = "", values = cols) +
  theme(legend.position = c(0.7, 0.2)) +
  xlab("Replicate 1") +
  ylab("Replicate 2") 
p
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/phenoage_svd_truncated.pdf", width = 6, height = 5)  

cols <- c("SVD denoised" = "steelblue", 
          "Original" = "gray")

p <- ggplot(data = combined_data) +
  geom_histogram(aes(x = abs(diff), fill = "Original"), bins = 10, alpha = 0.8) +
  geom_histogram(aes(x = abs(SVD_diff), fill = "SVD denoised"), bins = 10, alpha = 0.8) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8, 0.8)) +
  xlab("Absolute PhenoAge replicate difference") +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  scale_fill_manual(name = "", values = cols)
  
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/phenoage_hist.pdf", width = 6, height = 3)  


# examine variability contribution
summary_data_2 <- as.data.frame(pheno_age_2$Ax) %>%
  mutate(probe = rownames(pheno_age_2$Ax)) %>%
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

summary_data_combined <- summary_data %>%
  mutate(truncated_abs_deviation = summary_data_2$abs_deviation)


cols <- c("SVD denoised" = "circle",
          "Original" = "square")
p <- ggplot(data = summary_data_combined) +
  geom_point(aes(x = weight, y = abs_deviation, shape = "Original"), color = "gray", alpha = 0.6) +
  geom_point(aes(x = weight, y = truncated_abs_deviation, color = ICC, shape = "SVD denoised"), alpha = 0.9) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.9, 0.4)) +
  ylab("weight * mean replicate difference") +
  scale_color_viridis_c() +
  xlab("weight")  +
  scale_shape_manual(name = "", values = cols)
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/probe_variation_truncated.pdf", width = 7, height = 4.5)

p <- ggplot(data = summary_data_combined,
            aes(x = category, y = truncated_abs_deviation, color = ICC)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.9) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.1, 0.8)) +
  ylab("weight * mean replicate difference") +
  scale_color_viridis_c()
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/probe_variation_box_2.pdf", width = 7, height = 4.5)

