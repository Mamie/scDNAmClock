# evaluate the age prediction using the phenoage model
library(scDNAmClock)
samples <- readr::read_csv("data-raw/technical_rep_sample_lists.csv")
data <- readr::read_tsv("data-raw/tech_rep_data.tsv")
colnames(data) <- c("probe", samples$name)

# Let's fill out the matrix by nearest neighbor imputation
data_mat <- as.matrix(data[,-1])
data_mat <- impute::impute.knn(data_mat)$data
rownames(data_mat) <- data$probe
data[,2:ncol(data)] <- data_mat

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
res <- data.frame(
  probe = names(probes_data),
  ICC = purrr::map_dbl(probe_ICC, ~.x$value)
)
res$ICC[res$ICC < 0] <- 0



# phenoage prediction
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
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/phenoage.pdf", width = 6, height = 5)  
  

p <- ggplot(data = plot_data) +
  geom_histogram(aes(x = abs_diff), fill = "steelblue", bins = 10, alpha = 0.8) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("Absolute PhenoAge replicate difference") +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/phenoage_hist.pdf", width = 4, height = 2)  
res <- t(as.matrix(summary(plot_data$abs_diff)))
kableExtra::kable(res, digits = 2) %>% 
  kableExtra::kable_styling()

# examine variability contribution
View(PhenoAge)
summary_data <- as.data.frame(pheno_age$Ax) %>%
  mutate(probe = rownames(pheno_age$Ax)) %>%
  tidyr::gather(name, value, -probe) %>%
  left_join(samples, by = "name") %>%
  select(-c(name, geo_accession, gender, age)) %>%
  tidyr::spread(group, value) %>%
  mutate(abs_deviation = abs(`1` - `2`)) %>%
  group_by(probe) %>%
  summarize(abs_deviation = mean(abs_deviation)) %>%
  left_join(data.frame(probe = pheno_age_dat$CpG,
                       weight = pheno_age_dat$weight), by = "probe") %>%
  left_join(res, by = "probe")

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

res$category <- "excellent (> 0.9)"
res$category[res$ICC < 0.9] <- "good (0.75 - 0.9)"
res$category[res$ICC < 0.75] <- "moderate (0.5 - 0.75)"
res$category[res$ICC < 0.5] <- "poor (< 0.5)"

summary_data <- left_join(summary_data, res) %>%
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

head(summary_data)

summary_data$abs_dev_div_weight <- summary_data$abs_deviation/abs(summary_data$weight)
summary_data %>%
  arrange(desc(abs_dev_div_weight))

# use rvsd to denoise the dataset
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

