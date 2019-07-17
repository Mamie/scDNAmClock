library(dplyr)
library(ggplot2)

metadata <- readr::read_tsv( "data-raw/liver_samples.tsv")

# evaluate the missingness from the dataset
inner_joined <- readRDS("~/Downloads/mm9_mm10_inner_join.rds")
dim(inner_join) # 12689082      307
# 152 samples before filtering

# compute the methylation percentage
as.data.frame(inner_joined) %>%
  tidyr::gather(category, value, -c(chr, position, strand)) %>%
  tidyr::separate(category, c("type", "geo_accession"), sep = "[.]") %>%
  tidyr::spread(type, value) %>%
  mutate(`% methylated` = met_reads/coverage) %>%
  select(chr, position, geo_accession, `% methylated`) %>%
  tidyr::spread(geo_accession, `% methylated`)

data <- readRDS("~/Downloads/mm9_mm10_inner_join_methyl_perc.rds")


# filter sites until missingness is around 2 %

missing_proportion <- function(x) {
  mean(is.na(x[,-c(1,2)]))
}
missing_proportion(data) # 71.56 %

# evaluate sites missingness
site_missingness <- function(x) {
  apply(x[,-c(1,2)], 1, function (t) mean(is.na(t)))
}

patient_missingness <- function(x) {
  apply(x[,-c(1,2)], 2, function(t) mean(is.na(t)))
}

miss_sites <- site_missingness(data)
miss_patients <- patient_missingness(data)


p <- data.frame(geo_accession = colnames(data)[-c(1,2)], 
           missing_proportion = miss_patients) %>%
  left_join(metadata) %>%
  select(missing_proportion, dataset) %>%
  ggplot(data = .) +
  geom_histogram(aes(x = missing_proportion, fill = dataset), bins = 40) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/data_integration/figs/missing_proportion.pdf", width = 5, height = 3)

cutoffs <- seq(0.1, 0.9, by = 0.1)
sites_avail <- c()
missingness <- c()
for (cutoff in cutoffs) {
  data2 <- data[miss_sites < cutoff, ]
  sites_avail <- c(sites_avail, nrow(data2))
  missingness <- c(missingness, missing_proportion(data2))
}
dat <- data.frame(cutoff = cutoffs, sites_avail = sites_avail,
                  missingness = missingness)

p <- ggplot(data = dat, aes(x = missingness, y = sites_avail)) +
  geom_text(aes(label = cutoff)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_log10(labels = scDNAmClock:::fancy_scientific) +
  ylab("# sites available") +
  xlab("missing proportion")
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/data_integration/figs/missing_cutoff.png", width = 4, height = 3)

# choose 0.6 as cutoff for iteration 1
data2 <- data[miss_sites < 0.6, ]
miss_patients <- patient_missingness(data2)

# discard patients
cutoffs <- seq(0.1, 0.9, by = 0.1)
sites_avail <- c()
missingness <- c()
for (cutoff in cutoffs) {
  data3 <- data2[, c(T, T, miss_patients < cutoff)]
  sites_avail <- c(sites_avail, ncol(data3) - 2)
  missingness <- c(missingness, missing_proportion(data3))
}
dat <- data.frame(cutoff = cutoffs, sites_avail = sites_avail,
                  missingness = missingness)
p <- ggplot(data = dat, aes(x = missingness, y = sites_avail)) +
  geom_text(aes(label = cutoff)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_log10(labels = scDNAmClock:::fancy_scientific) +
  ylab("# patients available") +
  xlab("missing proportion")
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/data_integration/figs/missing_cutoff_2.png", width = 4, height = 3)

# choose 0.7
data4 <- data3[, c(T, T, miss_patients < 0.7)]

miss_sites <- site_missingness(data4)
cutoffs <- seq(0.1, 0.9, by = 0.1)
sites_avail <- c()
missingness <- c()
for (cutoff in cutoffs) {
  data5 <- data4[miss_sites < cutoff, ]
  sites_avail <- c(sites_avail, nrow(data5))
  missingness <- c(missingness, missing_proportion(data5))
}
dat <- data.frame(cutoff = cutoffs, sites_avail = sites_avail,
                  missingness = missingness)
p <- ggplot(data = dat, aes(x = missingness, y = sites_avail)) +
  geom_text(aes(label = cutoff)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_log10(labels = scDNAmClock:::fancy_scientific) +
  ylab("# sites available") +
  xlab("missing proportion")
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/data_integration/figs/missing_cutoff_3.png", width = 4, height = 3)

# choose 0.2 
data6 <- data5[miss_sites < 0.2, ]
miss_patients <- patient_missingness(data6)

cutoffs <- seq(0.1, 0.9, by = 0.1)
sites_avail <- c()
missingness <- c()
for (cutoff in cutoffs) {
  data7 <- data6[, c(T, T, miss_patients < cutoff)]
  sites_avail <- c(sites_avail, ncol(data7) - 2)
  missingness <- c(missingness, missing_proportion(data7))
}
dat <- data.frame(cutoff = cutoffs, sites_avail = sites_avail,
                  missingness = missingness)
p <- ggplot(data = dat, aes(x = missingness, y = sites_avail)) +
  geom_text(aes(label = cutoff)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(labels = scDNAmClock:::fancy_scientific) +
  ylab("# patients available") +
  xlab("missing proportion")
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/data_integration/figs/missing_cutoff_4.png", width = 4, height = 3)

# choose 0.2
data8 <- data7[, c(T, T, miss_patients < 0.2)]

miss_sites <- site_missingness(data8)
cutoffs <- seq(0.1, 0.9, by = 0.1)
sites_avail <- c()
missingness <- c()
for (cutoff in cutoffs) {
  data9 <- data8[miss_sites < cutoff, ]
  sites_avail <- c(sites_avail, nrow(data9))
  missingness <- c(missingness, missing_proportion(data9))
}
dat <- data.frame(cutoff = cutoffs, sites_avail = sites_avail,
                  missingness = missingness)
p <- ggplot(data = dat, aes(x = missingness, y = sites_avail)) +
  geom_text(aes(label = cutoff)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_log10(labels = scDNAmClock:::fancy_scientific) +
  ylab("# sites available") +
  xlab("missing proportion")
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/data_integration/figs/missing_cutoff_5.png", width = 4, height = 3)

# 44451 x 109 patients (each 20 % x 20 %, 6 % overall missingness)
missing_proportion(data8)
imputed <- impute::impute.knn(as.matrix(data8[,-c(1,2)]), k = 10)
dim(imputed)
imputed <- cbind(data8[, 1:2], as.data.frame(imputed$data))
saveRDS(imputed, file = "data-raw/imputed_liver_data.rds")
# summary statistics imputed
hist(as.numeric(imputed[,-c(1:2)]), breaks = 30)
avail_geo <- colnames(imputed)[-c(1,2)]

avail_metadata <- metadata %>%
  filter(geo_accession %in% avail_geo) 
