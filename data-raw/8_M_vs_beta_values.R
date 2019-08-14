# investigate the effect of M values vs Beta values
# M values = log2(beta/1-beta)
# beta = 2^M/(2^M + 1)

rep_data <- readRDS("data-raw/Lehne2015_M_Levine.rds")

dim(rep_data) # [1]  513 2712
colnames(rep_data)

samples <- readr::read_csv("data-raw/technical_rep_sample_lists.csv")
summ_dat <- rep_data %>%
  tidyr::gather(name, M, -ID_REF) %>%
  inner_join(samples[, c(1, 5, 6)], by = "name") %>%
  group_by(sample, ID_REF) %>%
  summarize(mean = mean(M), sd = sd(M))
head(summ_dat)
unique(summ_dat$sample)

ggplot(data = summ_dat) +
  geom_point(aes(x = mean, y = sd), size = 0.1, alpha = 0.2) +
  theme_classic() +
  theme(panel.grid = element_blank())
