library(ggplot2)

# reproducing the problem using technical replicates
samples <- readr::read_csv("data-raw/technical_rep_sample_lists.csv")
data <- readr::read_tsv("data-raw/tech_rep_data.tsv")
colnames(data) <- c("probe", samples$name)
mapped <- data %>%
  tidyr::gather(pid, beta, -probe) %>%
  left_join(samples[, -2], by = c("pid" = "name")) %>%
  select(-pid) %>%
  tidyr::spread(group, beta)
colnames(mapped)[5:6] <- c("replicate 1", "replicate 2")
probes_data <- split(mapped, mapped$probe)
probe_cor <- purrr::map(probes_data, ~cor.test(.x$`replicate 1`, .x$`replicate 2`, method = "pearson"))
res <- data.frame(
  probe = names(probes_data),
  cor = purrr::map_dbl(probe_cor, "estimate"),
  lwr = purrr::map_dbl(probe_cor, ~.x$conf.int[1]),
  upr = purrr::map_dbl(probe_cor, ~.x$conf.int[2])
)

mapped <- mapped %>%
  left_join(res, by = "probe")

# in this analysis, we want to check the absolute consistency of the each probe
# on two samples (probe x replicate x samples) 
# In a interrater reliability analysis, we have raters who evaluate subjects
# In this case, sample = subjet, replicate = rater => mean reliability of the raters (reflects probe reliability)

# scatterplot of the 513 sites
p <- ggplot(data = mapped) +
  geom_point(aes(x = `replicate 1`, y = `replicate 2`, color = abs(cor)), size = 0.1, alpha = 0.1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8, 0.3)) +
  scale_color_viridis_c()
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/abs_cor.png", width = 4, height = 4)
# compute ICC for the 513 sites
# intraclass correlation coefficients => inter rater agreement
# relationship with expression level
# relationship with variation in expression level
# interpretation of ICC
# Koo and Li (2016)[17]:
# below 0.50: poor
# between 0.50 and 0.75: moderate
# between 0.75 and 0.90: good
# above 0.90: excellent