library(ggplot2)
library(tidyverse)
library(dendextend)

# reproducing the problem using technical replicates
samples <- readr::read_csv("data-raw/technical_rep_sample_lists.csv")
data <- readr::read_tsv("data-raw/tech_rep_data.tsv")
colnames(data) <- c("probe", samples$name)
mapped <- data %>%
  tidyr::drop_na() %>%
  tidyr::gather(pid, beta, -probe) %>%
  left_join(samples[, -2], by = c("pid" = "name")) %>%
  select(-pid) %>%
  tidyr::spread(group, beta)
colnames(mapped)[5:6] <- c("replicate 1", "replicate 2")
readr::write_csv(mapped, path = "data-raw/tech_rep_data_annotated.csv")


probes_data <- split(mapped, mapped$probe)

probe_ICC <- purrr::map(probes_data, 
                        ~irr::icc(.x[, c(5, 6)], 
                             model = "oneway",
                             type = "agreement",
                             unit = "single"))
res <- data.frame(
  probe = names(probes_data),
  icc = purrr::map_dbl(probe_ICC, ~.x$value),
  lwr_icc = purrr::map_dbl(probe_ICC, ~.x$lbound),
  upr_icc = purrr::map_dbl(probe_ICC, ~.x$ubound)
)

mapped <- mapped %>%
  left_join(res, by = "probe")
mapped$icc_truncated <- mapped$icc
mapped$icc_truncated[mapped$icc_truncated < 0] <- 0
# in this analysis, we want to check the absolute consistency of the each probe
# on two samples (probe x replicate x samples) 
# In a interrater reliability analysis, we have raters who evaluate subjects
# In this case, sample = subjet, replicate = rater => mean reliability of the raters (reflects probe reliability)

# scatterplot of the 513 sites
p <- ggplot(data = mapped %>% rename(ICC = icc_truncated)) +
  geom_point(aes(x = `replicate 1`, y = `replicate 2`, color = ICC), size = 0.1, alpha = 0.1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8, 0.3)) +
  scale_color_viridis_c() 
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/abs_cor.png", width = 4, height = 4)

# low consistency probes
probe_cor <- mapped %>%
  distinct(probe, cor, icc_truncated)

p <- ggplot(data = probes_data$cg13631913) +
  geom_point(aes(x = `replicate 1`, y = `replicate 2`), size = 2, alpha = 0.3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8, 0.3),
        plot.title = element_text(hjust = 0.5, size = 9)) +
  ggtitle("cg13631913 (ICC = 0.0151)")
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/low_cor_1.png", width = 2.6, height = 2.6)


p <- ggplot(data = probes_data$cg26109803) +
  geom_point(aes(x = `replicate 1`, y = `replicate 2`), size = 2, alpha = 0.3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8, 0.3),
        plot.title = element_text(hjust = 0.5, size = 9)) +
  ggtitle("cg26109803 (ICC = 0.0068)")
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/low_cor_2.png", width = 2.6, height = 2.6)

min_max <- purrr::imap(probes_data, function(x, y) {
  data.frame(probe = y, 
             low = min(c(x$`replicate 1`, x$`replicate 2`)),
             high = max(c(x$`replicate 1`, x$`replicate 2`)),
             mean = mean(c(x$`replicate 1`, x$`replicate 2`)),
             sd = sd(c(x$`replicate 1`, x$`replicate 2`)))
}) %>%
  dplyr::bind_rows() %>%
  left_join(probe_cor, by = "probe")


p <- ggplot(data = min_max %>% rename(ICC = icc_truncated)) +
  geom_hline(aes(yintercept = 0.5), size = 0.4, linetype = "dashed") +
  geom_hline(aes(yintercept = 0.75), size = 0.4, linetype = "dashed") +
  geom_hline(aes(yintercept = 0.9), size = 0.4, linetype = "dashed") +
  geom_point(aes(x = sd, y = ICC, color = mean), size = 1, alpha = 0.9) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8, 0.3),
        plot.title = element_text(hjust = 0.5, size = 9)) +
  scale_color_viridis_c() 
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/var_cor.png", width = 4, height = 3.5)

# beta values as a mixture
plot_probe <- function(data, probe_id) {
  x <- data %>% 
    dplyr::filter(probe %in% probe_id) %>%
    tidyr::gather(sid, value, -probe) %>%
    mutate(probe = factor(probe, levels = probe_id))
  ggplot(data = x) +
    geom_histogram(aes(x = value), bins = 50) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_wrap(~probe) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    xlab("beta")
}

head(mapped %>% distinct(probe, icc_truncated) %>% arrange(icc_truncated) %>% .$probe)

probe_list <- c("cg18771300", "cg19566405", "cg23159337", "cg26357744", "cg17324128", "cg08424423", "cg00230271", "cg01254459", "cg01651821")
p <- plot_probe(data, probe_list)
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/histogram.png", width = 6, height = 5)


# filter for ICC < 0.5
probe_reliability <- mapped %>%
  distinct(probe, icc_truncated) 

probe_reliability$category <- "excellent (> 0.9)"
probe_reliability$category[probe_reliability$icc_truncated < 0.9] <- "good (0.75 - 0.9)"
probe_reliability$category[probe_reliability$icc_truncated < 0.75] <- "moderate (0.5 - 0.75)"
probe_reliability$category[probe_reliability$icc_truncated < 0.5] <- "poor (< 0.5)"
kableExtra::kable(table(probe_reliability$category)) %>%
  kableExtra::kable_styling()

poor_probes <- probe_reliability %>%
  filter(category == "poor (< 0.5)")
readr::write_csv(poor_probes, path = "data-raw/poor_reliability_probes.csv")

poor_probes_data <- mapped %>%
  filter(probe %in% poor_probes$probe) %>%
  group_by(probe, sample) %>%
  summarize(mean_beta = mean(`replicate 1`, `replicate 2`))

mat <- as.matrix(tidyr::spread(poor_probes_data, probe, mean_beta)[,-1])
bicor_dist <- function(x) {
  as.dist(1 - WGCNA::bicor(x))
} 

dissim <- 1 - bicor_dist(mat)
tree <- hclust(dissim, method = "ward.D2")
plot(tree, labels = F)

cor_mat <- WGCNA::bicor(mat) 
# png(file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/cor_plot.png", width = 600, height = 600)
# res <- corrplot::corrplot(cor_mat, diag = F, 
#                    col = colorRampPalette(c("#0047BB", "white", "#D1350F"))(10), 
#                    tl.pos = "n", order = "hclust", hclust.method = "ward.D2",
#                    is.corr = F, method = "color", type = "upper")
# dev.off()



res_clust <- hclust(dist(cor_mat), method = "ward.D2")
dend <- as.dendrogram(res_clust) %>%
  set("branches_k_color", k = 4) %>% 
  ladderize
pdf(file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/cor_dend.pdf", width = 6, height = 5)
gplots::heatmap.2(cor_mat, trace = "none", labRow = F, labCol = F, col = colorRampPalette(c("#0047BB", "white", "#D1350F"))(10),
                  Rowv = dend, Colv = dend, dendrogram = "column", symm = T,
                  margins = c(1, 1))
dev.off()

# examine the characteristics of the clusters
cuts <- cutree(dend, k = 4)
table(cuts)

# within cluster characteristics
cor_mat[lower.tri(cor_mat,diag=TRUE)] <- NA # put NA
cor_mat<-as.data.frame(as.table(cor_mat)) # as a dataframe
cor_mat<-na.omit(cor_mat) # remove NA
cor_mat<-cor_mat[with(cor_mat, order(-Freq)), ]

cluster_map <- hashmap::hashmap(names(cuts), cuts)
cor_mat$Var1_cluster <- purrr::map_chr(cor_mat$Var1, ~cluster_map[[.x]])
cor_mat$Var2_cluster <- purrr::map_chr(cor_mat$Var2, ~cluster_map[[.x]])
cor_mat %>%
  group_by(Var1_cluster, Var2_cluster) %>%
  rename(cluster_mean_cor = Var1_cluster) %>%
  summarize(mean_cor = mean(Freq)) %>%
  tidyr::spread(Var2_cluster, mean_cor) %>%
  kableExtra::kable(., digits = 2) %>%
  kableExtra::kable_styling()

res <- t(mat) %>%
  as.data.frame() %>%
  mutate(probe = colnames(mat)) 
poor_probes$cluster <- purrr::map_chr(poor_probes$probe, ~cluster_map[[.x]])
poor_probes$ICC <- poor_probes$icc_truncated
p <- ggplot(data = poor_probes ) + 
  geom_histogram(aes(x = ICC)) +
  facet_wrap(~cluster) + 
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/ICC_dist.pdf", width = 4, height = 3.5)
