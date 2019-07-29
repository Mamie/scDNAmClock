# modeling of the mixture model
library(rstan)
library(dplyr)
library(bayesplot)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = T)


mapped <- readr::read_csv("data-raw/tech_rep_data_annotated.csv")
set.seed(100)
probes <- unique(mapped$probe)  #%>% sample(., 20)

all_probes <- unique(mapped$probe)
dat <- mapped %>%
  filter(probe %in% probes) %>%
  select(probe, sample, `replicate 1`, `replicate 2`) %>%
  arrange(probe, sample)
all_dat <- mapped %>%
  filter(probe %in% all_probes) %>%
  select(probe, sample, `replicate 1`, `replicate 2`) %>%
  arrange(probe, sample)
sapply(split(all_dat$`replicate 1`, all_dat$probe), sd) %>%
  quantile(.) %>%
  kableExtra::kable(., digits = 4) %>%
  kableExtra::kable_styling()
apply(all_dat[,3:4], 1, sd) %>%
  quantile(.) %>%
  kableExtra::kable(., digits = 4) %>%
  kableExtra::kable_styling()
dat_split <- purrr::map(split(dat, dat$probe), 
                                 function(x) {
                                   samples <- x$sample
                                   x <- as.matrix(x[, -c(1,2)])
                                   rownames(x) <- samples
                                   x
                                   } )

test_dat <- list(
  I = length(dat_split),
  J = 36,
  K = 2,
  betahat = dat_split
)
probes <- names(dat_split)

fit <- rstan::stan(file = "data-raw/noise_model.stan", data = test_dat, chains = 4)

saveRDS(fit, file = "data-raw/3_modeling.rds")
fit <- readRDS("data-raw/3_modeling.rds")
fit <- readRDS("data-raw/20190709_noise_model.rds")

shinystan::launch_shinystan(fit)
summary_fit <- summary(fit, pars = c("beta"))
beta_true <- summary_fit$summary

beta_true_mat <- matrix(beta_true[,1], nrow = 36, byrow = F)
colnames(beta_true_mat) <- all_probes

poor_reliability_probes <- readr::read_csv("data-raw/poor_reliability_probes.csv")
cor_mat <- WGCNA::bicor(beta_true_mat[,poor_reliability_probes$probe])

corrplot::corrplot(cor_mat, diag = F, 
                   col = colorRampPalette(c("#0047BB", "white", "#D1350F"))(10), 
                   tl.pos = "n", order = "hclust", hclust.method = "ward.D2",
                   is.corr = T, method = "color", type = "upper")

mean_replicate <- purrr::imap(dat_split, 
                             function(x, y) data.frame(beta = apply(x, 1, mean)) %>%
                               mutate(probe = y)) %>%
  dplyr::bind_rows() %>% 
  mutate(mean_betahat = beta_true[,1])

labels <- unique(mean_replicate$probe)
mean_replicate$probe <- as.character(as.numeric(factor(mean_replicate$probe, levels = labels)))
summary_replicate <- mean_replicate %>%
  group_by(probe) %>%
  summarize_all(mean) 
p <- ggplot(data = mean_replicate) + 
  geom_point(aes(x = beta, y = mean_betahat, color = probe), size = 0.1, alpha = 0.3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position  = "none",
        plot.title = element_text(hjust = 0.5, size = 9)) +
  ylab(expression(bar(beta)[replicate])) +
  xlab(expression(beta)) +
  geom_text(data = summary_replicate, aes(x = beta, y = mean_betahat, label = probe), size = 1) 
  
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/beta_hat_beta_2.png", width = 3, height = 3)

# plot_title <- ggtitle("Posterior distributions",
#                       "with medians and 95% intervals")
# p <- bayesplot::mcmc_areas(as.matrix(fit),
#            pars = c("eta"),
#            prob = 0.95) + plot_title
# ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/eta.png", width = 3, height = 2)

