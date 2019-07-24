library(dplyr)
library(tidyr)
library(ggplot2)

# hyperparameter tuning
hyperparameters <- readr::read_csv("data-raw/hyperparameters.csv")
plot_data <- tidyr::gather(hyperparameters, parameters, values, -k) %>%
  mutate(parameters = case_when(
    parameters == "mae" ~ "MAE",
    parameters == "rmse" ~ "RMSE",
    parameters == "rsq" ~ "R squared",
    parameters == "non_zero" ~ "# nonzero coefficients"
  ))
p <- ggplot(data = plot_data) +
  geom_line(aes(x = k, y = values)) +
  facet_wrap(~parameters, ncol = 1, scale = "free_y") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_blank()) +
  geom_vline(aes(xintercept = 650), linetype = "dashed", size = 0.3) 
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/hyperparameters.png", width = 2, height = 5)

# changes of the coefficients chosen
k <- sapply(files, function(x) as.numeric(gsub(".csv", "", x)))
data <- purrr::map(hyperparameters$k[-19], ~readr::read_csv(file.path("data-raw/hyperparameters", paste0(.x, ".csv"))) %>%
                     mutate(k = .x)) %>%
  dplyr::bind_rows() %>%
  filter(term != "intercept") %>%
  tidyr::spread(k, estimate) 
data$n <- apply(data[-1], 1, function(x) sum(!is.na(x)))
data <- data %>%
  arrange(desc(n))
plot_mat <- as.matrix(data[,-c(1, 20)])
colnames(plot_mat) <- colnames(data)[-c(1,20)]
rownames(plot_mat) <- data$term

row_side_colors <- data$n
colours <- viridis::inferno(18)[row_side_colors]

pdf(file = "~/Dropbox/600 Presentations/Yale projects/low_intensity_probe_correction/figs/coeffs.pdf", width = 5, height = 6.5)
gplots::heatmap.2(plot_mat, trace = "none", 
                  labRow = F, labCol = colnames(plot_mat), 
                  col = color.palette(c("#0047BB", "white", "#D1350F"), c(2, 2), space = "rgb")(4),
                  Rowv = F,
                  Colv = F,
                  dendrogram = "none",
                  RowSideColors = colours)
dev.off()

library(magrittr)
data %>%
  filter(!is.na(`650`)) %>%
  .$term %>% 
  length(.) %T>%
  intersect(., scDNAmClock:::pheno_age_dat$CpG) %>%
  length
