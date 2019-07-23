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
