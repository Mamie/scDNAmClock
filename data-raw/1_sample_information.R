library(dplyr)
library(ggplot2)

# preprocess the metadata

# horvath 2018
dat <- readr::read_tsv("~/Downloads/GSE120132.tsv", col_names = F)
dat <- t(dat)
colnames(dat) <- dat[1,]
dat <- dat[-1,]
colnames(dat) <- c("name", "geo_accession", "source", "organism",
                   "strain", "tissue", "age", "sex")

dat <- as.data.frame(dat)
dat$strain <- gsub("strain: ", "", dat$strain)
dat$tissue <- gsub("tissue: ", "", dat$tissue)
dat$age <- as.numeric(gsub("age: |mo", "", dat$age))
dat$sex <- gsub("Sex: ", "", dat$sex)
readr::write_tsv(dat, path = "data-raw/Horvath2018_GSE120132.tsv")

# cole 2017

dat <- readr::read_tsv("~/Downloads/GSE89275.txt", col_names = F)
dat <- t(dat)
colnames(dat) <- dat[1,]
dat <- dat[-1,]
colnames(dat) <- c("name", "geo_accession", "source", "organism",
                   "strain", "characteristics")

dat <- as.data.frame(dat)
dat$strain <- gsub("strain: |sample type: ", "", dat$strain)
dat$characteristics <- gsub("genotype: |treatment: ", "", dat$characteristics)
dat$age <- as.numeric(purrr::map_chr(dat$name, ~stringr::str_match(.x, "([0-9]+) [Mm]onth")[2]))
dat$sex <- "Female"
readr::write_tsv(dat, path = "data-raw/Cole2017_GSE89275.tsv")

# Hahn 2017

dat <- readr::read_tsv("~/Downloads/GSE92486.txt", col_names = F)
dat <- t(dat)
colnames(dat) <- dat[1,]
dat <- dat[-1,]
colnames(dat) <- c("name", "geo_accession", "organism",
                   "strain", "sex", "characteristics", "age", "source")
dat <- as.data.frame(dat)
dat$strain <- gsub("strain: ", "", dat$strain)
dat$characteristics <- gsub("diet: ", "", dat$characteristics)
dat$age <- as.numeric(purrr::map_chr(dat$age, ~stringr::str_match(.x, "([0-9]+) [Mm]onth")[2]))
dat$sex <- "Female"
dat$source <- purrr::map_chr(dat$source, ~gsub("tissue: ", "", .x))
dat$data <- purrr::map_chr(dat$name, ~tail(strsplit(as.character(.x), split = " ", fixed = T)[[1]], 1))
readr::write_tsv(dat[, c("name", "geo_accession", "source", "organism", "strain", "age", "sex", "characteristics", "data")], path = "data-raw/Hahn2017_GSE92486.tsv")

# Zhou 2016

dat <- readr::read_tsv("~/Downloads/GSE85772.txt", col_names = F)
dat <- t(dat)
colnames(dat) <- dat[1,]
dat <- dat[-1,]
colnames(dat) <- c("name", "geo_accession", "source",
                   "organism", "strain", "age", "characteristics")
dat <- as.data.frame(dat)
dat$strain <- gsub("species: ", "", dat$strain)
dat$age <- as.numeric(purrr::map_chr(dat$age, ~stringr::str_match(.x, "([0-9]+) [Mm]onth")[2]))
dat$sex <- "Male"
dat$characteristics <- 
  purrr::map_chr(as.character(dat$name), ~strsplit(.x, split = "_")[[1]][2])
dat$characteristics[27] <- "beta-cat_TG"
dat$characteristics[28] <- "beta-cat_TG+DEN"

readr::write_tsv(dat[, c("name", "geo_accession", "source", "organism", "strain", "age", "sex", "characteristics")], path = "data-raw/Zhou2016_GSE85772.tsv")

# Stubbs 2017

dat <- readr::read_tsv("~/Downloads/GSE93957.txt", col_names = F)
dat <- t(dat)
colnames(dat) <- dat[1,]
dat <- dat[-1,]
colnames(dat) <- c("name", "geo_accession", 
                   "organism", "age", 
                   "source", "strain")
dat <- as.data.frame(dat)
dat$age <- as.numeric(purrr::map_chr(dat$age, ~stringr::str_match(.x, "([0-9]+)")[2])) * 7 / 30
dat$sex <- "Male"
dat$source <- gsub("tissue: ", "", dat$source)
dat$strain <- gsub("strain: ", "", dat$strain)
readr::write_tsv(dat[, c("name", "geo_accession", "source", "organism", "strain", "age", "sex")], path = "data-raw/Stubbs2017_GSE93957.tsv")


# combine the datasets together
all <- purrr::map(c("data-raw/Horvath2018_GSE120132.tsv", 
  "data-raw/Cole2017_GSE89275.tsv",
  "data-raw/Hahn2017_GSE92486.tsv",
  "data-raw/Zhou2016_GSE85772.tsv",
  "data-raw/Stubbs2017_GSE93957.tsv"), 
    ~readr::read_tsv(.x) %>% 
    mutate(dataset = strsplit(as.character(.x), split = "[/_]", fixed = F)[[1]][2])) %>%
  dplyr::bind_rows() %>%
  filter(source %in% c("liver", "Liver", "Hepatocytes")) %>%
  select(-tissue) %>%
  mutate(source = "Liver")
all$characteristics[is.na(all$characteristics)] <- "Untreated"
all$characteristics[grepl("MO", all$characteristics)] <- "Untreated"
all$data[is.na(all$data)] <- "BS-seq"
all <- all %>%
  filter(data == "BS-seq") %>%
  select(-data)
all$characteristics[all$characteristics %in% c("control-fed (AL)", "NCC")] <- "Untreated"
all$characteristics[all$characteristics %in% c("dietary-restricted (DR)", "Calorie Restriction Treated")] <- "Calorie restriction (CR)"
all$characteristics[all$characteristics %in% c("FFC")] <- "Fast food"
all$characteristics[all$characteristics %in% c("FFE")] <- "Fast food + exercise"
all$characteristics[all$characteristics %in% c("NCE")] <- "Exercise"
all$characteristics[all$characteristics %in% c("DEN")] <- "Liver damage"
all$characteristics <- factor(all$characteristics, levels = c("Prop1 Het", "Prop1 Hom", "Calorie restriction (CR)", "Rapamycin Treated", "Exercise", "Fast food", "Fast food + exercise", "Liver damage", "Untreated"))
readr::write_tsv(all, path = "data-raw/liver_samples.tsv")
p <- all %>%
  ggplot(data = .) +
  geom_histogram(aes(x = age, fill = characteristics), bins = 30) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank()) +
  scale_fill_brewer(palette = "Set3") +
  xlab("age (month)")
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/data_integration/figs/sample_information.png", width = 5, height = 2.5)

kableExtra::kable(table(all$characteristics)) %>%
  kableExtra::kable_styling()
