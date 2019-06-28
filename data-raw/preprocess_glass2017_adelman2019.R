library(ggplot2)
library(dplyr)
library(scDNAmClock)

# preprocess Glass2017 and Adelman2019 ERRBS dataset
sample_info <- readr::read_csv("data-raw/glass2017_adelman2019_sample_info.csv")
table(sample_info$sex, sample_info$cell_type)


p <- sample_info %>%
  mutate(age = as.numeric(age)) %>%
ggplot(data = .) +
  geom_histogram(aes(x = age, fill = sex)) +
  facet_wrap(~cell_type, ncol = 1) +
  theme_bw() +
  theme(panel.background = element_blank(), 
        strip.background = element_blank(),
        panel.grid = element_blank())
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/scRNAseq_RRBS/figs/eRRBS_sample.png",
       width = 5, height = 4)

########## Adelman 2019 #################################################

adelman2019 <- list.files("/gpfs/ysm/project/mw957/data/public/Adelman2019",
                          pattern = "myCpG.txt.gz")

# format 
# chrBase	chr	base	strand	coverage	freqC	freqT
# chr1.10631	chr1	10631	F	112	95.536	4.464
# chr1.10542	chr1	10542	F	80	100.000	0.000
#
# C represents methylated cytosine
# T represents unmethylated cytosine

id <- purrr::map_chr(adelman2019, ~strsplit(.x, split = "_|/")[[1]][3])

data <- read_meth(
  files = cov_files,
  chr_idx = 2,
  pos_idx = 3,
  met_idx = 6, 
  unmet_idx = 7,
  strand_idx = 4,
  coverage_idx = 5,
  id = id,
  deduplicate = F,
  header = T
)
saveRDS(data, file = "/gpfs/ysm/project/mw957/data/processed/glass2017_adelman2019/adelman2019_meth_call.rds")

adelman2019_sample_info <- sample_info %>%
  filter(cell_type == "HSCe") 
adelman2019_sample_info$category <- purrr::map_chr(adelman2019_sample_info$name,
                                                   ~strsplit(.x, split = "_")[[1]][1])

task_list <- data.frame(gene = c("ELOVL2", "PDE4C"),
                        position = c(11044877, 18343901))

for (i in seq(nrow(task_list))) {
  data <- readRDS(paste0("data-raw/adelman2019/adelman2019_", task_list$gene[i], ".rds"))
  data <- purrr::imap(data, ~mutate(as_df(.x), name = .y)) %>%
    dplyr::bind_rows() 
  p <- data %>%
    left_join(adelman2019_sample_info, by = c("name" = "id")) %>%
    ggplot(data = .) +
    geom_jitter(aes(x = position, y = met_reads, color = category), size = 1) +
    theme_bw() + 
    theme(panel.background = element_blank(), 
          strip.background = element_blank(),
          panel.grid = element_blank()) +
    scale_color_viridis_d() +
    geom_vline(aes(xintercept = task_list$position[i])) +
    geom_vline(aes(xintercept = task_list$position[i] - 5), size = 0.1) +
    geom_vline(aes(xintercept = task_list$position[i] + 5), size = 0.1) +
    ylab("% methylated")
  ggsave(p, filename = paste0("~/Dropbox/600 Presentations/Yale projects/scRNAseq_RRBS/figs/adelman2019_", task_list$gene[i], ".png"), width = 6, height = 3)
}

# average using a bin of +/- 5 bp

for (i in seq(nrow(task_list))) {
  data <- readRDS(paste0("data-raw/adelman2019/adelman2019_", task_list$gene[i], ".rds"))
  p <- purrr::imap(data, ~mutate(as_df(.x), name = .y)) %>%
    dplyr::bind_rows() %>%
    filter(position > task_list$position[i] - 5 & position < task_list$position[i] + 5) %>%
    group_by(name) %>%
    summarize(mean_meth = mean(met_reads)) %>%
    left_join(adelman2019_sample_info, by = c("name" = "id")) %>%
    mutate(age = as.numeric(age)) %>%
    ggplot(data = .) +
    geom_point(aes(x = age, y = mean_meth)) +
    theme_bw() + 
    theme(panel.background = element_blank(), 
          strip.background = element_blank(),
          panel.grid = element_blank()) +
    ylab("mean % methylation")
  ggsave(p, filename = paste0("~/Dropbox/600 Presentations/Yale projects/scRNAseq_RRBS/figs/adelman2019_", task_list$gene[i], "_line.png"), width = 3, height = 3)
}


########## Glass 2017 #################################################

glass2017_sample_info <- sample_info %>%
  filter(cell_type == "blast") %>%
  mutate(age = as.numeric(age))

glass2017 <- list.files("/gpfs/ysm/project/mw957/data/public/Glass2017", 
                        pattern = "mincov10.txt.gz")

# format
# chrBase	chr	base	strand	coverage	freqC	freqT
# chr1.10497	chr1	10497	F	501	97.01	2.99
# chr1.10563	chr1	10563	F	11	100.00	0.00

id <- purrr::map_chr(glass2017, ~strsplit(.x, split = "_|/")[[1]][3])

data <- read_meth(
  files = glass2017,
  chr_idx = 2,
  pos_idx = 3,
  met_idx = 6, 
  unmet_idx = 7,
  strand_idx = 4,
  coverage_idx = 5,
  id = id,
  deduplicate = F,
  header = T
)
saveRDS(data, file = "/gpfs/ysm/project/mw957/data/processed/glass2017_adelman2019/glass2017_meth_call.rds")

kableExtra::kable(sample_info %>%
                    filter(cell_type == "HSCe")) %>%
  kableExtra::kable_styling()

glass2017 <- readRDS("data-raw/adelman2019.rds")
purrr::map_int(glass2017, ~min(.x@data$coverage))
# whole genome coverage

clock_info <- readr::read_csv("/gpfs/ysm/project/mw957/data/public/450K_manifest/probes_of_interest.csv") # genome build 37/hg19
regions <- paste0("chr", clock_info$CHR, ":", clock_info$MAPINFO - 500,
                  "-", clock_info$MAPINFO + 500)
# ERRBS is also using hg 19 genome reference
info <- data.frame(
  probe = clock_info$IlmnID,
  chr = clock_info$CHR,
  position = clock_info$MAPINFO,
  gene = purrr::map_chr(clock_info$UCSC_RefGene_Name, ~strsplit(.x, split = ";")[[1]][1]),
  region = regions,
  clock = c("Blocklandt 2011", "Garangnani 2012", rep("Weidner 2014", 3)),
  stringsAsFactors = F)
info$avail <- c(1, 12, 2, 0, 12)
info$avail <- c(4, 119, 11, 1, 119)
kableExtra::kable(info) %>%
  kableExtra::kable_styling()
task_list <- info %>%
  select(gene, region)

task_list$chr <- purrr::map_chr(task_list$region, ~strsplit(.x, split = "chr|:")[[1]][2])
task_list$start <- purrr::map_chr(task_list$region, ~strsplit(.x, split = ":|-")[[1]][2])
task_list$end <- purrr::map_chr(task_list$region, ~strsplit(.x, split = ":|-")[[1]][3])
readr::write_csv(task_list, 
                 path = "/gpfs/ysm/project/mw957/data/public/450K_manifest/probe_of_interest_hg19.csv")

task_list <- readr::read_csv("/gpfs/ysm/project/mw957/data/public/450K_manifest/probe_of_interest_hg19.csv")

task_list %>% 
  mutate(dataset = "adelman2019") %>%
  mutate(dir = "/gpfs/ysm/project/mw957/data/processed/glass2017_adelman2019/adelman2019/") %>%
  rbind(
    task_list %>% 
      mutate(dataset = "glass2017") %>%
      mutate(dir = "/gpfs/ysm/project/mw957/data/processed/glass2017_adelman2019/glass2017/")
  ) %>%
  readr::write_csv(., path = "/gpfs/ysm/project/mw957/data/processed/glass2017_adelman2019/region_task_list.csv")

for (i in seq(nrow(task_list))) {
  data <- readRDS(paste0("data-raw/glass2017/glass2017_", task_list$gene[i], ".rds"))
  data <- purrr::imap(data, ~mutate(as_df(.x), name = .y)) %>%
    dplyr::bind_rows() 
  p <- data %>%
    left_join(glass2017_sample_info, by = c("name" = "id")) %>%
    ggplot(data = .) +
    geom_vline(aes(xintercept = task_list$position[i]), size = 0.1) +
    geom_vline(aes(xintercept = task_list$position[i]- 5) , size = 0.1) +
    geom_vline(aes(xintercept = task_list$position[i]+ 5) , size = 0.1) +
    geom_jitter(aes(x = position, y = met_reads, color = age), size = 0.2) +
    theme_bw() + 
    theme(panel.background = element_blank(), 
          strip.background = element_blank(),
          panel.grid = element_blank()) +
    scale_color_viridis_c() +
    ylab("% methylated")
  ggsave(p, filename = paste0("~/Dropbox/600 Presentations/Yale projects/scRNAseq_RRBS/figs/glass2017_", task_list$gene[i], ".png"), width = 6, height = 3)
}

# average using a bin of +/- 5 bp

for (i in seq(nrow(task_list))) {
  data <- readRDS(paste0("data-raw/glass2017/glass2017_", task_list$gene[i], ".rds"))
  p <- purrr::imap(data, ~mutate(as_df(.x), name = .y)) %>%
    dplyr::bind_rows() %>%
    filter(position > task_list$position[i] - 5 & position < task_list$position[i] + 5) %>%
    group_by(name) %>%
    summarize(mean_meth = mean(met_reads)) %>%
    left_join(glass2017_sample_info, by = c("name" = "id")) %>%
    mutate(age = as.numeric(age)) %>%
    ggplot(data = .) +
    geom_point(aes(x = age, y = mean_meth)) +
    theme_bw() + 
    theme(panel.background = element_blank(), 
          strip.background = element_blank(),
          panel.grid = element_blank()) +
    ylab("mean % methylation")
  ggsave(p, filename = paste0("~/Dropbox/600 Presentations/Yale projects/scRNAseq_RRBS/figs/glass2017_", task_list$gene[i], "_line.png"), width = 3, height = 3)
}
