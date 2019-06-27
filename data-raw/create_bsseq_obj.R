#' ---
#' title: "Mouse ESC"
#' author: "Mamie Wang"
#' date: "June 20 2019"
#' ---

#' We are testing the utility functions to preprocess the single cell bisulfite
#' sequencing data. We used scM&T-seq dataset from [Angermueller 2016 Nature Methods](https://www.nature.com/articles/nmeth.3728) as an example.
#' The dataset (GEO accession: [GSE68642](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68642)) provides 81 mouse embryonic stem cells that passed the quality filter, 
#' 16 of which was cultured in '2i' media and 65 in serum. The paper mentioned 
#' that in serum, ESCs will be metastable and heterogenous, switching between  
#' transcriptional states, whereas 2i ESC will have hypomethylation. 

library(scDNAmClock)

in_dir <- "/gpfs/ysm/project/mw957/data/public/scMT-seq"
cov_files <- list.files(in_dir, pattern = "cov.gz", full.names = T)
id <- purrr::map_chr(cov_files, ~strsplit(.x, split = "[/.]")[[1]][9])

data <- read_meth(
  files = cov_files,
  chr_idx = 1,
  pos_idx = 2,
  met_idx = 5, 
  unmet_idx = 6,
  strand_idx = NULL,
  id = id,
  deduplicate = T)
saveRDS(data, file = "/gpfs/ysm/project/mw957/data/processed/scMT-seq/deduplicated_meth_call.rds")

data <- readRDS("/gpfs/ysm/project/mw957/data/processed/Angermueller2016/deduplicated_meth_call.rds")

# summary statistics on the data

# plot sites covered
sites_covered <- purrr::map_int(data@listData, nrow)
p <- ggplot(data = data.frame(n_sites = sites_covered, 
                              type = purrr::map_chr(names(data), ~strsplit(.x, split = "_")[[1]][1]))) +
  geom_histogram(aes(x = n_sites, fill = type), bins = 30) +
  theme_classic() +
  xlab("Number of sites covered") +
  scale_x_continuous(labels = fancy_scientific) +
  theme(axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12), 
        legend.position = c(0.8, 0.8), 
        legend.title = element_blank())

ggsave(p, 
       filename = "~/Dropbox/600 Presentations/Yale projects/scMethylseq/figs/20190624_sites_covered.png",
       width = 5, height = 3)

summary(sites_covered)

metadata <- data.frame(id = names(data),
           n_sites = sites_covered) 
readr::write_csv(metadata, "~/Dropbox/600 Presentations/Yale projects/scMethylseq/figs/data_information.csv")

# coverage of each site
p_list <- plot_coverage(data@listData[c("Serum_G02", "2I_A11", "Serum_A07", "2I_E12", "2I_E11", "Serum_D02")])
p <- cowplot::plot_grid(plotlist = p_list, nrow = 2, align = "v")
ggsave(p, filename = "~/Dropbox/600 Presentations/Yale projects/scMethylseq/figs/coverage.png", width = 9, height = 5)

joined <- join_meth_list(data, all.x = FALSE, all.y = FALSE) 
saveRDS(joined, "/gpfs/ysm/project/mw957/data/processed/Angermueller2016/inner_joined_meth.rds")

joined <- readRDS("data-raw/inner_joined_meth.rds")
dim(joined)

# list of tasks to evaluate the sites
# 81 cells
# let's select for most covered 10, 20, ..., 80 cells
metadata %>%
  arrange(n_sites) %>%
  top_n(10, n_sites)

# The script are located in https://github.com/Mamie/scDNAmClock_scripts/blob/master/20190624/common_sites.sh
res <- data.frame(
  top_n = c(2, 4, 6, 8, seq(10, 80, by = 10)),
  overlap = c(2004397, 247978, 47547, 15466, 7996, 2922, 2257, 1705, 1473, 1325, 1199, 1199)
)
readr::write_csv(res, path = "~/Dropbox/600 Presentations/Yale projects/scMethylseq/figs/overlapping_sites.csv")

p <- ggplot(data = res, aes(x = top_n, y = overlap)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 2) +
  theme_classic() +
  scale_y_log10(labels = fancy_scientific) +
  xlab("# cells with highest coverage") +
  ylab("Overlapping sites") +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) 
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/scMethylseq/figs/overlap_n.png", width = 5, height = 3)

# find the regions for probe
# EDARADD (cg09809672)
# ELOVL2 (cg16867657)
# cluster CG locus database
# https://support.illumina.com/content/dam/illumina-marketing/documents/products/technotes/technote_cpg_loci_identification.pdf
clock_info <- readr::read_csv("/gpfs/ysm/project/mw957/data/public/450K_manifest/HumanMethylation450_15017482_v1-2.csv",
                skip = 7) %>%
  filter(IlmnID %in% c("cg09809672", 
                       "cg16867657", 
                       "cg02228185",
                       "cg25809905",
                       "cg17861230"))
readr::write_csv(clock_info, "/gpfs/ysm/project/mw957/data/public/450K_manifest/probes_of_interest.csv")

clock_info %>% select(IlmnID, CHR, MAPINFO, UCSC_RefGene_Name) %>%
  mutate(UCSC_RefGene_Name = purrr::map_chr(UCSC_RefGene_Name, ~strsplit(.x, split = ";")[[1]][1])) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

# the genomic coordinate is human genome build37/hg19, 
# the genome build in the scM&T-seq paper is GRCm38/mm10,
# lifting the human genome to mouse genome 
# http://genome.ucsc.edu/cgi-bin/hgLiftOver
clock_info <- readr::read_csv("/gpfs/ysm/project/mw957/data/public/450K_manifest/probes_of_interest.csv")
regions <- paste0("chr", clock_info$CHR, ":", clock_info$MAPINFO - 500,
                 "-", clock_info$MAPINFO + 500)

info <- data.frame(
  gene = purrr::map_chr(clock_info$UCSC_RefGene_Name, ~strsplit(.x, split = ";")[[1]][1]),
  CGI = clock_info$UCSC_CpG_Islands_Name,
  region = regions,
  stringsAsFactors = F)

info$CGI_mm10 <- c("chr13:12519445-12520370",
                   "chr13:41220215-41220833",
                   NA,
                   NA,
                   "chr8:70730146-70730302")
info$region_mm10 <- c("chr13:12520677-12522012",
                      "chr13:41220278-41220833",
                      "chr11:73323854-73324841",
                      "chr11:102470087-102471373",
                      "chr8:70730146-70730302")
task_list <- info %>%
  select(gene, CGI_mm10, region_mm10) %>%
  tidyr::gather(type, region, -gene) %>%
  tidyr::drop_na() %>%
  tidyr::unite(id, gene, type)

# now subset for cells within these regions, and see the numbers of hits within
# each cell
task_list$chr <- purrr::map_chr(task_list$region, ~strsplit(.x, split = "chr|:")[[1]][2])
task_list$start <- purrr::map_chr(task_list$region, ~strsplit(.x, split = ":|-")[[1]][2])
task_list$end <- purrr::map_chr(task_list$region, ~strsplit(.x, split = ":|-")[[1]][3])
readr::write_csv(task_list, 
                 path = "/gpfs/ysm/project/mw957/data/public/450K_manifest/probe_of_interest_mm10.csv")

library(Sushi)

info <- task_list %>%
  dplyr::select(chr, start, end, id) %>%
  mutate(avail = c(63, 62, 37, 61),
         type = c(1, 1, 2, 2))

chrom = "13"
chromstart = 12519445
chromend = 12522012
plotBed(beddata = info[c(1,3),], chrom = chrom, chromstart = chromstart,
        chromend = chromend, colorby = info$type[c(1,3)], 
        colorbycol = viridis::viridis, splitstrand=TRUE)
labelgenome(chrom,chromstart,chromend,n=2,scale="Kb")
legend("topright", legend=c("CGI","500 bp"),fill=viridis::viridis(2),
       border=viridis::viridis(2), text.font=2,cex=0.6)

to_df <- function(methCall_obj, name) {
  if(n_row(methCall_obj) == 0) return()
  df <- data.frame(
    chr = as.character(methCall_obj@data$chr),
    start = methCall_obj@data$position,
    end = methCall_obj@data$position,
    methylated = methCall_obj@data$unmet_reads < methCall_obj@data$met_reads,
    name = name,
    stringsAsFactors = F
  )
  df
}

task_list$avail <- c(63, 62, 42, 37, 61, 24, 45, 42)
task_list %>%
  select(-region) %>%
  arrange(desc(avail)) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

for (i in seq_along(task_list$id)) {
  object <- readRDS(paste0("data-raw/", task_list$id[i], ".rds"))
  object <- purrr::imap(object, to_df)
  object <- dplyr::bind_rows(object)
  p <- ggplot(data = object) +
    geom_histogram(aes(x = start, fill = methylated), bins = 100) +
    theme_classic() +
    theme(legend.position = c(0.8, 0.8)) +
    xlab(paste0("chr", task_list$chr[i]))
  ggsave(p, filename = paste0("~/Dropbox/600 Presentations/Yale projects/scMethylseq/figs/", task_list$id[i], ".png"), width = 5, height = 3)
}

chrom = "13"
chromstart = 41220215
chromend = 41220833
plotBed(beddata = info[c(2,4),], chrom = chrom, chromstart = chromstart,
        chromend = chromend)
labelgenome(chrom,chromstart,chromend,n=2,scale="Kb")

# find the number of cells having overlapping marker


# are cells with no-methylation at certain sites 
data <- data.frame()
for (i in seq(4, 8)) {
  object <- readRDS(paste0("data-raw/", task_list$id[i], ".rds"))
  object <- purrr::imap(object, to_df)
  object <- dplyr::bind_rows(object)
  object$gene <- strsplit(task_list$id[i], split = "_")[[1]][1]
  data <- dplyr::bind_rows(data, object)
}

# Let's bin each region into 100 bins as the result of ggplot
data %<>% group_by(gene, name) %>%
  summarize(mean_methyl = mean(methylated)) 
data$condition <- purrr::map_chr(data$name, ~strsplit(.x, split = "_")[[1]][1])
data_wide <- data %>%
  tidyr::spread(gene, mean_methyl)

corrplot::corrplot(cor(as.matrix(data_wide[,-1]), use = "pairwise.complete.obs"),
                   type = "lower", diag = F, order = "hclust")

p <- ggplot(data = data) +
  geom_histogram(aes(x = mean_methyl, fill = condition), bins = 10) +
  facet_wrap(~gene, scale = "free", ncol = 5) +
  theme_classic() +
  theme(legend.position = c(0.95, 0.8),
        strip.background = element_blank()) +
  xlab("mean methylation")
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/scMethylseq/figs/mean_methylation.png", width = 10, height = 2.5)

p <- GGally::ggpairs(data_wide, mapping = aes(color = condition, alpha = 0.6), columns = 2:7) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
ggsave(p, file = "~/Dropbox/600 Presentations/Yale projects/scMethylseq/figs/correlation.png", width = 12, height = 6)

# add the sample size for each 
data_wide %>%
  select(PDE4C, ITGA2B) %>%
  tidyr::drop_na()

test <- t(read.delim("~/Downloads/raw.txt", header = F))[-1, c(1,2,4,5,3,6)]
colnames(test) <- c("id", "name", "age", "sex", "cell_type", "source")
test <- as.data.frame(test)
test[, "age"] <- purrr::map_int(test[, "age"], ~as.integer(gsub("age: ", "", .x)))
test[, "sex"] <- purrr::map_chr(test[, "sex"], ~gsub("gender: ", "", .x))
test[, "cell_type"] <- purrr::map_chr(test[, "cell_type"], ~gsub("cell type: ", "", .x))
readr::write_csv(test, path = '~/Downloads/sampl2_info.csv')
