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

# sites covered
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
  filter(IlmnID %in% c("cg09809672", "cg16867657"))
readr::write_csv(clock_info, "/gpfs/ysm/project/mw957/data/public/450K_manifest/EDARADD_ELOVL2.csv")

clock_info <- readr::read_csv("data-raw/EDARADD_ELOVL2.csv")

# the genomic coordinate is human genome build37/hg19, 
# the genome build in the scM&T-seq paper is GRCm38/mm10,
# lifting the human genome to mouse genome 
# http://genome.ucsc.edu/cgi-bin/hgLiftOver
clock_info$UCSC_CpG_Islands_Name
#  "chr1:236558459-236559336" "chr6:11043913-11045206"
#  "chr13:12519445-12520370" "chr13:41220215-41220833"
paste0("chr", clock_info$CHR, ":", clock_info$MAPINFO)
# expand for 500 bp up and downstream
# "chr1:236557682" "chr6:11044877" => "chr1:236557182-236558182" "chr6:11044377-11045377" 
# "chr13:12520677-12522012" "chr13:41220278-41220833"

# now subset for cells within these regions, and see the numbers of hits within
# each cell
task_list <- data.frame(
  region = c("chr13:12520677-12522012",
             "chr13:41220278-41220833"),
  id = c("EDARADD", 
         "ELOVL2"),
  stringsAsFactors = F)
task_list$chr <- strsplit(task_list$region, split = "chr|:")[[1]][2]
task_list$start <- purrr::map_chr(task_list$region, ~strsplit(.x, split = ":|-")[[1]][2])
task_list$end <- purrr::map_chr(task_list$region, ~strsplit(.x, split = ":|-")[[1]][3])

readr::write_csv(task_list, 
                 path = "/gpfs/ysm/project/mw957/data/public/450K_manifest/EDARADD_ELOVL2_mm10.csv")


class(as.data.frame(data@listData[1]))
