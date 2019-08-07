load("~/Downloads/Rat_CpGs.RData") # map rn6 to mm10
library(rtracklayer)
library(GenomicRanges)
library(magrittr)
library(R.utils)
library(dplyr)

# download the chain file
chain_file <- "http://hgdownload.cse.ucsc.edu/goldenPath/rn6/liftOver/rn6ToMm10.over.chain.gz"
dest_file <- "~/Downloads/rn6ToMm10.over.chain.gz"
download.file(chain.file, dest_file)
gunzip(dest_file, remove = F)
chainObject <- rtracklayer::import.chain(gsub(".gz", "", dest_file, fixed = T))

# construct a GRanges object for the rat data
chr <- purrr::map_chr(Rat_CpGs, ~strsplit(.x, ":")[[1]][1])
chr[!grepl("AABR|KL|MT", chr, fixed = F) ] %<>% paste0("chr", .)
pos <- purrr::map_int(Rat_CpGs, ~as.integer(strsplit(.x, ":")[[1]][2]))
Rat_CpGs_Granges <- GRanges(seqnames = Rle(chr), 
                            ranges = IRanges(start = pos, end = pos),
                            rn6_coord = Rat_CpGs)

# liftOver from rn6 to mm10
mm10_df <- as.data.frame(rtracklayer::liftOver(Rat_CpGs_Granges, chainObject)) 
mm10_df %<>%
  dplyr::select(seqnames, start, rn6_coord) %>%
  mutate(mm10_coord = paste0(gsub("chr", "", seqnames), ":", start)) %>%
  dplyr::select(rn6_coord, mm10_coord)
readr::write_tsv(mm10_df, path = "~/Downloads/rn6_mm10_map.tsv")
