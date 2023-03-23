library(tidyverse)
library(readxl)
library(stringr)
library(lubridate)


source("R/microhaplot-genos-funcs.R")


#### Call genos from the microhaplot rds files ####

# this next step uses a little data frame that has the names and GTSeq run numbers
# I got that file like this:
# dsb:rds_files dianabaetscher$ pwd
# /Users/dianabaetscher/Desktop/NOAA_grad/git-repos/rockfish-species-id/new_baseline_data/feather_files
# dsb:feather_files dianabaetscher$ ls -l | awk 'BEGIN {print "gtseq_run", "file"} NR >1 && !/posinfo/ {num = $NF; gsub(/gtseq/, "", num); gsub(/[._].*$/, "", num);  print num, $NF}' > ../../data/rds-file-list.txt 

# get the names of the files
fdf <- read.table("data/rds_file_list.txt", stringsAsFactors = FALSE, header = TRUE) %>%
  tbl_df()
dir <- "rds_files"


# cycle over them, read them and add the gtseq_run column on each.
# at the end, bind them together.
genos_long <- lapply(1:nrow(fdf), function(i) {
  message("Working on ", fdf$file[i])
  call_genos_from_haplotRDS(path = file.path(dir, fdf$file[i])) %>%
    mutate(gtseq_run = fdf$gtseq_run[i]) %>%
    select(gtseq_run, everything())
}) %>%
  bind_rows()

# fix the syntax
genos_long$locus <- gsub("_PCR_Product", "", genos_long$locus)
genos_long$locus <- gsub("_extraction", "", genos_long$locus)


# we go ahead and save it in data/processed, with xz compression
saveRDS(genos_long, file = "data/processed/called_genos.rds", compress = "xz")


# 
#### In the end, let us get a data frame that includes genotypes for all the individuals  ####
# and which explicitly has NAs in places where data are missing

genos_long_explicit_NAs <- genos_long %>%
  select(gtseq_run, id) %>%
  unique() %>%
  unite(col = gid, sep = "_", gtseq_run, id) %>%
  select(gid) %>%
  unlist() %>%
  unname() %>%
  expand.grid(gid = ., locus = unique(genos_long$locus), gene_copy = 1:2, stringsAsFactors = FALSE) %>%
  tbl_df() %>%
  separate(gid, into = c("gtseq_run", "id"), convert = TRUE) %>%
  left_join(., genos_long) %>%
  arrange(gtseq_run, id, locus, gene_copy)

# and then save that
saveRDS(genos_long_explicit_NAs, file = "data/processed/called_genos_na_explicit.rds", compress = "xz")

