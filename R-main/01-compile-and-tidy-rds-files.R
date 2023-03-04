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

#### Read in the label.txt files ####

# Same drill here.  First we make a file that holds a data frame of file names:
# dsb:label_files dianabaetscher$ pwd
# /Users/dianabaetscher/Desktop/NOAA_grad/git-repos/rockfish-species-id/new_baseline_data/label_files
# dsb:label_files dianabaetscher$ ls -l | awk 'BEGIN {print "gtseq_run", "file"} NR >1 && !/posinfo/ {num = $NF; gsub(/gtseq/, "", num); gsub(/[._].*$/, "", num);  print num, $NF}' > ../../new_baseline_data/label-file-list.txt 

# fdf <- read.table("new_baseline_data/label-file-list.txt", stringsAsFactors = FALSE, header = TRUE) %>%
#   tbl_df()
# dir <- "new_baseline_data/label_files"
# 
# # trying to break this down to make it work!
# labels <- lapply(1:nrow(fdf), function(i) {
#   message("Working on ", fdf$file[i]) }) %>%
#   bind_rows()
#   
# # I'm not sure why this isn't working! 
# labels <- lapply(1:nrow(fdf), function(i) {
#   message("Working on ", fdf$file[i])
# read_tsv(path = file.path(dir, fdf$file[i]), col_names = c("file", "id", "species")) %>%
#   select(-file) %>%
#   mutate(gtseq_run = fdf$gtseq_run[i]) %>%
#   select(gtseq_run, id, everything())
# }) %>%
#   bind_rows()
# 
# # Now sure why the function way wasn't working, but here I'll do it manually just to get a move on.
# run11 <- read_tsv("new_baseline_data/label_files/gtseq11_label.txt", col_names = c("file", "id", "species")) %>%
#   select(-file) %>%
#   mutate(gtseq_run = "11") %>%
#   select(gtseq_run, everything()) %>%
#   bind_rows()
# 
# run19 <- read_tsv("new_baseline_data/label_files/gtseq19_label.txt", col_names = c("file", "id", "species")) %>%
#   select(-file) %>%
#   mutate(gtseq_run = "19") %>%
#   select(gtseq_run, everything()) %>%
#   bind_rows()
#   
# run48 <- read_tsv("new_baseline_data/label_files/gtseq48_label.txt", col_names = c("file", "id", "species")) %>%
#   select(-file) %>%
#   mutate(gtseq_run = "48") %>%
#   select(gtseq_run, everything()) %>%
#   bind_rows()
# 
# run54 <- read_tsv("new_baseline_data/label_files/gtseq54_label.txt", col_names = c("file", "id", "species")) %>%
#   select(-file) %>%
#   mutate(gtseq_run = "54") %>%
#   select(gtseq_run, everything()) %>%
#   bind_rows()
# 
# run55 <- read_tsv("new_baseline_data/label_files/gtseq55_label.txt", col_names = c("file", "id", "species")) %>%
#   select(-file) %>%
#   mutate(gtseq_run = "55") %>%
#   select(gtseq_run, everything()) %>%
#   bind_rows()
# 
# run62 <- read_tsv("new_baseline_data/label_files/gtseq62_label.txt", col_names = c("file", "id", "species")) %>%
#   select(-file) %>%
#   mutate(gtseq_run = "62") %>%
#   select(gtseq_run, everything()) %>%
#   bind_rows()
# 
# run65 <- read_tsv("new_baseline_data/label_files/gtseq65_label.txt", col_names = c("file", "id", "species")) %>%
#   select(-file) %>%
#   mutate(gtseq_run = "65") %>%
#   select(gtseq_run, everything()) %>%
#   bind_rows()
# 
# run66 <- read_tsv("new_baseline_data/label_files/gtseq66_label.txt", col_names = c("file", "id", "species")) %>%
#   select(-file) %>%
#   mutate(gtseq_run = "66") %>%
#   select(gtseq_run, everything()) %>%
#   bind_rows()
# 
# # combine all samples
# labels <- bind_rows(run11, run19, run48, run54, run55, run62, run65, run66)
# 
# # save that.
# 
# saveRDS(labels, "new_baseline_data/processed/label-tibble.rds", compress = "xz")
# 
# 
# #### Read the sample sheets from the Excel files ####
# 
# # Same drill here.  First we make a file that holds a data frame of file names:
# # dsb:sample_sheets dianabaetscher$ pwd
# # /Users/dianabaetscher/Desktop/NOAA_grad/git-repos/rockfish-species-id/new_baseline_data/sample_sheets
# # dsb:sample_sheets dianabaetscher$ ls -l | awk 'BEGIN {print "gtseq_run", "file"} NR > 1 {num = $NF; gsub(/GTseq/, "", num); gsub(/_.*$/, "", num);  print num, $NF}' > ../../new_baseline_data/sample-sheet-file-list.txt
# 
# 
# fdf <- read.table("new_baseline_data/sample-sheet-file-list.txt", stringsAsFactors = FALSE, header = TRUE) %>%
#   tbl_df()
# dir <- "new_baseline_data/sample_sheets"
# 
# sample_sheets <- lapply(1:nrow(fdf), function(i) {
#   message("Working on ", fdf$file[i])
#   read_excel(path = file.path(dir, fdf$file[i]), sheet = "sample_sheet", skip = 21) %>%
#     filter(!is.na(Sample_Plate)) %>%  # read_excel sometimes reads a lot of blank lines...
#     tidyr::separate(Sample_Plate, into = c("NMFS_DNA_ID", "ssBOX_ID", "ssBOX_POSITION")) %>%
#     mutate(id = str_replace(Sample_ID, "s_0*", "s")) %>%
#     mutate(gtseq_run = fdf$gtseq_run[i]) %>%
#     select(gtseq_run, id, everything())
# }) %>%
#   bind_rows()
# 
# # save that.
# 
# saveRDS(sample_sheets, "new_baseline_data/processed/sample-sheet-tibble.rds", compress = "xz")
# 
# 
# 
# #### And finally, let's get the meta-data read in and cleaned up (if it needs it) ####
# meta1 <- read_xls("new_baseline_data/metadata/rockfish_spp_metadata.xls") %>%
#   mutate(BATCH_ID = as.character(BATCH_ID),
#          WEIGHT = as.numeric(WEIGHT)) %>%
#   mutate(SAMPLE_ID = as.character(SAMPLE_ID)) %>%
#   mutate(LENGTH = as.numeric(LENGTH)) %>% 
#   mutate(COLLECTION_DATE = ymd(COLLECTION_DATE)) %>%
#   mutate(PICK_DATE = ymd(PICK_DATE))
# 
# # when we read that in, we lose the "None."s in the LEFTOVER_SAMPLE fields.  That 
# # is OK for now.  
# meta2 <- read_xls("new_baseline_data/metadata/R375.xls", sheet = "Repository") %>%
#   mutate(SAMPLE_ID = as.character(SAMPLE_ID)) %>%
# mutate(LENGTH = as.numeric(LENGTH))
# meta3 <- read_xls("new_baseline_data/metadata/R376.xls", sheet = "Repository")%>%
#   mutate(SAMPLE_ID = as.character(SAMPLE_ID))  %>%
#   mutate(LENGTH = as.numeric(LENGTH))
# meta4 <- read_xls("new_baseline_data/metadata/R377.xls", sheet = "Repository")%>%
#   mutate(SAMPLE_ID = as.character(SAMPLE_ID)) %>%
#   mutate(LENGTH = as.numeric(LENGTH))
# 
# more_meta <- bind_rows(meta2, meta3, meta4) %>%
#   mutate(COLLECTION_DATE = ymd(COLLECTION_DATE)) %>% 
#   mutate(PICK_DATE = ymd(PICK_DATE))  %>%
#   mutate(LENGTH = as.numeric(LENGTH))
# 
# # and finally, to get the meta data for the NSF samples from GTseq run 19,
# # we read in the kelp meta data
# 
# meta5 <- readRDS("new_baseline_data/metadata/kelp-meta-data-tibble.rds")
# 
# meta6 <- read_xlsx("new_baseline_data/metadata/Sebastes_metadata_HNuetzel.xlsx") %>%
#   mutate(COLLECTION_DATE = ymd(Collect_date)) %>%
#   select(-Collect_date) %>%
#   mutate(LENGTH = as.numeric(LENGTH))
# 
# meta7 <- read_csv("new_baseline_data/metadata/R370_372.csv") %>%
#   mutate(BATCH_ID = as.character(BATCH_ID),
#          WEIGHT = as.numeric(WEIGHT)) %>%
#   mutate(SAMPLE_ID = as.character(SAMPLE_ID)) %>%
#   mutate(LENGTH = as.numeric(LENGTH)) %>% 
#   mutate(COLLECTION_DATE = ymd(COLLECTION_DATE)) %>%
#   mutate(PICK_DATE = ymd(PICK_DATE)) %>%
#   select(-HAUL)
#   
# most_meta <- bind_rows(meta1, more_meta) %>% 
#   mutate(PICK_DATE = ymd(PICK_DATE))
# 
# meta <- bind_rows(most_meta, meta5, meta6, meta7) 
# 
# saveRDS(meta, "new_baseline_data/processed/meta-data-tibble.rds", compress = "xz")
# 
# # for re-doing this analysis, I just read in the existing sample sheet tibble and meta data tibble
# #meta <- readRDS("data/processed/meta-data-tibble.rds")
# #sample_sheets <- readRDS("data/processed/sample-sheet-tibble.rds")
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

