01-metadata
================
Diana Baetscher
2023-03-30

30 March 2023

Cleaning up final set of metadata before R&R final report and passing
off this repo to Jessie.

## Metadata for bycatch and reference samples

``` r
library(tidyverse)
```

    ## -- Attaching packages --------------------------------------- tidyverse 1.3.1 --

    ## v ggplot2 3.4.1      v purrr   0.3.4 
    ## v tibble  3.1.2      v dplyr   1.0.10
    ## v tidyr   1.2.0      v stringr 1.4.0 
    ## v readr   1.4.0      v forcats 0.5.1

    ## Warning: package 'ggplot2' was built under R version 4.1.3

    ## Warning: package 'tidyr' was built under R version 4.1.3

    ## Warning: package 'dplyr' was built under R version 4.1.3

    ## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(data.table)
```

    ## 
    ## Attaching package: 'data.table'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last

    ## The following object is masked from 'package:purrr':
    ## 
    ##     transpose

``` r
library(readxl)
```

Metadata from Jessie for the bycatch and banded bycatch birds

``` r
# while I'm at it, let's clean up location info for this
# add metadata to that based on sample id
new_meta <- read_xlsx("../data/metadata/Oikonos_NOAA_BFAL_Metadata_updated 3.15.23.xlsx")


# all possible location info syntax
meta_fixed_locs <- new_meta %>%
  #select(Location, `Loc-Abbr`) %>%
  #unique() %>%
   mutate(Location = ifelse(Location == "French Frigate Shoals (Tern Island)", "Tern Island, French Frigate Shoals, HI", Location)) %>%
    mutate(Location = ifelse(Location == "Laysan Island", "Laysan Island, Honolulu County, Hawaii", Location)) %>%
      mutate(Location = ifelse(Location == "Kure Attol (Green Island)", "Green Island, Kure Atoll, Honolulu County, Hawaii", Location)) %>%
       mutate(`Loc-Abbr` = ifelse(Location == "Tern Island, French Frigate Shoals, HI", "NWHI, FFS", `Loc-Abbr`)) %>%
        mutate(`Loc-Abbr` = ifelse(Location == "Laysan Island, Honolulu County, Hawaii", "NWHI, Laysan", `Loc-Abbr`)) %>%
          mutate(`Loc-Abbr` = ifelse(Location == "French Frigate Shoals", "NWHI, FFS", `Loc-Abbr`)) %>%
            mutate(`Loc-Abbr` = ifelse(Location == "Green Island, Kure Atoll, Honolulu County, Hawaii", "NWHI, Kure", `Loc-Abbr`)) %>%
   unique() #%>% # there are three metadata entries that are exact duplicates. Remove those here
  # arrange(Location)
```

846 bycatch, some of which were banded and used as reference birds

Metadata from the reference birds (non-banded)

``` r
# from the MiSeq sample sheet
ref_samplelist <- read_csv("../data/gtseq5_samplelist.csv") %>%
  mutate(gtseq_run = "gtseq5") %>%
  rename(id = sample)
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   Sample_ID = col_character(),
    ##   Sample_Plate = col_character(),
    ##   Sample_Well = col_character(),
    ##   population = col_character(),
    ##   sample = col_character()
    ## )

``` r
# metadata for reference samples used for lcWGS
ref_pops <- readxl::read_xlsx("../data/BFAL_WGS_01_finalPlateMap.xlsx", sheet = "combinedPlate1And2SamplePops")


reference_meta_wo_banded <- ref_pops %>%
  left_join(., ref_samplelist, by = c("ind_ID" = "Sample_ID")) %>%
  filter(str_detect(ind_ID, "BFAL-")) %>% # just the reference samples that came from non-Jessie sources (not banded birds)
  rename(Location = population) %>%
  mutate(Location = ifelse(Location == "Kure", "Kure Atoll, Honolulu County, Hawaii", Location)) %>%
   mutate(Location = ifelse(Location == "Midway", "Midway Atoll, NWHI, HI", Location)) %>%
      mutate(Location = ifelse(Location == "Torishima", "Torishima Island, Japan", Location)) %>%
       mutate(`Loc-Abbr` = ifelse(Location == "Midway Atoll, NWHI, HI", "NWHI, Midway", NA)) %>%
        mutate(`Loc-Abbr` = ifelse(Location == "Torishima Island, Japan", "Japan, Torishima", `Loc-Abbr`)) %>%
            mutate(`Loc-Abbr` = ifelse(Location == "Kure Atoll, Honolulu County, Hawaii", "NWHI, Kure", `Loc-Abbr`)) %>%
  rename(`MWVCRC#` = ind_ID)
  
reference_meta_wo_banded
```

    ## # A tibble: 92 x 8
    ##    `MWVCRC#`    pop   Sample_Plate Sample_Well Location    id    gtseq~1 Loc-A~2
    ##    <chr>        <chr> <chr>        <chr>       <chr>       <chr> <chr>   <chr>  
    ##  1 BFAL-KURE-01 Kure  BFAL002      A07         Kure Atoll~ s337  gtseq5  NWHI, ~
    ##  2 BFAL-KURE-02 Kure  BFAL002      B07         Kure Atoll~ s338  gtseq5  NWHI, ~
    ##  3 BFAL-KURE-03 Kure  BFAL002      C07         Kure Atoll~ s339  gtseq5  NWHI, ~
    ##  4 BFAL-KURE-04 Kure  BFAL002      D07         Kure Atoll~ s340  gtseq5  NWHI, ~
    ##  5 BFAL-KURE-05 Kure  BFAL002      E07         Kure Atoll~ s341  gtseq5  NWHI, ~
    ##  6 BFAL-KURE-06 Kure  BFAL002      F07         Kure Atoll~ s342  gtseq5  NWHI, ~
    ##  7 BFAL-KURE-07 Kure  BFAL002      G07         Kure Atoll~ s343  gtseq5  NWHI, ~
    ##  8 BFAL-KURE-08 Kure  BFAL002      H07         Kure Atoll~ s344  gtseq5  NWHI, ~
    ##  9 BFAL-KURE-09 Kure  BFAL002      A08         Kure Atoll~ s345  gtseq5  NWHI, ~
    ## 10 BFAL-KURE-10 Kure  BFAL002      B08         Kure Atoll~ s346  gtseq5  NWHI, ~
    ## # ... with 82 more rows, and abbreviated variable names 1: gtseq_run,
    ## #   2: `Loc-Abbr`

92 reference samples without banded birds.

``` r
full_metadata_all_samples <- reference_meta_wo_banded %>%
  full_join(., meta_fixed_locs)  %>%
  select(-Sample_Plate, -Sample_Well, -id, -gtseq_run) %>%
  rename(sampleID = `MWVCRC#`, population = pop)
```

    ## Joining, by = c("MWVCRC#", "Location", "Loc-Abbr")

``` r
full_metadata_all_samples %>%
  write_rds("../data/processed/metadata_bycatch_and_reference_20230330.rds")
```

## Samplesheets

Now organize the sample info associated with the birds we genotyped at
ABL

``` r
#create a list of the files from your target directory
file_list <- list.files(path="../data/samplesheets/")

#initiate a blank data frame, each iteration of the loop will append the data from the given file to this variable
samplesheets <- data.frame()

# read in sample sheets
for (i in 1:length(file_list)){
  temp_data <- read_csv(paste0("../data/samplesheets/",file_list[i]), skip = 19) #read in files using the fread function from the data.table package
  samplesheets <- rbindlist(list(samplesheets, temp_data), use.names = T) #for each iteration, bind the new data to the building dataset
}
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   Sample_ID = col_character(),
    ##   Sample_Plate = col_character(),
    ##   Sample_Well = col_character(),
    ##   I7_Index_ID = col_character(),
    ##   index = col_character(),
    ##   I5_Index_ID = col_character(),
    ##   index2 = col_character(),
    ##   Sample_Project = col_character(),
    ##   Description = col_character(),
    ##   s_number = col_character()
    ## )
    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   Sample_ID = col_character(),
    ##   Sample_Plate = col_character(),
    ##   Sample_Well = col_character(),
    ##   I7_Index_ID = col_character(),
    ##   index = col_character(),
    ##   I5_Index_ID = col_character(),
    ##   index2 = col_character(),
    ##   Sample_Project = col_character(),
    ##   Description = col_character(),
    ##   s_number = col_character()
    ## )
    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   Sample_ID = col_character(),
    ##   Sample_Plate = col_character(),
    ##   Sample_Well = col_character(),
    ##   I7_Index_ID = col_character(),
    ##   index = col_character(),
    ##   I5_Index_ID = col_character(),
    ##   index2 = col_character(),
    ##   Sample_Project = col_character(),
    ##   Description = col_character(),
    ##   s_number = col_character()
    ## )

``` r
samplesheets1 <- samplesheets %>%
  rename(id = s_number, gtseq_run = Description) 



samplesheets1$Sample_ID <- gsub("_", "-", samplesheets1$Sample_ID)
samplesheets1$Sample_ID <- gsub("banded", "", samplesheets1$Sample_ID)

samplesheets2 <- samplesheets1 %>%  # one ID to fix based on Claire's description of mislabeling
  mutate(Sample_ID = ifelse(Sample_ID == "BFAL08-1144", "BFAL08-1844", Sample_ID))
```
