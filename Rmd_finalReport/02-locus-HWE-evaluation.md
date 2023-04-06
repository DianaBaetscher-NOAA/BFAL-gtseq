02-locus-HWE-evaluation
================
Diana Baetscher
2023-03-30

Using the complete dataset to evaluate HWE across loci and populations
in the baseline.

Updates: use the NAs explicit genos file. Don’t filter loci based on
missing data at this point. Just look at reference pops for evaluating
HWE.

Evaluating the loci that we’re using for self-assignment and bycatch
assignment. These were originally derived from lcWGS data and the
microhaplotype results suggest that the real variation present is
potentially quite different. Initially, I should probably just focus on
the baseline genotypes rather than including everything.

To evaluate the markers, I need to subset the reference baseline
genotypes for just the minimal amount of missing data.

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
library(adegenet)
```

    ## Loading required package: ade4

    ## Registered S3 method overwritten by 'spdep':
    ##   method   from
    ##   plot.mst ape

    ## 
    ##    /// adegenet 2.1.3 is loaded ////////////
    ## 
    ##    > overview: '?adegenet'
    ##    > tutorials/doc/questions: 'adegenetWeb()' 
    ##    > bug reports/feature requests: adegenetIssues()

``` r
library(radiator)
library(DescTools)
```

    ## Warning: package 'DescTools' was built under R version 4.1.3

    ## Registered S3 method overwritten by 'DescTools':
    ##   method         from 
    ##   reorder.factor gdata

``` r
# read in rds file with genotypes
genos_long <- read_rds("../data/processed/called_genos_na_explicit.rds")

# just select the reference samples based on metadata
samplesheets <- read_rds("../data/processed/metadata_bycatch_and_reference_20230330.rds")

reference_samples <- samplesheets %>%
  filter(Sample_Plate %in% c("BFAL001", "BFAL002") &
           !is.na(Location),
         !`Loc-Abbr` %in% c("Oahu", "Kauai, Kilauea")) %>%
  arrange(sampleID) %>%
  select(sampleID, gtseq_run, id, Location, `Loc-Abbr`)

# combine the genotypes with the reference
ref_genos <- reference_samples %>%
  left_join(., genos_long, by = c("gtseq_run", "id")) 
```

Deal explicitly with duplicate samples:

## Some initial filters

### Take highest read-depth call for multiply-genotyped DNA_IDs

I’m not sure if there are any of these, but best to leave it in here…

Now, here is a harder operation: if an individual is multiply-genotyped,
take the genotype with the highest total read depth.

``` r
# slow-ish function to get the total read depth column
tdepth <- function(a, d) {
  if(any(is.na(a))) {
    return(NA)
  }
  if(a[1]==a[2]) {
    return(d[1])
  } else {
    return(d[1] + d[2])
  }
  
}
# this takes the highest read-depth instance of each duplicately-genotyped individual.
geno_one_each <- ref_genos %>%
  group_by(sampleID, locus, gtseq_run) %>%
  mutate(total_depth = tdepth(allele, depth)) %>%
  ungroup() %>%
  arrange(sampleID, locus, total_depth, gtseq_run, depth) %>%
  group_by(sampleID, locus) %>%
  mutate(rank = 1:n()) %>% # this ranks all alleles for a given individual by depth
  ungroup() %>%
  filter(rank <= 2) # this takes the top ranked alleles and should remove duplicates
  

# how many samples now?
geno_one_each %>%
  filter(!is.na(sampleID)) %>%
  group_by(sampleID) %>% 
  select(sampleID, gtseq_run, id) %>%
  unique() %>%
  tally() 
```

    ## # A tibble: 151 x 2
    ##    sampleID     n
    ##    <chr>    <int>
    ##  1 12-0155      1
    ##  2 12-0270      1
    ##  3 12-0271      1
    ##  4 12-0371      1
    ##  5 12-0372      1
    ##  6 12-0381      1
    ##  7 13-0991      1
    ##  8 13-1187      1
    ##  9 13-1198      1
    ## 10 14-0454      1
    ## # ... with 141 more rows

Take a look at missing data overall:

``` r
geno_one_each %>%
  ggplot(aes(x = reorder(locus, depth), y = reorder(sampleID, depth), fill = log10(depth))) +
  geom_tile()
```

![](02-locus-HWE-evaluation_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

Deal with missing data:

## Missing data in loci

How many loci and how many alleles?

``` r
# alleles
geno_one_each %>%
  filter(!is.na(allele)) %>%
  select(locus, allele) %>%
  unique()
```

    ## # A tibble: 458 x 2
    ##    locus                allele
    ##    <chr>                <chr> 
    ##  1 scaffold_0_4543438   GGG   
    ##  2 scaffold_0_4543438   GGA   
    ##  3 scaffold_1_5556176   CT    
    ##  4 scaffold_10_1333563  G     
    ##  5 scaffold_10_1333563  A     
    ##  6 scaffold_1050_36334  T     
    ##  7 scaffold_1053_177405 A     
    ##  8 scaffold_1078_116737 GC    
    ##  9 scaffold_11_3057945  C     
    ## 10 scaffold_115_2200290 C     
    ## # ... with 448 more rows

``` r
# loci
geno_one_each %>%
  filter(!is.na(allele)) %>%
  select(locus, allele) %>%
  unique() %>%
  group_by(locus) %>%
  tally() %>%
  arrange(desc(n))
```

    ## # A tibble: 189 x 2
    ##    locus                    n
    ##    <chr>                <int>
    ##  1 scaffold_1226_115783     6
    ##  2 scaffold_1589_18141      5
    ##  3 scaffold_316_437138      5
    ##  4 scaffold_33_2811656      5
    ##  5 scaffold_48_3128896      5
    ##  6 scaffold_75_1184494      5
    ##  7 scaffold_0_4543438       4
    ##  8 scaffold_1164_119437     4
    ##  9 scaffold_127_1105814     4
    ## 10 scaffold_15_2923659      4
    ## # ... with 179 more rows

In the reference, 458 alleles across 189 loci with 1-6 alleles per
locus.

Some individuals are missing all their data, and some loci are missing
all their data. Beginning with loci:

``` r
# missing data across loci
locs_to_toss <- geno_one_each %>%
  group_by(locus) %>%
  mutate(missingness = ifelse(is.na(allele), 1, 0)) %>%
  summarise(sum(missingness)) %>% # 189 loci x2 = total
  filter(`sum(missingness)`>151) %>% # more than 50% missing data (given 151 individuals)
  select(locus) # drop those loci for now and see how the assignment goes

# just the keepers
genos_locs_filtered <- geno_one_each %>%
  anti_join(., locs_to_toss)
```

    ## Joining, by = "locus"

Remove these 4 loci because of too much missing data in the reference
samples.

## Missing data in individuals

Total number of loci = 185

Total number of samples = 151

``` r
inds_to_toss <- genos_locs_filtered %>%
  group_by(sampleID) %>%
  mutate(missingness = ifelse(is.na(allele), 1, 0)) %>%
  summarise(sum(missingness)) %>%
  arrange(desc(`sum(missingness)`)) %>%
  #filter(`sum(missingness)` > 47) # remove samples with >25% missing data
 filter(`sum(missingness)` > 19) # remove samples with >10% missing data

# just the keepers
genos_locs_ind_filtered <- genos_locs_filtered %>%
  anti_join(., inds_to_toss)
```

    ## Joining, by = "sampleID"

23 samples removed with \>10% missing data.

Now look at missing data:

``` r
genos_locs_ind_filtered %>%
  ggplot(aes(x = reorder(locus, -depth), y = reorder(sampleID, -depth), fill = log10(depth))) +
  geom_tile()
```

![](02-locus-HWE-evaluation_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Much better. Now move forward from here.

## Prep data for analysis

``` r
# first make integers of the alleles
alle_idxs <- genos_locs_ind_filtered %>% 
  filter(!is.na(sampleID)) %>%
  dplyr::select(sampleID, locus, gene_copy, allele) %>%
  group_by(locus) %>%
  mutate(alleidx = as.integer(factor(allele, levels = unique(allele)))) %>%
  ungroup() %>%
  arrange(sampleID, locus, alleidx) # rubias can handle NA's, so no need to change them to 0's

# reformat
reference <- alle_idxs %>%
  inner_join(., reference_samples) %>%
  select(-allele, -gtseq_run, -id) %>%
  select(`Loc-Abbr`, sampleID, everything()) %>%
  rename(collection = `Loc-Abbr`, indiv = sampleID)
```

    ## Joining, by = "sampleID"

I need to finagle the dataset into a genind object.

``` r
# read in data
reference_genos <- reference

# number of missing alleles (max = 185*2)
inds_to_toss_missing_data <- reference_genos %>%
  group_by(indiv) %>%
  filter(is.na(alleidx)) %>%
  group_by(indiv, collection) %>%
  tally() %>%
  arrange(desc(n)) %>%
  filter(n > 9)
  
# if I removed all the indivs with > 10 missing alleles, that's only 5 indivs.
# what would the distribution of indivs per population be? 
# still ok, I'm pretty sure.
reference_genos %>%
  anti_join(., inds_to_toss_missing_data) %>%
  group_by(collection, indiv) %>%
  tally() %>%
  tally()
```

    ## Joining, by = c("collection", "indiv")

    ## # A tibble: 6 x 2
    ##   collection           n
    ##   <chr>            <int>
    ## 1 Japan, Torishima    40
    ## 2 Kauai, Lehua         2
    ## 3 NWHI, FFS           35
    ## 4 NWHI, Kure          22
    ## 5 NWHI, Laysan         8
    ## 6 NWHI, Midway         7

I feel fine about that. Let’s see if that makes any of the HWE
calculations easier because less missing data?

``` r
# SKIP THIS FOR NOW.
# ref_genos_low_missing_data <- reference_genos %>%
#   anti_join(., inds_to_toss_missing_data)

ref_genos_low_missing_data <- reference_genos
```

Make the df match the requirements for tidy_genomic_data

``` r
long_df <- ref_genos_low_missing_data  %>%
  select(-gene_copy) %>%
  select(collection, everything()) %>%
  rename(INDIVIDUALS = indiv, STRATA = collection, MARKERS = locus, GT = alleidx)
```

Genotypes should be coded with 3 integers for each alleles. 6 integers
in total for the genotypes. e.g. 001002 or 111333 (for heterozygote
individual). 6 integers WITH separator: e.g. 001/002 or 111/333 (for
heterozygote individual). The separator can be any of these: “/”, “:”,
“\_“,”-“,”.”, and will be removed.

``` r
# create 3 digit integers from the genotypes
long_df$GT3 <- Format(long_df$GT, ldigits = 3, digits = 0)

# fix NAs
long_df0s <- long_df %>%
  mutate(GT3 = ifelse(is.na(GT3), "000", GT3)) # I don't love that this creates potential artifacts!
```

Now combine the GT3 column per indiv/marker:

``` r
# make the genos characters and then try pasting them as strings
long_df0s$GT3 <- as.character(long_df0s$GT3)

long_df3digit <- long_df0s %>%
  group_by(INDIVIDUALS, MARKERS) %>% 
  arrange(GT3, .by_group = TRUE) %>% 
  summarise(GENOTYPE = toString(GT3))
```

    ## `summarise()` has grouped output by 'INDIVIDUALS'. You can override using the
    ## `.groups` argument.

``` r
# paste strings together
long_df3digit$GENOTYPE <- gsub(", ","",long_df3digit$GENOTYPE)


# add back on species identity as strata
df_for_conversion <- long_df0s %>% 
  select(-GT, -GT3) %>%
  left_join(., long_df3digit) %>%
  unique() %>%
  rename(GT = GENOTYPE) %>%
  mutate(GT = ifelse(GT == "000000", NA, GT)) %>%
  filter(!str_detect(STRATA, "Kauai"))
```

    ## Joining, by = c("INDIVIDUALS", "MARKERS")

``` r
df_for_conversion$STRATA <- as.factor(df_for_conversion$STRATA)
```

Double check - how many samples per population?

``` r
df_for_conversion %>%
  select(INDIVIDUALS, STRATA) %>%
  unique() %>%
  group_by(STRATA) %>%
  tally()
```

    ## # A tibble: 5 x 2
    ##   STRATA               n
    ##   <fct>            <int>
    ## 1 Japan, Torishima    43
    ## 2 NWHI, FFS           36
    ## 3 NWHI, Kure          24
    ## 4 NWHI, Laysan         9
    ## 5 NWHI, Midway        14

Let’s remove Lehua and Kilauea at this point because of so few samples.

``` r
convert_df_wo_2pops <- df_for_conversion 
```

``` r
# use the radiator package for this conversion
genind_df <- write_genind(convert_df_wo_2pops)
```

Basic population genetic evaluation of the markers, following, in part,
this: <https://bookdown.org/hhwagner1/LandGenCourse_book/WE_3.html>

``` r
sum_output <- summary(genind_df)

names(sum_output)
```

    ## [1] "n"         "n.by.pop"  "loc.n.all" "pop.n.all" "NA.perc"   "Hobs"     
    ## [7] "Hexp"

``` r
expected_hz <- as.data.frame(sum_output$Hexp)
observed_hz <- as.data.frame(sum_output$Hobs)
```

Expected heterozygosity (here: Hexp) is the heterozygosity expected in a
population under HWE, and observed heterozygosity (here: Hobs) is the
observed number of heterozygotes at a locus divided by the total number
of genotyped individuals. Here are the global values (pooled across all
populations):

``` r
# expected heterozygosity per population
adegenet::Hs(genind2genpop(genind_df))
```

    ## 
    ##  Converting data from a genind to a genpop object... 
    ## 
    ## ...done.

    ## Japan,_Torishima        NWHI,_FFS       NWHI,_Kure     NWHI,_Laysan 
    ##        0.3155225        0.3627827        0.3531330        0.3728518 
    ##     NWHI,_Midway 
    ##        0.3470202

``` r
Hobs <- t(sapply(seppop(genind_df), function(ls) summary(ls)$Hobs))
  Hexp <- t(sapply(seppop(genind_df), function(ls) summary(ls)$Hexp))
  {cat("Expected heterozygosity (Hexp):", "\n")
  round(Hexp, 2)
  cat("\n", "Observed heterozygosity (Hobs):", "\n")
  round(Hobs, 2)}
```

    ## Expected heterozygosity (Hexp): 
    ## 
    ##  Observed heterozygosity (Hobs):

    ##                  scaffold_0_4543438 scaffold_102_694013 scaffold_1050_36334
    ## Japan,_Torishima               0.28                0.19                0.14
    ## NWHI,_FFS                      0.53                0.56                0.28
    ## NWHI,_Kure                     0.58                0.46                0.46
    ## NWHI,_Laysan                   0.78                0.44                0.44
    ## NWHI,_Midway                   0.50                0.36                0.29
    ##                  scaffold_1053_177405 scaffold_1078_116737 scaffold_10_1333563
    ## Japan,_Torishima                 0.37                 0.49                0.49
    ## NWHI,_FFS                        0.31                 0.56                0.36
    ## NWHI,_Kure                       0.50                 0.46                0.50
    ## NWHI,_Laysan                     0.33                 0.56                0.44
    ## NWHI,_Midway                     0.50                 0.50                0.43
    ##                  scaffold_115_2200290 scaffold_1164_119437 scaffold_116_1932131
    ## Japan,_Torishima                    0                 0.02                 0.33
    ## NWHI,_FFS                           0                 0.06                 0.28
    ## NWHI,_Kure                          0                 0.12                 0.50
    ## NWHI,_Laysan                        0                 0.22                 0.22
    ## NWHI,_Midway                        0                 0.14                 0.50
    ##                  scaffold_116_386790 scaffold_11_3057945 scaffold_121_1829963
    ## Japan,_Torishima                0.05                   0                    0
    ## NWHI,_FFS                       0.25                   0                    0
    ## NWHI,_Kure                      0.12                   0                    0
    ## NWHI,_Laysan                    0.11                   0                    0
    ## NWHI,_Midway                    0.14                   0                    0
    ##                  scaffold_121_565678 scaffold_1226_115783 scaffold_123_1115977
    ## Japan,_Torishima                0.47                 0.09                 0.72
    ## NWHI,_FFS                       0.50                 0.00                 0.72
    ## NWHI,_Kure                      0.54                 0.04                 0.71
    ## NWHI,_Laysan                    0.33                 0.11                 0.67
    ## NWHI,_Midway                    0.29                 0.79                 0.43
    ##                  scaffold_1273_1937 scaffold_127_1105814 scaffold_127_901882
    ## Japan,_Torishima               0.09                 0.51                   0
    ## NWHI,_FFS                      0.06                 0.53                   0
    ## NWHI,_Kure                     0.00                 0.58                   0
    ## NWHI,_Laysan                   0.00                 0.33                   0
    ## NWHI,_Midway                   0.00                 0.57                   0
    ##                  scaffold_12_5002902 scaffold_1327_103421 scaffold_140_259881
    ## Japan,_Torishima                0.47                 0.09                0.30
    ## NWHI,_FFS                       0.19                 0.31                0.42
    ## NWHI,_Kure                      0.46                 0.29                0.33
    ## NWHI,_Laysan                    0.33                 0.44                0.22
    ## NWHI,_Midway                    0.29                 0.21                0.43
    ##                  scaffold_143_1756997 scaffold_146_16700 scaffold_148_727111
    ## Japan,_Torishima                 0.00               0.35                0.51
    ## NWHI,_FFS                        0.08               0.44                0.47
    ## NWHI,_Kure                       0.08               0.50                0.50
    ## NWHI,_Laysan                     0.22               0.44                0.22
    ## NWHI,_Midway                     0.00               0.50                0.43
    ##                  scaffold_14_2621198 scaffold_14_5075430 scaffold_155_1707144
    ## Japan,_Torishima                0.21                0.40                 0.40
    ## NWHI,_FFS                       0.25                0.50                 0.25
    ## NWHI,_Kure                      0.42                0.54                 0.21
    ## NWHI,_Laysan                    0.44                0.44                 0.22
    ## NWHI,_Midway                    0.50                0.36                 0.29
    ##                  scaffold_155_468225 scaffold_157_815403 scaffold_1589_18141
    ## Japan,_Torishima                0.28                0.28                0.53
    ## NWHI,_FFS                       0.28                0.22                0.44
    ## NWHI,_Kure                      0.25                0.38                0.38
    ## NWHI,_Laysan                    0.22                0.33                0.56
    ## NWHI,_Midway                    0.36                0.21                0.43
    ##                  scaffold_15_2923659 scaffold_15_5063205 scaffold_164_1134575
    ## Japan,_Torishima                0.30                0.23                 0.14
    ## NWHI,_FFS                       0.19                0.28                 0.06
    ## NWHI,_Kure                      0.62                0.42                 0.21
    ## NWHI,_Laysan                    0.56                0.56                 0.22
    ## NWHI,_Midway                    0.29                0.50                 0.50
    ##                  scaffold_166_410622 scaffold_16_1050955 scaffold_16_1254862
    ## Japan,_Torishima                0.02                0.19                0.37
    ## NWHI,_FFS                       0.14                0.25                0.36
    ## NWHI,_Kure                      0.12                0.29                0.38
    ## NWHI,_Laysan                    0.33                0.00                0.44
    ## NWHI,_Midway                    0.07                0.29                0.29
    ##                  scaffold_16_31673 scaffold_16_32172 scaffold_177_673090
    ## Japan,_Torishima              0.44              0.44                0.12
    ## NWHI,_FFS                     0.50              0.33                0.44
    ## NWHI,_Kure                    0.42              0.42                0.33
    ## NWHI,_Laysan                  0.44              0.56                0.22
    ## NWHI,_Midway                  0.57              0.50                0.43
    ##                  scaffold_184_724429 scaffold_184_734991 scaffold_185_1133045
    ## Japan,_Torishima                0.53                0.72                 0.35
    ## NWHI,_FFS                       0.56                0.58                 0.19
    ## NWHI,_Kure                      0.67                0.67                 0.46
    ## NWHI,_Laysan                    0.44                0.44                 0.33
    ## NWHI,_Midway                    0.64                0.57                 0.50
    ##                  scaffold_185_1154507 scaffold_190_1605668 scaffold_199_875998
    ## Japan,_Torishima                 0.35                 0.14                0.58
    ## NWHI,_FFS                        0.25                 0.22                0.44
    ## NWHI,_Kure                       0.46                 0.33                0.38
    ## NWHI,_Laysan                     0.33                 0.33                0.67
    ## NWHI,_Midway                     0.64                 0.14                0.57
    ##                  scaffold_1_5556176 scaffold_204_1432955 scaffold_204_239685
    ## Japan,_Torishima               0.30                 0.19                0.23
    ## NWHI,_FFS                      0.50                 0.28                0.47
    ## NWHI,_Kure                     0.42                 0.12                0.33
    ## NWHI,_Laysan                   0.56                 0.22                0.22
    ## NWHI,_Midway                   0.29                 0.29                0.50
    ##                  scaffold_209_721065 scaffold_20_1133858 scaffold_210_1478805
    ## Japan,_Torishima                0.28                0.47                 0.09
    ## NWHI,_FFS                       0.47                0.39                 0.39
    ## NWHI,_Kure                      0.33                0.38                 0.25
    ## NWHI,_Laysan                    0.22                0.11                 0.11
    ## NWHI,_Midway                    0.29                0.36                 0.21
    ##                  scaffold_214_606303 scaffold_223_94277 scaffold_224_319624
    ## Japan,_Torishima                0.07                  0                0.19
    ## NWHI,_FFS                       0.25                  0                0.39
    ## NWHI,_Kure                      0.21                  0                0.21
    ## NWHI,_Laysan                    0.33                  0                0.44
    ## NWHI,_Midway                    0.43                  0                0.14
    ##                  scaffold_227_1219626 scaffold_229_770334 scaffold_234_1146621
    ## Japan,_Torishima                 0.21                0.70                 0.00
    ## NWHI,_FFS                        0.11                0.64                 0.17
    ## NWHI,_Kure                       0.08                0.54                 0.12
    ## NWHI,_Laysan                     0.00                0.56                 0.11
    ## NWHI,_Midway                     0.00                0.57                 0.29
    ##                  scaffold_237_341989 scaffold_238_668548 scaffold_245_674441
    ## Japan,_Torishima                0.19                0.58                   0
    ## NWHI,_FFS                       0.31                0.31                   0
    ## NWHI,_Kure                      0.21                0.50                   0
    ## NWHI,_Laysan                    0.22                0.33                   0
    ## NWHI,_Midway                    0.14                0.36                   0
    ##                  scaffold_246_1244055 scaffold_246_31712 scaffold_247_950885
    ## Japan,_Torishima                 0.40               0.21                0.58
    ## NWHI,_FFS                        0.33               0.69                0.53
    ## NWHI,_Kure                       0.50               0.38                0.54
    ## NWHI,_Laysan                     0.33               0.56                0.44
    ## NWHI,_Midway                     0.64               0.50                0.50
    ##                  scaffold_249_340732 scaffold_24_1490422 scaffold_253_684619
    ## Japan,_Torishima                0.09                0.30                0.23
    ## NWHI,_FFS                       0.19                0.11                0.44
    ## NWHI,_Kure                      0.12                0.29                0.33
    ## NWHI,_Laysan                    0.11                0.22                0.22
    ## NWHI,_Midway                    0.14                0.29                0.57
    ##                  scaffold_256_454799 scaffold_275_1006390 scaffold_27_1646617
    ## Japan,_Torishima                0.40                 0.51                0.65
    ## NWHI,_FFS                       0.47                 0.44                0.50
    ## NWHI,_Kure                      0.25                 0.42                0.46
    ## NWHI,_Laysan                    0.78                 0.56                0.56
    ## NWHI,_Midway                    0.64                 0.50                0.71
    ##                  scaffold_284_808709 scaffold_286_314167 scaffold_287_423519
    ## Japan,_Torishima                0.21                0.26                0.28
    ## NWHI,_FFS                       0.17                0.53                0.50
    ## NWHI,_Kure                      0.12                0.46                0.25
    ## NWHI,_Laysan                    0.11                0.44                0.44
    ## NWHI,_Midway                    0.29                0.36                0.50
    ##                  scaffold_28_4260845 scaffold_28_4262705 scaffold_298_359540
    ## Japan,_Torishima                0.44                0.51                0.26
    ## NWHI,_FFS                       0.33                0.33                0.28
    ## NWHI,_Kure                      0.50                0.54                0.12
    ## NWHI,_Laysan                    0.44                0.44                0.22
    ## NWHI,_Midway                    0.50                0.43                0.07
    ##                  scaffold_298_460712 scaffold_29_3347800 scaffold_306_328368
    ## Japan,_Torishima                0.05                0.44                0.63
    ## NWHI,_FFS                       0.22                0.53                0.50
    ## NWHI,_Kure                      0.17                0.50                0.67
    ## NWHI,_Laysan                    0.00                0.22                0.22
    ## NWHI,_Midway                    0.14                0.21                0.43
    ##                  scaffold_306_944429 scaffold_310_880276 scaffold_316_437138
    ## Japan,_Torishima                0.30                0.19                0.67
    ## NWHI,_FFS                       0.36                0.22                0.64
    ## NWHI,_Kure                      0.29                0.00                0.33
    ## NWHI,_Laysan                    0.22                0.11                0.56
    ## NWHI,_Midway                    0.57                0.07                0.71
    ##                  scaffold_31_750839 scaffold_324_135439 scaffold_32_3673235
    ## Japan,_Torishima               0.47                0.63                0.35
    ## NWHI,_FFS                      0.47                0.44                0.42
    ## NWHI,_Kure                     0.33                0.42                0.38
    ## NWHI,_Laysan                   0.33                0.00                0.33
    ## NWHI,_Midway                   0.43                0.43                0.29
    ##                  scaffold_32_803437 scaffold_32_811415 scaffold_335_662765
    ## Japan,_Torishima               0.05               0.05                0.23
    ## NWHI,_FFS                      0.47               0.42                0.25
    ## NWHI,_Kure                     0.42               0.42                0.17
    ## NWHI,_Laysan                   0.22               0.33                0.00
    ## NWHI,_Midway                   0.43               0.21                0.07
    ##                  scaffold_336_143440 scaffold_33_2797456 scaffold_33_2811656
    ## Japan,_Torishima                0.42                0.40                0.28
    ## NWHI,_FFS                       0.31                0.31                0.39
    ## NWHI,_Kure                      0.29                0.67                0.21
    ## NWHI,_Laysan                    0.56                0.33                0.44
    ## NWHI,_Midway                    0.14                0.50                0.64
    ##                  scaffold_342_417561 scaffold_343_778064 scaffold_343_863766
    ## Japan,_Torishima                0.09                0.30                0.00
    ## NWHI,_FFS                       0.33                0.39                0.00
    ## NWHI,_Kure                      0.12                0.33                0.00
    ## NWHI,_Laysan                    0.11                0.56                0.11
    ## NWHI,_Midway                    0.21                0.14                0.00
    ##                  scaffold_345_694649 scaffold_34_2385714 scaffold_351_703875
    ## Japan,_Torishima                0.16                0.42                0.40
    ## NWHI,_FFS                       0.42                0.33                0.53
    ## NWHI,_Kure                      0.42                0.33                0.54
    ## NWHI,_Laysan                    0.22                0.44                0.44
    ## NWHI,_Midway                    0.57                0.29                0.43
    ##                  scaffold_356_86112 scaffold_360_269450 scaffold_369_134745
    ## Japan,_Torishima               0.33                   0                0.23
    ## NWHI,_FFS                      0.17                   0                0.25
    ## NWHI,_Kure                     0.29                   0                0.12
    ## NWHI,_Laysan                   0.22                   0                0.44
    ## NWHI,_Midway                   0.00                   0                0.14
    ##                  scaffold_381_698673 scaffold_397_13764 scaffold_40_612943
    ## Japan,_Torishima                0.35               0.28               0.63
    ## NWHI,_FFS                       0.19               0.36               0.50
    ## NWHI,_Kure                      0.25               0.46               0.67
    ## NWHI,_Laysan                    0.00               0.33               0.56
    ## NWHI,_Midway                    0.21               0.64               0.64
    ##                  scaffold_411_530837 scaffold_417_428513 scaffold_41_907611
    ## Japan,_Torishima                0.51                0.00               0.26
    ## NWHI,_FFS                       0.42                0.31               0.50
    ## NWHI,_Kure                      0.29                0.17               0.62
    ## NWHI,_Laysan                    0.67                0.44               0.44
    ## NWHI,_Midway                    0.36                0.29               0.43
    ##                  scaffold_429_778790 scaffold_437_192045 scaffold_45_3079339
    ## Japan,_Torishima                0.16                0.07                0.40
    ## NWHI,_FFS                       0.44                0.33                0.69
    ## NWHI,_Kure                      0.42                0.54                0.62
    ## NWHI,_Laysan                    0.44                0.22                0.44
    ## NWHI,_Midway                    0.29                0.36                0.43
    ##                  scaffold_460_237213 scaffold_461_674999 scaffold_472_10182
    ## Japan,_Torishima                0.37                0.58               0.51
    ## NWHI,_FFS                       0.33                0.47               0.81
    ## NWHI,_Kure                      0.54                0.62               0.62
    ## NWHI,_Laysan                    0.33                0.67               0.56
    ## NWHI,_Midway                    0.21                0.29               0.50
    ##                  scaffold_47_2495309 scaffold_487_133339 scaffold_48_3128896
    ## Japan,_Torishima                0.47                0.70                0.47
    ## NWHI,_FFS                       0.25                0.47                0.64
    ## NWHI,_Kure                      0.17                0.38                0.54
    ## NWHI,_Laysan                    0.11                0.44                0.56
    ## NWHI,_Midway                    0.50                0.36                0.64
    ##                  scaffold_491_362387 scaffold_491_382261 scaffold_4_4005861
    ## Japan,_Torishima                0.40                0.00               0.44
    ## NWHI,_FFS                       0.58                0.19               0.22
    ## NWHI,_Kure                      0.46                0.08               0.17
    ## NWHI,_Laysan                    0.44                0.00               0.22
    ## NWHI,_Midway                    0.50                0.29               0.21
    ##                  scaffold_500_214884 scaffold_519_190199 scaffold_51_1220763
    ## Japan,_Torishima                0.16                0.44                0.28
    ## NWHI,_FFS                       0.39                0.31                0.11
    ## NWHI,_Kure                      0.29                0.54                0.17
    ## NWHI,_Laysan                    0.78                0.56                0.78
    ## NWHI,_Midway                    0.64                0.36                0.29
    ##                  scaffold_526_334093 scaffold_52_723164 scaffold_532_590089
    ## Japan,_Torishima                0.21               0.37                0.44
    ## NWHI,_FFS                       0.50               0.42                0.47
    ## NWHI,_Kure                      0.38               0.25                0.50
    ## NWHI,_Laysan                    0.56               0.44                0.33
    ## NWHI,_Midway                    0.29               0.36                0.36
    ##                  scaffold_547_47246 scaffold_54_1630779 scaffold_552_154281
    ## Japan,_Torishima               0.47                0.16                0.42
    ## NWHI,_FFS                      0.47                0.22                0.50
    ## NWHI,_Kure                     0.38                0.42                0.33
    ## NWHI,_Laysan                   0.33                0.33                0.11
    ## NWHI,_Midway                   0.43                0.36                0.43
    ##                  scaffold_555_302525 scaffold_557_489027 scaffold_55_2977720
    ## Japan,_Torishima                0.40                0.33                0.42
    ## NWHI,_FFS                       0.39                0.28                0.28
    ## NWHI,_Kure                      0.42                0.38                0.58
    ## NWHI,_Laysan                    0.22                0.00                0.44
    ## NWHI,_Midway                    0.14                0.14                0.21
    ##                  scaffold_565_253439 scaffold_56_1290372 scaffold_572_19499
    ## Japan,_Torishima                0.47                0.23                  0
    ## NWHI,_FFS                       0.44                0.67                  0
    ## NWHI,_Kure                      0.58                0.33                  0
    ## NWHI,_Laysan                    0.44                0.78                  0
    ## NWHI,_Midway                    0.36                0.29                  0
    ##                  scaffold_57_1788671 scaffold_582_107987 scaffold_582_222480
    ## Japan,_Torishima                0.26                0.09                0.07
    ## NWHI,_FFS                       0.25                0.39                0.33
    ## NWHI,_Kure                      0.38                0.46                0.29
    ## NWHI,_Laysan                    0.44                0.11                0.78
    ## NWHI,_Midway                    0.29                0.43                0.71
    ##                  scaffold_5_697809 scaffold_605_114686 scaffold_608_211809
    ## Japan,_Torishima              0.07                   0                   1
    ## NWHI,_FFS                     0.36                   0                   1
    ## NWHI,_Kure                    0.29                   0                   1
    ## NWHI,_Laysan                  0.67                   0                   1
    ## NWHI,_Midway                  0.14                   0                   1
    ##                  scaffold_60_341016 scaffold_612_363793 scaffold_621_290581
    ## Japan,_Torishima               0.42                0.30                0.33
    ## NWHI,_FFS                      0.53                0.44                0.42
    ## NWHI,_Kure                     0.42                0.21                0.33
    ## NWHI,_Laysan                   0.33                0.11                0.44
    ## NWHI,_Midway                   0.64                0.21                0.36
    ##                  scaffold_62_2806526 scaffold_633_500454 scaffold_64_2598599
    ## Japan,_Torishima                0.40                0.72                0.30
    ## NWHI,_FFS                       0.61                0.64                0.36
    ## NWHI,_Kure                      0.25                0.58                0.50
    ## NWHI,_Laysan                    0.33                0.33                0.33
    ## NWHI,_Midway                    0.50                0.71                0.64
    ##                  scaffold_65_986791 scaffold_670_51777 scaffold_67_2416699
    ## Japan,_Torishima               0.16               0.07                0.49
    ## NWHI,_FFS                      0.58               0.47                0.50
    ## NWHI,_Kure                     0.29               0.29                0.54
    ## NWHI,_Laysan                   0.56               0.67                0.33
    ## NWHI,_Midway                   0.43               0.50                0.64
    ##                  scaffold_684_229342 scaffold_691_412074 scaffold_694_285663
    ## Japan,_Torishima                0.12                0.56                0.33
    ## NWHI,_FFS                       0.25                0.50                0.31
    ## NWHI,_Kure                      0.42                0.46                0.33
    ## NWHI,_Laysan                    0.11                0.44                0.78
    ## NWHI,_Midway                    0.50                0.14                0.43
    ##                  scaffold_698_186739 scaffold_700_166185 scaffold_71_2209803
    ## Japan,_Torishima                0.21                   1                0.14
    ## NWHI,_FFS                       0.25                   1                0.19
    ## NWHI,_Kure                      0.25                   1                0.25
    ## NWHI,_Laysan                    0.22                   1                0.00
    ## NWHI,_Midway                    0.57                   1                0.36
    ##                  scaffold_72_1463271 scaffold_738_79656 scaffold_756_67966
    ## Japan,_Torishima                0.23                  0               0.37
    ## NWHI,_FFS                       0.17                  0               0.61
    ## NWHI,_Kure                      0.25                  0               0.50
    ## NWHI,_Laysan                    0.11                  0               0.11
    ## NWHI,_Midway                    0.14                  0               0.36
    ##                  scaffold_75_1184494 scaffold_760_180618 scaffold_768_353388
    ## Japan,_Torishima                   1                0.12                0.30
    ## NWHI,_FFS                          1                0.11                0.44
    ## NWHI,_Kure                         1                0.46                0.33
    ## NWHI,_Laysan                       1                0.11                0.11
    ## NWHI,_Midway                       1                0.14                0.43
    ##                  scaffold_76_2683423 scaffold_774_111773 scaffold_786_264628
    ## Japan,_Torishima                0.33                0.42                0.12
    ## NWHI,_FFS                       0.42                0.42                0.61
    ## NWHI,_Kure                      0.21                0.33                0.42
    ## NWHI,_Laysan                    0.67                0.44                0.33
    ## NWHI,_Midway                    0.29                0.29                0.64
    ##                  scaffold_78_1735985 scaffold_793_76520 scaffold_7_15896
    ## Japan,_Torishima                0.19               0.56             0.77
    ## NWHI,_FFS                       0.19               0.56             0.97
    ## NWHI,_Kure                      0.21               0.29             0.92
    ## NWHI,_Laysan                    0.33               0.67             1.00
    ## NWHI,_Midway                    0.14               0.50             0.64
    ##                  scaffold_7_4109482 scaffold_820_286874 scaffold_822_90287
    ## Japan,_Torishima               0.63                0.47               0.21
    ## NWHI,_FFS                      0.42                0.67               0.25
    ## NWHI,_Kure                     0.38                0.58               0.29
    ## NWHI,_Laysan                   0.11                0.44               0.67
    ## NWHI,_Midway                   0.21                0.43               0.21
    ##                  scaffold_82_879210 scaffold_834_252344 scaffold_84_2655661
    ## Japan,_Torishima               0.63                0.19                0.44
    ## NWHI,_FFS                      0.42                0.17                0.47
    ## NWHI,_Kure                     0.25                0.21                0.42
    ## NWHI,_Laysan                   0.44                0.11                0.67
    ## NWHI,_Midway                   0.43                0.14                0.21
    ##                  scaffold_854_86476 scaffold_857_326525 scaffold_869_275845
    ## Japan,_Torishima               0.00                0.16                0.60
    ## NWHI,_FFS                      0.22                0.53                0.56
    ## NWHI,_Kure                     0.25                0.58                0.46
    ## NWHI,_Laysan                   0.22                0.56                0.22
    ## NWHI,_Midway                   0.21                0.43                0.43
    ##                  scaffold_88_684287 scaffold_8_4227715 scaffold_91_2132606
    ## Japan,_Torishima               0.44               0.37                0.40
    ## NWHI,_FFS                      0.42               0.14                0.19
    ## NWHI,_Kure                     0.50               0.17                0.12
    ## NWHI,_Laysan                   0.33               0.11                0.22
    ## NWHI,_Midway                   0.29               0.14                0.29
    ##                  scaffold_91_555426 scaffold_925_195188 scaffold_94_1302246
    ## Japan,_Torishima                  0                0.53                0.30
    ## NWHI,_FFS                         0                0.75                0.19
    ## NWHI,_Kure                        0                0.92                0.42
    ## NWHI,_Laysan                      0                0.67                0.00
    ## NWHI,_Midway                      0                0.86                0.50
    ##                  scaffold_958_223818 scaffold_95_2059721 scaffold_97_1719991
    ## Japan,_Torishima                0.44                0.14                0.49
    ## NWHI,_FFS                       0.42                0.39                0.53
    ## NWHI,_Kure                      0.46                0.42                0.46
    ## NWHI,_Laysan                    0.56                0.33                0.67
    ## NWHI,_Midway                    0.29                0.36                0.36
    ##                  scaffold_984_180851 scaffold_990_59711 scaffold_9_4398937
    ## Japan,_Torishima                0.26                  0               0.28
    ## NWHI,_FFS                       0.17                  0               0.56
    ## NWHI,_Kure                      0.21                  0               0.42
    ## NWHI,_Laysan                    0.56                  0               0.33
    ## NWHI,_Midway                    0.36                  0               0.36

``` r
par(mar=c(5.5, 4.5, 1, 1))
  Hobs.pop <- apply(Hobs, MARGIN = 1, FUN = mean)
  Hexp.pop <- apply(Hexp, 1, mean) 
  barplot(Hexp.pop, ylim=c(0,1), las=3, ylab="Expected heterozygosity")
```

![](02-locus-HWE-evaluation_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
  barplot(Hobs.pop, ylim=c(0,1), las=3, ylab="Observed heterozygosity")
```

![](02-locus-HWE-evaluation_files/figure-gfm/unnamed-chunk-18-2.png)<!-- -->

``` r
bfal_diversity <- data.frame(Pop = names(Hobs.pop),
                              n_samples = sum_output$n.by.pop,
                              Hobs = Hobs.pop,
                              Hexp = Hexp.pop)
                              #Ar = Richness$mean.richness)
as.tibble(bfal_diversity) %>%
  rename(Population = Pop, H_obs = Hobs, H_exp = Hexp) #%>%
```

    ## Warning: `as.tibble()` was deprecated in tibble 2.0.0.

    ## Warning: Please use `as_tibble()` instead.

    ## Warning: The signature and semantics have changed, see `?as_tibble`.

    ## # A tibble: 5 x 4
    ##   Population       n_samples H_obs H_exp
    ##   <chr>                <dbl> <dbl> <dbl>
    ## 1 Japan,_Torishima        43 0.315 0.316
    ## 2 NWHI,_FFS               36 0.366 0.363
    ## 3 NWHI,_Kure              24 0.358 0.353
    ## 4 NWHI,_Laysan             9 0.349 0.373
    ## 5 NWHI,_Midway            14 0.358 0.347

``` r
  #write_csv("csv_outputs/reference_pop_hz_summary.csv")
```

chi^2: value of the classical chi-squared test statistic df: degrees of
freedom of the chi-squared test Pr(chi^2 \>): p-value of the chi-squared
test (‘\>’ indicates that the alternative is ‘greater,’ which is always
the case for a chi-squared test) Pr.exact: p-value from an exact test
based on Monte Carlo permutation of alleles (for diploids only). The
default is B = 1000 permutations (set B = 0 to skip this test). Here we
use the function ‘round’ with argument ‘digits = 3’ to round all values
to 3 decimals.

<https://rdrr.io/cran/pegas/man/hw.test.html>

In this case, the degrees of freedom depend on the number alleles for a
given locus.

``` r
# HWE test for all loci - but really what we want is all loci in each population (farther below)
hwe_output <- round(pegas::hw.test(genind_df, B = 1000), digits = 3)
```

    ## Registered S3 method overwritten by 'pegas':
    ##   method      from
    ##   print.amova ade4

    ## Warning in hw.test.loci(x = x, B = B, ...): The following loci were ignored: scaffold_1164_119437
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_1226_115783
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_143_1756997
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_148_727111
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_14_2621198
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_166_410622
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_184_734991
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_204_1432955
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_20_1133858
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_227_1219626
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_229_770334
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_238_668548
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_28_4260845
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_336_143440
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_33_2811656
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_360_269450
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_369_134745
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_411_530837
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_472_10182
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_4_4005861
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_555_302525
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_64_2598599
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_72_1463271
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_738_79656
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_7_15896
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_84_2655661
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_91_2132606
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_94_1302246
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_95_2059721
    ## (not the same ploidy for all individuals, or too many missing data)

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

``` r
# which ones are statistically sign.
rownames_to_column(data.frame(hwe_output)) %>%
  filter(Pr.exact < 0.01)
```

    ##                 rowname   chi.2 df Pr.chi.2... Pr.exact
    ## 1    scaffold_0_4543438  26.921  6       0.000    0.000
    ## 2  scaffold_116_1932131   9.560  1       0.002    0.005
    ## 3  scaffold_123_1115977  14.317  3       0.003    0.000
    ## 4  scaffold_127_1105814  19.087  6       0.004    0.005
    ## 5  scaffold_1327_103421 143.357  3       0.000    0.000
    ## 6  scaffold_164_1134575  39.078  1       0.000    0.000
    ## 7  scaffold_185_1133045  11.926  1       0.001    0.001
    ## 8  scaffold_185_1154507   7.980  1       0.005    0.001
    ## 9   scaffold_306_944429   9.828  1       0.002    0.004
    ## 10  scaffold_345_694649  13.944  1       0.000    0.000
    ## 11  scaffold_429_778790   7.503  1       0.006    0.007
    ## 12  scaffold_51_1220763 253.446  6       0.000    0.004
    ## 13   scaffold_572_19499 126.000  1       0.000    0.000
    ## 14  scaffold_582_107987  22.894  1       0.000    0.000
    ## 15  scaffold_582_222480  17.441  1       0.000    0.000
    ## 16  scaffold_608_211809 126.000  1       0.000    0.000
    ## 17  scaffold_612_363793 126.019  3       0.000    0.008
    ## 18  scaffold_621_290581 126.010  3       0.000    0.006
    ## 19   scaffold_65_986791   8.564  1       0.003    0.005
    ## 20  scaffold_684_229342 126.006  3       0.000    0.005
    ## 21  scaffold_698_186739  23.069  1       0.000    0.000
    ## 22  scaffold_700_166185 126.000  3       0.000    0.000
    ## 23  scaffold_75_1184494 126.000 10       0.000    0.000
    ## 24  scaffold_760_180618  28.854  1       0.000    0.000
    ## 25  scaffold_78_1735985  45.655  1       0.000    0.000
    ## 26   scaffold_7_4109482 126.939  3       0.000    0.004
    ## 27  scaffold_925_195188  19.592  6       0.003    0.000

28 loci that are out of HWE globally.

After looking at the loci x population HWE info, I can remove any that
are out of HWE in the majority of populations.

``` r
# Chi-squared test: p-value
HWE.test <- data.frame(sapply(seppop(genind_df), 
                              function(ls) pegas::hw.test(ls, B=0)[,3]))
```

    ## Warning in hw.test.loci(x = x, B = B, ...): The following loci were ignored: scaffold_1164_119437
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_1226_115783
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_143_1756997
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_166_410622
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_20_1133858
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_336_143440
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_33_2811656
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_360_269450
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_369_134745
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_411_530837
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_555_302525
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_738_79656
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_84_2655661
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_94_1302246
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_95_2059721
    ## (not the same ploidy for all individuals, or too many missing data)

    ## Warning in hw.test.loci(x = x, B = B, ...): The following loci were ignored: scaffold_1164_119437
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_143_1756997
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_166_410622
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_184_734991
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_20_1133858
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_28_4260845
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_360_269450
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_738_79656
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_84_2655661
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_95_2059721
    ## (not the same ploidy for all individuals, or too many missing data)

    ## Warning in hw.test.loci(x = x, B = B, ...): The following loci were ignored: scaffold_1164_119437
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_143_1756997
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_14_2621198
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_166_410622
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_20_1133858
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_238_668548
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_336_143440
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_360_269450
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_555_302525
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_738_79656
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_84_2655661
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_95_2059721
    ## (not the same ploidy for all individuals, or too many missing data)

    ## Warning in hw.test.loci(x = x, B = B, ...): The following loci were ignored: scaffold_143_1756997
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_148_727111
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_20_1133858
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_227_1219626
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_360_269450
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_472_10182
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_738_79656
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_84_2655661
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_95_2059721
    ## (not the same ploidy for all individuals, or too many missing data)

    ## Warning in hw.test.loci(x = x, B = B, ...): The following loci were ignored: scaffold_1164_119437
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_1226_115783
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_143_1756997
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_166_410622
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_204_1432955
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_20_1133858
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_229_770334
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_238_668548
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_336_143440
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_360_269450
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_369_134745
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_4_4005861
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_555_302525
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_64_2598599
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_72_1463271
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_738_79656
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_7_15896
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_84_2655661
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_91_2132606
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_95_2059721
    ## (not the same ploidy for all individuals, or too many missing data)

``` r
HWE.test.chisq <- t(data.matrix(HWE.test))
{cat("Chi-squared test (p-values):", "\n")
round(HWE.test.chisq,3)}
```

    ## Chi-squared test (p-values):

    ##                  scaffold_0_4543438 scaffold_102_694013 scaffold_1050_36334
    ## Japan._Torishima              0.026               0.501               0.035
    ## NWHI._FFS                     0.273               0.064               1.000
    ## NWHI._Kure                    0.650               0.145               0.432
    ## NWHI._Laysan                  0.853               0.764               0.764
    ## NWHI._Midway                  0.424               0.416               0.533
    ##                  scaffold_1053_177405 scaffold_1078_116737 scaffold_10_1333563
    ## Japan._Torishima                0.781                0.939               0.879
    ## NWHI._FFS                       0.279                0.492               0.096
    ## NWHI._Kure                      0.540                0.689               0.973
    ## NWHI._Laysan                    0.317                0.249               0.764
    ## NWHI._Midway                    0.584                0.584               0.803
    ##                  scaffold_115_2200290 scaffold_1164_119437 scaffold_116_1932131
    ## Japan._Torishima                    1                   NA                0.022
    ## NWHI._FFS                           1                   NA                0.239
    ## NWHI._Kure                          1                   NA                0.744
    ## NWHI._Laysan                        1                0.038                0.099
    ## NWHI._Midway                        1                   NA                0.985
    ##                  scaffold_116_386790 scaffold_11_3057945 scaffold_121_1829963
    ## Japan._Torishima               0.876                   1                    1
    ## NWHI._FFS                      0.391                   1                    1
    ## NWHI._Kure                     0.744                   1                    1
    ## NWHI._Laysan                   0.860                   1                    1
    ## NWHI._Midway                   0.773                   1                    1
    ##                  scaffold_121_565678 scaffold_1226_115783 scaffold_123_1115977
    ## Japan._Torishima               0.771                   NA                0.089
    ## NWHI._FFS                      1.000                1.000                0.099
    ## NWHI._Kure                     0.202                0.000                0.333
    ## NWHI._Laysan                   0.612                0.029                0.172
    ## NWHI._Midway                   0.119                   NA                0.611
    ##                  scaffold_1273_1937 scaffold_127_1105814 scaffold_127_901882
    ## Japan._Torishima              0.992                0.584                   1
    ## NWHI._FFS                     0.864                0.044                   1
    ## NWHI._Kure                    1.000                0.000                   1
    ## NWHI._Laysan                  1.000                0.317                   1
    ## NWHI._Midway                  1.000                0.727                   1
    ##                  scaffold_12_5002902 scaffold_1327_103421 scaffold_140_259881
    ## Japan._Torishima               0.708                0.063               0.243
    ## NWHI._FFS                      0.505                0.050               0.960
    ## NWHI._Kure                     0.432                0.200               0.959
    ## NWHI._Laysan                   0.317                0.019               0.099
    ## NWHI._Midway                   0.571                0.033               0.308
    ##                  scaffold_143_1756997 scaffold_146_16700 scaffold_148_727111
    ## Japan._Torishima                   NA              0.102               0.408
    ## NWHI._FFS                          NA              0.518               0.824
    ## NWHI._Kure                         NA              0.532               0.973
    ## NWHI._Laysan                       NA              0.391                  NA
    ## NWHI._Midway                       NA              0.687               0.308
    ##                  scaffold_14_2621198 scaffold_14_5075430 scaffold_155_1707144
    ## Japan._Torishima               0.443               0.953                0.953
    ## NWHI._FFS                      0.007               0.453                0.146
    ## NWHI._Kure                        NA               0.367                0.569
    ## NWHI._Laysan                   0.764               0.764                0.708
    ## NWHI._Midway                   0.212               0.416                0.262
    ##                  scaffold_155_468225 scaffold_157_815403 scaffold_1589_18141
    ## Japan._Torishima               0.153               0.606               0.705
    ## NWHI._FFS                      0.120               0.453               0.516
    ## NWHI._Kure                     0.236               0.763               0.501
    ## NWHI._Laysan                   0.284               0.370               0.878
    ## NWHI._Midway                   0.498               0.109               0.803
    ##                  scaffold_15_2923659 scaffold_15_5063205 scaffold_164_1134575
    ## Japan._Torishima               0.574               0.051                0.000
    ## NWHI._FFS                      0.518               0.017                0.005
    ## NWHI._Kure                     0.026               0.759                0.121
    ## NWHI._Laysan                   0.022               0.249                0.134
    ## NWHI._Midway                   0.435               0.857                0.985
    ##                  scaffold_166_410622 scaffold_16_1050955 scaffold_16_1254862
    ## Japan._Torishima                  NA               0.011               0.134
    ## NWHI._FFS                         NA               0.391               0.449
    ## NWHI._Kure                        NA               0.403               0.258
    ## NWHI._Laysan                   0.549               1.000               1.000
    ## NWHI._Midway                      NA               0.533               0.262
    ##                  scaffold_16_31673 scaffold_16_32172 scaffold_177_673090
    ## Japan._Torishima             0.453             0.453               0.686
    ## NWHI._FFS                    0.285             0.198               1.000
    ## NWHI._Kure                   0.586             0.586               0.959
    ## NWHI._Laysan                 0.764             0.613               0.284
    ## NWHI._Midway                 0.593             0.985               0.308
    ##                  scaffold_184_724429 scaffold_184_734991 scaffold_185_1133045
    ## Japan._Torishima               0.478               0.426                0.166
    ## NWHI._FFS                      0.453                  NA                0.505
    ## NWHI._Kure                     0.014               0.921                0.838
    ## NWHI._Laysan                   0.764               1.000                0.612
    ## NWHI._Midway                   0.193               0.994                0.212
    ##                  scaffold_185_1154507 scaffold_190_1605668 scaffold_199_875998
    ## Japan._Torishima                0.166                0.623               0.257
    ## NWHI._FFS                       0.837                0.453               0.362
    ## NWHI._Kure                      0.736                0.327               0.290
    ## NWHI._Laysan                    0.612                0.317               0.134
    ## NWHI._Midway                    0.193                0.773               0.533
    ##                  scaffold_1_5556176 scaffold_204_1432955 scaffold_204_239685
    ## Japan._Torishima              0.051                0.501               0.388
    ## NWHI._FFS                     0.864                0.120               0.391
    ## NWHI._Kure                    0.586                0.744               0.327
    ## NWHI._Laysan                  0.739                0.099               0.708
    ## NWHI._Midway                  0.262                   NA               0.584
    ##                  scaffold_209_721065 scaffold_20_1133858 scaffold_210_1478805
    ## Japan._Torishima               0.288                  NA                0.749
    ## NWHI._FFS                      0.196                  NA                0.453
    ## NWHI._Kure                     0.959                  NA                0.624
    ## NWHI._Laysan                   0.708                  NA                0.072
    ## NWHI._Midway                   0.571                  NA                0.653
    ##                  scaffold_214_606303 scaffold_223_94277 scaffold_224_319624
    ## Japan._Torishima               0.813                  1               0.534
    ## NWHI._FFS                      0.320                  1               0.824
    ## NWHI._Kure                     0.955                  1               0.569
    ## NWHI._Laysan                   0.317                  1               0.764
    ## NWHI._Midway                   0.791                  1               0.773
    ##                  scaffold_227_1219626 scaffold_229_770334 scaffold_234_1146621
    ## Japan._Torishima                0.443               0.744                1.000
    ## NWHI._FFS                       0.102               0.752                0.585
    ## NWHI._Kure                      0.831               0.103                0.744
    ## NWHI._Laysan                       NA               0.878                0.860
    ## NWHI._Midway                    1.000                  NA                0.533
    ##                  scaffold_237_341989 scaffold_238_668548 scaffold_245_674441
    ## Japan._Torishima               0.534               0.278                   1
    ## NWHI._FFS                      0.658               0.020                   1
    ## NWHI._Kure                     0.569                  NA                   1
    ## NWHI._Laysan                   0.708               0.549                   1
    ## NWHI._Midway                   0.773                  NA                   1
    ##                  scaffold_246_1244055 scaffold_246_31712 scaffold_247_950885
    ## Japan._Torishima                0.456              0.899               0.278
    ## NWHI._FFS                       0.701              0.611               0.324
    ## NWHI._Kure                      0.102              0.734               0.202
    ## NWHI._Laysan                    0.549              0.318               0.391
    ## NWHI._Midway                    0.076              0.823               0.985
    ##                  scaffold_249_340732 scaffold_24_1490422 scaffold_253_684619
    ## Japan._Torishima               0.749               0.758               0.836
    ## NWHI._FFS                      0.022               0.009               0.697
    ## NWHI._Kure                     0.744               0.403               0.327
    ## NWHI._Laysan                   0.860               0.708               0.134
    ## NWHI._Midway                   0.994               0.533               0.360
    ##                  scaffold_256_454799 scaffold_275_1006390 scaffold_27_1646617
    ## Japan._Torishima               0.282                0.577               0.028
    ## NWHI._FFS                      0.965                0.298               0.005
    ## NWHI._Kure                     0.624                0.532               0.242
    ## NWHI._Laysan                   0.302                0.970               0.515
    ## NWHI._Midway                   0.274                0.956               0.384
    ##                  scaffold_284_808709 scaffold_286_314167 scaffold_287_423519
    ## Japan._Torishima               0.443               0.336               0.153
    ## NWHI._FFS                      0.016               0.735               0.755
    ## NWHI._Kure                     0.744               0.838               0.624
    ## NWHI._Laysan                   0.072               0.391               0.764
    ## NWHI._Midway                   0.533               0.347               0.584
    ##                  scaffold_28_4260845 scaffold_28_4262705 scaffold_298_359540
    ## Japan._Torishima               0.453               0.738               0.336
    ## NWHI._FFS                         NA               0.198               0.333
    ## NWHI._Kure                     0.973               0.622               0.744
    ## NWHI._Laysan                   0.391               0.391               0.708
    ## NWHI._Midway                   0.212               0.852               0.019
    ##                  scaffold_298_460712 scaffold_29_3347800 scaffold_306_328368
    ## Japan._Torishima               0.876               0.517               0.511
    ## NWHI._FFS                      0.453               0.324               0.711
    ## NWHI._Kure                     0.656               0.102               0.578
    ## NWHI._Laysan                   1.000               0.134               0.284
    ## NWHI._Midway                   0.773               0.057               0.611
    ##                  scaffold_306_944429 scaffold_310_880276 scaffold_316_437138
    ## Japan._Torishima               0.014               0.534               0.560
    ## NWHI._FFS                      0.168               0.230               0.887
    ## NWHI._Kure                     0.042               1.000               0.987
    ## NWHI._Laysan                   0.284               0.860               0.795
    ## NWHI._Midway                   0.134               0.890               0.229
    ##                  scaffold_31_750839 scaffold_324_135439 scaffold_32_3673235
    ## Japan._Torishima              0.976               0.093               0.930
    ## NWHI._FFS                     0.391               0.825               0.114
    ## NWHI._Kure                    0.102               0.431               0.763
    ## NWHI._Laysan                  0.370               0.003               0.612
    ## NWHI._Midway                  0.852               0.308               0.533
    ##                  scaffold_32_803437 scaffold_32_811415 scaffold_335_662765
    ## Japan._Torishima              0.876              0.876               0.388
    ## NWHI._FFS                     0.196              0.352               0.391
    ## NWHI._Kure                    0.197              0.197               0.243
    ## NWHI._Laysan                  0.708              0.948               1.000
    ## NWHI._Midway                  0.308              0.653               0.890
    ##                  scaffold_336_143440 scaffold_33_2797456 scaffold_33_2811656
    ## Japan._Torishima                  NA               0.282                  NA
    ## NWHI._FFS                      0.658               0.658               0.094
    ## NWHI._Kure                        NA               0.102               0.129
    ## NWHI._Laysan                   0.613               0.549               0.261
    ## NWHI._Midway                      NA               0.584               0.594
    ##                  scaffold_342_417561 scaffold_343_778064 scaffold_343_863766
    ## Japan._Torishima               0.063               0.758                1.00
    ## NWHI._FFS                      0.230               0.569                1.00
    ## NWHI._Kure                     0.744               0.811                1.00
    ## NWHI._Laysan                   0.072               0.739                0.86
    ## NWHI._Midway                   0.313               0.031                1.00
    ##                  scaffold_345_694649 scaffold_34_2385714 scaffold_351_703875
    ## Japan._Torishima               0.389               0.606               0.184
    ## NWHI._FFS                      0.960               0.310               0.450
    ## NWHI._Kure                     0.967               0.959               0.676
    ## NWHI._Laysan                   0.099               0.764               0.764
    ## NWHI._Midway                   0.360               0.533               0.852
    ##                  scaffold_356_86112 scaffold_360_269450 scaffold_369_134745
    ## Japan._Torishima              0.623                  NA                  NA
    ## NWHI._FFS                     0.585                  NA               0.837
    ## NWHI._Kure                    0.403                  NA               0.106
    ## NWHI._Laysan                  0.099                  NA               0.764
    ## NWHI._Midway                  1.000                  NA                  NA
    ##                  scaffold_381_698673 scaffold_397_13764 scaffold_40_612943
    ## Japan._Torishima               0.166              0.811              0.050
    ## NWHI._FFS                      0.505              0.626              0.242
    ## NWHI._Kure                     0.053              0.548              0.394
    ## NWHI._Laysan                   1.000              0.478              0.878
    ## NWHI._Midway                   0.653              0.274              0.883
    ##                  scaffold_411_530837 scaffold_417_428513 scaffold_41_907611
    ## Japan._Torishima                  NA               1.000              0.464
    ## NWHI._FFS                      0.367               0.845              0.835
    ## NWHI._Kure                     0.056               0.656              0.073
    ## NWHI._Laysan                   0.134               0.391              1.000
    ## NWHI._Midway                   0.347               0.533              0.803
    ##                  scaffold_429_778790 scaffold_437_192045 scaffold_45_3079339
    ## Japan._Torishima               0.389               0.813               0.590
    ## NWHI._FFS                      0.697               0.830               0.013
    ## NWHI._Kure                     0.431               0.676               0.186
    ## NWHI._Laysan                   0.391               0.708               0.764
    ## NWHI._Midway                   0.571               0.859               0.308
    ##                  scaffold_460_237213 scaffold_461_674999 scaffold_472_10182
    ## Japan._Torishima               0.182               0.286              0.632
    ## NWHI._FFS                      0.046               0.606              0.109
    ## NWHI._Kure                     0.367               0.186              0.065
    ## NWHI._Laysan                   0.370               0.294                 NA
    ## NWHI._Midway                   0.653               0.262              0.148
    ##                  scaffold_47_2495309 scaffold_487_133339 scaffold_48_3128896
    ## Japan._Torishima               0.501               0.007               0.175
    ## NWHI._FFS                      0.391               0.939               0.273
    ## NWHI._Kure                     0.656               0.400               0.194
    ## NWHI._Laysan                   0.860               0.731               0.562
    ## NWHI._Midway                   0.212               0.882               0.604
    ##                  scaffold_491_362387 scaffold_491_382261 scaffold_4_4005861
    ## Japan._Torishima               0.590               1.000              0.438
    ## NWHI._FFS                      0.981               0.518              0.384
    ## NWHI._Kure                     0.685               0.831              0.978
    ## NWHI._Laysan                   0.764               1.000              0.708
    ## NWHI._Midway                   0.857               0.571                 NA
    ##                  scaffold_500_214884 scaffold_519_190199 scaffold_51_1220763
    ## Japan._Torishima               0.561               0.868               0.606
    ## NWHI._FFS                      0.187               0.211               0.724
    ## NWHI._Kure                     0.056               0.625               0.050
    ## NWHI._Laysan                   0.056               0.613               0.002
    ## NWHI._Midway                   0.274               0.291               0.533
    ##                  scaffold_526_334093 scaffold_52_723164 scaffold_532_590089
    ## Japan._Torishima               0.443              0.882               0.785
    ## NWHI._FFS                      0.046              0.505               0.912
    ## NWHI._Kure                     0.377              0.484               0.889
    ## NWHI._Laysan                   0.249              0.391               0.612
    ## NWHI._Midway                   0.571              0.859               0.347
    ##                  scaffold_547_47246 scaffold_54_1630779 scaffold_552_154281
    ## Japan._Torishima              0.698               0.953               0.960
    ## NWHI._FFS                     0.965               0.670               0.755
    ## NWHI._Kure                    0.804               0.586               0.327
    ## NWHI._Laysan                  0.370               0.023               0.860
    ## NWHI._Midway                  0.852               0.416               0.308
    ##                  scaffold_555_302525 scaffold_557_489027 scaffold_55_2977720
    ## Japan._Torishima                  NA               0.623               0.606
    ## NWHI._FFS                      0.548               0.497               0.024
    ## NWHI._Kure                        NA               0.804               0.044
    ## NWHI._Laysan                   0.284               0.003               0.764
    ## NWHI._Midway                      NA               0.773               0.057
    ##                  scaffold_565_253439 scaffold_56_1290372 scaffold_572_19499
    ## Japan._Torishima               0.324               0.388                  1
    ## NWHI._FFS                      0.681               0.046                  0
    ## NWHI._Kure                     0.536               0.157                  1
    ## NWHI._Laysan                   1.000               0.056                  1
    ## NWHI._Midway                   0.698               0.109                  1
    ##                  scaffold_57_1788671 scaffold_582_107987 scaffold_582_222480
    ## Japan._Torishima               0.983               0.749               0.017
    ## NWHI._FFS                      0.352               0.148               0.134
    ## NWHI._Kure                     0.533               0.432               0.116
    ## NWHI._Laysan                   1.000               0.072               0.056
    ## NWHI._Midway                   0.571               0.308               0.109
    ##                  scaffold_5_697809 scaffold_605_114686 scaffold_608_211809
    ## Japan._Torishima             0.813                   1               0.000
    ## NWHI._FFS                    0.672                   1               0.000
    ## NWHI._Kure                   0.393                   1               0.000
    ## NWHI._Laysan                 0.134                   1               0.003
    ## NWHI._Midway                 0.773                   1               0.000
    ##                  scaffold_60_341016 scaffold_612_363793 scaffold_621_290581
    ## Japan._Torishima              0.400               0.744               0.564
    ## NWHI._FFS                     0.705               0.518               0.352
    ## NWHI._Kure                    0.033               0.422               0.327
    ## NWHI._Laysan                  0.612               0.029               0.029
    ## NWHI._Midway                  0.766               0.653               0.859
    ##                  scaffold_62_2806526 scaffold_633_500454 scaffold_64_2598599
    ## Japan._Torishima               0.048               0.348               0.080
    ## NWHI._FFS                      0.440               0.962               0.570
    ## NWHI._Kure                     0.229               0.745               0.303
    ## NWHI._Laysan                   0.382               0.317               0.317
    ## NWHI._Midway                   0.742               0.772                  NA
    ##                  scaffold_65_986791 scaffold_670_51777 scaffold_67_2416699
    ## Japan._Torishima              0.561              0.813               0.106
    ## NWHI._FFS                     0.085              0.606               0.864
    ## NWHI._Kure                    0.116              0.393               0.069
    ## NWHI._Laysan                  0.739              0.294               0.612
    ## NWHI._Midway                  0.640              0.212               0.193
    ##                  scaffold_684_229342 scaffold_691_412074 scaffold_694_285663
    ## Japan._Torishima               0.145               0.272               0.915
    ## NWHI._FFS                      0.352               0.140               0.279
    ## NWHI._Kure                     0.197               0.736               0.327
    ## NWHI._Laysan                   0.029               1.000               0.096
    ## NWHI._Midway                   0.212               0.773               0.308
    ##                  scaffold_698_186739 scaffold_700_166185 scaffold_71_2209803
    ## Japan._Torishima               0.000               0.000               0.623
    ## NWHI._FFS                      0.011               0.000               0.505
    ## NWHI._Kure                     0.102               0.000               0.484
    ## NWHI._Laysan                   0.134               0.029               0.003
    ## NWHI._Midway                   0.593               0.000               0.416
    ##                  scaffold_72_1463271 scaffold_738_79656 scaffold_756_67966
    ## Japan._Torishima               0.836                 NA              0.317
    ## NWHI._FFS                      0.349                 NA              0.086
    ## NWHI._Kure                     0.236                 NA              0.540
    ## NWHI._Laysan                   0.072                 NA              0.022
    ## NWHI._Midway                      NA                 NA              0.859
    ##                  scaffold_75_1184494 scaffold_760_180618 scaffold_768_353388
    ## Japan._Torishima               0.000               0.000               0.243
    ## NWHI._FFS                      0.000               0.009               0.777
    ## NWHI._Kure                     0.000               0.689               0.221
    ## NWHI._Laysan                   0.174               0.860               0.860
    ## NWHI._Midway                   0.003               0.010               0.852
    ##                  scaffold_76_2683423 scaffold_774_111773 scaffold_786_264628
    ## Japan._Torishima               0.564               0.321               0.686
    ## NWHI._FFS                      0.802               0.367               0.086
    ## NWHI._Kure                     0.422               0.102               0.586
    ## NWHI._Laysan                   0.321               0.764               0.370
    ## NWHI._Midway                   0.571               0.109               0.076
    ##                  scaffold_78_1735985 scaffold_793_76520 scaffold_7_15896
    ## Japan._Torishima               0.000              0.245            0.011
    ## NWHI._FFS                      0.000              0.969            0.000
    ## NWHI._Kure                     0.004              0.615            0.007
    ## NWHI._Laysan                   0.370              0.558            0.212
    ## NWHI._Midway                   0.015              0.956               NA
    ##                  scaffold_7_4109482 scaffold_820_286874 scaffold_822_90287
    ## Japan._Torishima              0.080               0.698              0.443
    ## NWHI._FFS                     0.802               0.046              0.837
    ## NWHI._Kure                    0.533               0.414              0.403
    ## NWHI._Laysan                  0.004               0.391              0.294
    ## NWHI._Midway                  0.039               0.803              0.653
    ##                  scaffold_82_879210 scaffold_834_252344 scaffold_84_2655661
    ## Japan._Torishima              0.093               0.501                  NA
    ## NWHI._FFS                     0.352               0.585                  NA
    ## NWHI._Kure                    0.484               0.422                  NA
    ## NWHI._Laysan                  0.391               0.860                  NA
    ## NWHI._Midway                  0.803               0.773                  NA
    ##                  scaffold_854_86476 scaffold_857_326525 scaffold_869_275845
    ## Japan._Torishima              1.000               0.561               0.112
    ## NWHI._FFS                     0.670               0.647               0.492
    ## NWHI._Kure                    0.484               0.414               0.145
    ## NWHI._Laysan                  0.708               0.613               0.284
    ## NWHI._Midway                  0.653               0.640               0.308
    ##                  scaffold_88_684287 scaffold_8_4227715 scaffold_91_2132606
    ## Japan._Torishima              0.639              0.134               0.330
    ## NWHI._FFS                     0.352              0.654               0.518
    ## NWHI._Kure                    0.744              0.243               0.744
    ## NWHI._Laysan                  0.612              0.860               0.099
    ## NWHI._Midway                  0.571              0.773                  NA
    ##                  scaffold_91_555426 scaffold_925_195188 scaffold_94_1302246
    ## Japan._Torishima                  1               0.268                  NA
    ## NWHI._FFS                         1               0.178               0.505
    ## NWHI._Kure                        1               0.025               0.197
    ## NWHI._Laysan                      1               0.558               0.003
    ## NWHI._Midway                      1               0.104               0.985
    ##                  scaffold_958_223818 scaffold_95_2059721 scaffold_97_1719991
    ## Japan._Torishima               0.639                  NA               0.605
    ## NWHI._FFS                      0.319                  NA               0.593
    ## NWHI._Kure                     0.993                  NA               0.672
    ## NWHI._Laysan                   0.613                  NA               0.946
    ## NWHI._Midway                   0.571                  NA               0.180
    ##                  scaffold_984_180851 scaffold_990_59711 scaffold_9_4398937
    ## Japan._Torishima               0.464                  1              0.267
    ## NWHI._FFS                      0.585                  1              0.771
    ## NWHI._Kure                     0.422                  1              0.149
    ## NWHI._Laysan                   0.739                  1              0.261
    ## NWHI._Midway                   0.416                  1              0.398

``` r
# these are the p-values for the per-pop calcs
loci_to_toss_hwe <- rownames_to_column(as.data.frame(HWE.test.chisq), var = "population") %>%
  pivot_longer(cols= 2:length(rownames_to_column(as.data.frame(HWE.test.chisq))), names_to = "locus", values_to = "p_val") %>%
  filter(p_val < 0.05) %>% # which p-values are significant?
  group_by( locus) %>%
  tally() %>% # how many loci are out of HWE in how many populations?
  arrange(desc(n)) %>%
  filter(n>3)
```

Remove loci that are statistically out of HWE in \>3 populations (the
majority).

That’s only 4 loci.

``` r
loci_to_toss_hwe %>%
  write_csv("csv_outputs/loci_to_toss_hwe.csv")

loci_to_keep <- rownames_to_column(as.data.frame(HWE.test.chisq), var = "population") %>%
  pivot_longer(cols= 2:length(rownames_to_column(as.data.frame(HWE.test.chisq))), names_to = "locus", values_to = "p_val") %>%
  select(locus) %>%
  unique() %>%
  anti_join(., loci_to_toss_hwe)
```

    ## Joining, by = "locus"

``` r
loci_to_keep %>%
  write_csv("csv_outputs/loci_to_keep_hwe.csv")
# from here, I can generate a "final" dataset.
```
