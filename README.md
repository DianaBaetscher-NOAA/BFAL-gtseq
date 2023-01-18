# BFAL-gtseq

December 2022


Population assignment of Black-footed Albatross using GTseq and rubias. 

Primers were designed from lcWGS data and will be used for:\
(1) self-assignment of known-colony samples, and \
(2) assignment of bycatch to breeding colonies.

Data are generated on a MiSeq at the AFSC Genetics Program lab in Juneau, AK, \
bioinformatic analyses are performed on the NMFS HPCC, Sedna, based in Seattle, WA, \
and finally, analysis of microhaplotypes, genotype data, and locus-fidelity is performed on Diana's laptop.

gtseq_test1 - 353 loci in a single primer pool, each at 0.25uM \
gtseq_test2 - 351 loci in two separate primer pools, each at 0.25uM \
gtseq_test3 - 284 loci in a single primer pool, with variable concentrations based on read depth results from test2 \


### Approach

For identifying which primer pairs to keep in the pool and which ones to remove, I am looking at two primary characteristics:

1. On-target reads - this is obtained through the GTscore analysis of merged vs. unmerged reads (R1 vs. Flashed reads) \

2. Read-depth - overamplifying primer sets should be removed or the concentration reduced.


