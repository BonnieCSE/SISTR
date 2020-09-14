# Running SISTR

## Note: Before running SISTR, lookup tables must first be generated. 
Custom lookup tables can be created by the user (see `simulations/` for further details) or the user may use the available lookup tables found here (used in Mitra, et. al. 2020): https://drive.google.com/drive/folders/1g70y6z6sU5DVpF6ZGzosNDVFFoXxnqmn?usp=sharing

## SISTR Usage
Example command:
```
python SISTR_v1.py \
     --inFile allele_freqs_test.txt \
     --outFile test_results.txt 
```

Required parameters:  
* `--inFile` name of input file containing allele frequency data (see below for format specifications)  
* `--outFile` name of output file 

Optional parameters (default values are those used in Mitra, et. al. 2020):
* `--eps-het-numerator` constant in numerator for heterozygosity error tolerance formula (default value: 0.005) 
* `--eps-het-denominator` constant in denominator for heterozygosity error tolerance formula (default value: 3)
* `--eps_bins` error tolerance for bins summary statistic (default value: 0.3)
* `--num_bins` number of bins to use (default value: 5)
* `--lrt_num_sims`: number of simulations used to generate LRT lookup tables (default value: 2000) 
* `--abc-lookup-folder` path to folder containing ABC lookup tables (default path: './../lookup_tables/abc/')
* `--lrt-lookup-folder` path to folder containing LRT lookup tables (default path: './../lookup_tables/lrt/')

## Input file format:
The input file should be a tab delimited file with 7 columns: chrom, start, end, allele_freqs, total, period, motif  
Each row should be a separate locus.  
See example input file 'allele_freqs_example.txt'.

### Explanation of each column:
* chrom: chromosome number  
* start: starting position of locus  
* end: ending position of locus  
* allele_freqs: allele frequency data  
* total: sample size (number of alleles) of empirical allele frequencies  
* period: number of base pairs in the repeat motif (e.g. period 2 refers to dinucleotides)  
* motif: repeat motif  

Note: Required format for allele frequency data  
Data for each allele present in the population is represented by two numbers separated by a colon. The first number is the allele size in base pairs relative to the reference allele size (which is calculated as end - start + 1). The second number is the number of copies of the allele in the population. Each allele is separated by a comma.  
   
Example: 
| chrom | start | end | allele_freqs | total | period | motif |
|---|---|---|---|---|---|---|
|1|1000|1025|-2:100,0:6000,4:400|6500|2|AT| 
   
At this locus, the reference allele length is 26 base pairs or 13 repeat units. In the population, there are 3 alleles present: 100 copies of 12 repeat units, 6000 copies of 13 repeat units, and 400 copies of 15 repeat units.

## Output file format:
The output file will be a tab delimited file with 14 columns: chrom, start, end, total, period, optimal_ru, motif, ABC_s_95_CI, Percent_s_accepted, Likelihood_0, Likelihood_s, LR, LogLR, LRT_p_value    
Each row is a separate locus.  
See example output file 'test_results.txt'.

### Explanation of each column:
* optimal_ru: number of repeat units in optimal allele (most frequent allele length in population)
* ABC_s_95_CI: 95% confidence interval of posterior distribution of s obtained from ABC
* Percent_s_accepted: percent of s values accepted from ABC
* Likelihood_0: likelihood s equals 0
* Likelihood_s: likelihood s equals ABC_s_median
* LR: Likelihood_0/Likelihood_s
* LogLR: -2ln(LR)
* LRT_p_value: p value testing whether a model with selection fits better than a model without selection (obtained using a likelihood ratio test)