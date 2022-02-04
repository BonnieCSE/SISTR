# Running SISTR

## Note: Before running SISTR, lookup tables must first be generated. 
**We recommend using the available lookup tables found [here](https://drive.google.com/drive/folders/1nar_UJ_jyS97J_tV2vnwDgNmqLSQCRkm?usp=sharing) by downloading the entire folder `sistr_resources/` and saving it at the same directory where you are running SISTR.** If you wish to create your own custom lookup tables, see `simulations/` for further details. 

## SISTR Usage
Example commands:
```
# Using alternate format for input allele frequencies
python SISTR_v1.py \
     --in-file allele_freqs_test.txt \
     --out-file test_results.txt 
     
# Using motif format for input allele frequencies
python SISTR_v1.py \
     --in-file allele_freqs_test_motif_format.txt \
     --out-file test_results_motif_format.txt \
     --motif-format  
```

Required parameters:  
* __`--in-file <string>`__ name of input file containing allele frequency data (see below for format specifications)  
* __`--out-file <string>`__ name of output file 

Optional parameters (default values are those used in Mitra, et. al. 2020):
* __`--eps-het-numerator <float>`__ constant in numerator for heterozygosity error tolerance formula (default value: 0.005) 
* __`--eps-het-denominator <int>`__ constant in denominator for heterozygosity error tolerance formula (default value: 3)
* __`--eps-bins <float>`__ error tolerance for bins summary statistic (default value: 0.3)
* __`--num-bins <int>`__ number of bins to use (default value: 5)
* __`--lrt-num-sims <int>`__ number of simulations used to generate LRT lookup tables (default value: 2000) 
* __`--abc-lookup-folder <string>`__ path to folder containing ABC lookup tables (default path: `sistr_resources/abc_lookup/`)
* __`--lrt-lookup-folder <string>`__ path to folder containing LRT lookup tables (default path: `sistr_resources/lrt_lookup/`)
* __`--motif-format`__ whether allele frequencies of the input data is in motif format

## Input file format:
The input file should be a tab delimited file with 7 columns: chrom, start, end, allele_freqs, total, period, motif  
Each row should be a separate locus.  
SISTR supports input allele frequencies in two different formats: motif format and alternate format.  
See example input file `allele_freqs_test_motif_format.txt` for the motif format.  
See example input file `allele_freqs_test.txt` for the alternate format.  
**Note: If using motif format for the input allele frequencies, make sure to use the flag `--motif-format` when running SISTR.**  

### Explanation of each column:
* chrom: chromosome number  
* start: starting position of locus in a reference genome  
* end: ending position of locus in a reference genome  
* allele_freqs: allele frequency data (see note below for required format)
* total: sample size (number of alleles) of empirical allele frequencies  
* period: number of base pairs in the repeat motif (e.g. period 2 refers to dinucleotides)  
* motif: repeat motif  

#### Note: SISTR supports two different formats for the allele frequency data: motif format and alternate format. 

**Motif format:**   
Data for each allele present in the population is represented by the STR allele and the number of copies of the allele in the population separated by a colon. Each allele is separated by a comma.  

Example: 
| chrom | start | end | allele_freqs | total | period | motif |
|---|---|---|---|---|---|---|
|1|1000|1025|ATATATATATATATATATATATAT:100,ATATATATATATATATATATATATAT:6000,ATATATATATATATATATATATATATAT:400|6500|2|AT| 

**Alternate format:**  

Data for each allele present in the population is represented by two numbers separated by a colon. The first number is the allele size in base pairs relative to a reference allele size (which is calculated as `end - start + 1`). The second number is the number of copies of the allele in the population. Each allele is separated by a comma.  
   
Example: 
| chrom | start | end | allele_freqs | total | period | motif |
|---|---|---|---|---|---|---|
|1|1000|1025|-2:100,0:6000,4:400|6500|2|AT| 
   
At this locus, the reference allele length is 26 base pairs or 13 repeat units. In the population, there are 3 alleles present: 100 copies of 12 repeat units, 6000 copies of 13 repeat units, and 400 copies of 15 repeat units.  

### How to generate a SISTR input file from a GangSTR VCF file:  

From a GangSTR VCF file (https://github.com/gymreklab/GangSTR), allele frequencies for each STR can be extracted using statSTR from TRTools (https://github.com/gymreklab/TRTools). Please use the `--acount` and `--num-called` flags when running statSTR to get the `allele_freqs` and `total` columns for the SISTR input file. The period and motif for each STR can be found from the reference file used when running GangSTR (https://github.com/gymreklab/GangSTR#gangstr-reference-files).  

## Output file format:
The output file is a tab delimited file with 14 columns: chrom, start, end, total, period, optimal_ru, motif, ABC_s_95_CI, Percent_s_accepted, Likelihood_0, Likelihood_s, LR, LogLR, LRT_p_value    
Each row is a separate locus.  
See example output file `test_results.txt`.

### Explanation of each column:
* optimal_ru: number of repeat units in optimal allele (most frequent allele length in population)
* ABC_s_95_CI: 95% confidence interval of posterior distribution of s obtained from ABC
* Percent_s_accepted: percent of s values accepted from ABC
* Likelihood_0: likelihood s equals 0
* Likelihood_s: likelihood s equals ABC_s_median
* LR: Likelihood_0/Likelihood_s
* LogLR: -2ln(LR)
* LRT_p_value: p value testing whether a model with selection fits better than a model without selection (obtained using a likelihood ratio test)