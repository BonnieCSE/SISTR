# Usage of scripts for preprocessing to generate lookup tables

## Note: Lookup tables used in Mitra et. al. 2020 and available for community use can be found at:
https://drive.google.com/drive/folders/1g70y6z6sU5DVpF6ZGzosNDVFFoXxnqmn?usp=sharing

## Command to generate ABC lookup table (generate one ABC lookup table for each period/optimal allele combination)
python ABC_lookup.py per opt num_sims k_param theta_param filenum out_folder use_TMRCA

### Command line arguments explanations:
per: period (repeat unit length) -> Note: SISTR currently only handles TRs with repeat unit lengths 2-4 bp  
opt: optimal allele -> Note: List of optimal alleles per period to simulate: (per 2, opt 11-20; per 3, opt 5-13; per 4, opt 7-10)  
num_sims: number of simulations  
a_param: shape parameter of prior gamma distribution (use a_param = -1 for a uniform distribution)  
b_param: scale parameter of prior gamma distribution  
filenum: number used in output file name  
out_folder: name of folder to which to output lookup table  
use_TMRCA: 0 for variable number of generations; 1 for constant number of generations (55920)  

### Example command:
python ABC_lookup.py 2 11 5000 0.0881 0.2541 0 prior2 0

### Output file:
The output file will be a tab delimited file with 5 columns: s (s value used), het (heterozygosity), common (number of common alleles, defined as frequency > 0.05), bins (1x5 vector of binned alleles), freqs (allele frequencies), gens (number of generations used).  
Each row is a separate simulation.

## Command to generate LRT lookup table (generate one ABC lookup table for each period/optimal allele combination)
python LRT_lookup.py per opt filenum s_vals num_sims use_TMRCA

### Command line arguments explanations:
per: period (see above)  
opt: optimal allele (see above)  
filenum: number used in output file name  
s_vals: list of s values, separated by commas, for which to generate lookup table  
num_sims: number of simulations per s value  
use_TMRCA: 0 for variable number of generations; 1 for constant number of generations (55920)  

### Example command:
python LRT_lookup.py 2 11 0 0,1e-06,1e-05 2000 0

### Output files:
The output will be a 5 tab delimited files. Each file has 2 columns: s (s value used), information about each simulation separated by colons. The 5 files contain information about heterozygosity, number of common alleles, vector of binned alleles, allele frequencies, and number of generations used.  
Each row is a separate s value. 