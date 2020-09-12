# Usage of scripts for preprocessing to generate lookup tables

## Command to run SISTR
python SISTR.py constant_het denom_het constant_common denom_common eps_bins inFile use_het use_common use_bins num_bins abc_model lrt_model LRT_num_sims 

### Command line arguments explanations:
constant_het: constant in numerator for heterozygosity error tolerance formula  
denom_het: constant in denominator for heterozygosity error tolerance formula  
constant_common: constant in numerator for number of common alleles error tolerance formula  
denom_common: constant in denominator for number of common alleles error tolerance formula  
eps_bins: error tolerance for bins summary statistic  
inFile: file with allele frequencies (see below for format specifications)  
use_het: whether to use heterozygosity as summary statistic (y or n)  
use_common: whether to use number of common alleles as summary statistic (y or n)  
use_bins: whether to use bins as summary statistic (y or n)  
num_bins: number of bins to use  
abc_model: folder name for ABC lookup tables  
lrt_model: folder name for LRT lookup tables  
LRT_num_sims: number of simulations used for LRT in preprocessing step  

### Example command (with parameters used in Mitra, et al. 2020):
python SISTR_v1.py 0.005 3 1 4 0.3 allele_freqs_example.txt y n y 5 eurodem_prior2 eurodem_0810 2000 

### Input file specifications:
The input file should be a tab delimited file with 14 columns: chrom, start, end, allele_freqs, total, period, motif, coding, intron, UTR5, UTR3, promoter5kb, intergenic, tss_gene.  
Each row should be a separate locus. See example input file 'allele_freqs_example.txt'.

### Output file:
The output file will be a tab delimited file with 24 columns: chrom, start, end, total, period, optimal_ru, motif, coding, intron, UTR5, UTR3, promoter5kb, intergenic, tss_gene, het, common, bins, ABC_s_95_CI, Num_s_accepted_(max_10000), Likelihood_0, Likelihood_s, LR, LogLR, LRT_p_value.  
Each row is a separate locus.