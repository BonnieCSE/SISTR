# Selection Inference on Short Tandem Repeats (SISTR)
SISTR is a method to measure negative selection at short tandem repeats.  

It takes allele frequency data (per-locus frequencies of each allele length) as input and outputs a posterior estimate of the selection coefficient at each locus.  

For questions on usage, please contact Bonnie Huang (bbhuang@ucsd.edu)

# Dependencies
SISTR uses Python3 and the following libraries: Python Standard Library, SciPy, NumPy.

# Pipeline
1. The first step is a preprocessing step which involves simulating allele frequencies to generate lookup tables used later in SISTR (see step 2).   

 See `simulations/` for further details. 
 
 Note: Alternatively, the user can skip this step and use lookup tables already generated for community use: (https://drive.google.com/drive/folders/1g70y6z6sU5DVpF6ZGzosNDVFFoXxnqmn?usp=sharing)
 
 Example command to run simulations:  
 ```
 python ABC_lookup.py per opt num_sims a_param b_param filenum out_folder use_TMRCA  
 ```

2. The second step is to run SISTR, which requires (1) allele frequency data and (2) precomputed lookup tables generated in step 1. It outputs a posterior estimate of the selection coefficient at each locus with a 95% confidence interval. 

 See `sistr/` for further details.

 Example command to run SISTR:  
 ```
 python SISTR_v1.py constant_het denom_het constant_common denom_common eps_bins inFile use_het use_common use_bins num_bins abc_model lrt_model LRT_num_sims 
 ```

# Folder Contents
* `helper_functions/`: Contains helper functions used in scripts in `simulations/` and `SISTR/`
* `simulations/`: Contains scripts to generate lookup tables
* `sistr/`: Contains scripts to run SISTR