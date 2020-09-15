# SISTR: Selection Inference on Short Tandem Repeats 
SISTR is a method to measure negative selection at short tandem repeats.  

It takes allele frequency data (per-locus frequencies of each allele length) as input and outputs a posterior estimate of the selection coefficient at each locus.  

For questions on usage, please contact Bonnie Huang (bbhuang@ucsd.edu)  

# Dependencies
SISTR uses Python3 and the following libraries in addition to the Python Standard Library: SciPy, NumPy

# Pipeline
1. The first step is a preprocessing step which involves simulating allele frequencies to generate lookup tables used later in SISTR (see step 2).   

   See `simulations/` for further details. 
   
   Example command to run simulations to generate an ABC lookup table:    
   ```
   python ABC_lookup.py \
     --out-folder ./../sistr_resources_test/abc_lookup/ \
     --per 3 \
     --opt-allele 5 \
     --num-sims 5
   ```
 
   Note: Alternatively, you can skip this step and use lookup tables already generated for community use: https://drive.google.com/drive/folders/1p_QoSQ7gzs7hVEwfJyhGMT-ELXZcWoA_?usp=sharing

2. The second step is to run SISTR, which requires (1) allele frequency data and (2) precomputed lookup tables generated in step 1. It outputs a posterior estimate of the selection coefficient at each locus. 

   See `sistr/` for further details.

   Example command to run SISTR:  
   ```
   python SISTR_v1.py \
     --in-file allele_freqs_test.txt \
     --out-file test_results.txt 
   ```

# Folder Contents
* `simulations/`: Contains scripts to run simulations to generate lookup tables
* `sistr/`: Contains scripts to run SISTR