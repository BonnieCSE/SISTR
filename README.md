# SISTR: Selection Inference on Short Tandem Repeats 
SISTR is a method to measure negative selection at short tandem repeats.  

It takes allele frequency data (per-locus frequencies of each allele length) as input and outputs a posterior estimate of the selection coefficient at each locus. 

[![DOI](https://zenodo.org/badge/294813574.svg)](https://zenodo.org/badge/latestdoi/294813574)

Preprint: https://www.biorxiv.org/content/10.1101/2020.03.04.974170v2

For questions on SISTR usage or setup, please contact Bonnie Huang (bbhuang@ucsd.edu)  

# Installation
SISTR uses Python3 and the following libraries in addition to the Python Standard Library: SciPy, NumPy

You can obtain SISTR by cloning the Github repository:

```
git clone https://github.com/BonnieCSE/SISTR
cd SISTR
```

# Pipeline
1. The first step is a preprocessing step which involves simulating allele frequencies to generate lookup tables used later in step 2. **Due to the compute time required to generate the lookup tables, we recommend using the following precomputed tables found [here](https://drive.google.com/drive/folders/1p_QoSQ7gzs7hVEwfJyhGMT-ELXZcWoA_?usp=sharing) by downloading the entire folder `sistr_resources/` and saving it in the same directory where you are running SISTR.** However, if you wish to create your own custom lookup tables, see `simulations/` for further details. 
   
   Example command to run simulations to generate a custom ABC lookup table:    
   ```
   python ./simulations/ABC_lookup.py \
     --out-folder sistr_resources_test/abc_lookup/ \
     --per 3 \
     --opt-allele 5 \
     --num-sims 5
   ```

2. The second step is to run SISTR, which requires (1) allele frequency data and (2) the lookup tables generated in step 1 as input. It outputs a posterior estimate of the selection coefficient at each locus. 

   See `sistr/` for further details.

   Example command to run SISTR:  
   ```
   python sistr/SISTR_v1.py \
     --in-file sistr/allele_freqs_test.txt \
     --out-file test_results.txt 
   ```

# Folder Contents
* `simulations/`: Contains scripts to run simulations to generate lookup tables
* `sistr/`: Contains scripts to run SISTR

# FAQ
Q: If I decide to create my own lookup tables, what are the main features I can customize ?  

A: Currently, for the ABC lookup tables, you can customize the number of simulations, the number of generations the simulations are run for, and the a (shape) and b (scale) parameters of the prior gamma distrbution from which s is drawn.  
For the LRT lookup tables, you can customize the number of simulations as well as the number of generations the simulations are run for.