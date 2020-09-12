# Selection Inference on Short Tandem Repeats (SISTR): A method to measure negative selection at short tandem repeats

# Folder Contents
* `helper_functions/`: Contains helper functions used in scripts in `simulations/` and `SISTR/`
* `simulations/`: Contains scripts to generate lookup tables
* `sistr/`: Contains scripts to run SISTR

# Pipeline
1. Preprocessing:  
Simulate allele frequencies using mutation model, selection model, and forward simulation algorithm to generate lookup tables used later in SISTR.   
See `simulations/`  

2. SISTR:  
Input -> allele frequencies  
Output -> posterior estimate of selection coefficient  
See `sistr/`