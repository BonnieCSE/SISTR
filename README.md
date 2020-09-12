# Selection Inference on Short Tandem Repeats (SISTR): A method to measure negative selection at short tandem repeats

# Folder Contents
* `helper_functions/`: Contains helper functions used in `simulations/` and `SISTR/` scripts
* `simulations/`: Contains scripts to generate lookup tables
* `sistr/`: Contains scripts to run SISTR

# Pipline
1. Preprocessing:  
Simulate allele frequencies using mutation model, selection model, forward simulation algorithm and generate lookup tables used later in SISTR.   
See `simulations/`  

2. SISTR:  
Input -> allele frequencies  
Output -> posterior estimate of selection coefficient  
See `sistr/`