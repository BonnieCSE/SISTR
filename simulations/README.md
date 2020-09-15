# Running simulations to generate lookup tables

## Note: Lookup tables used in Mitra et. al. 2020 are available for community use and can be found at:
https://drive.google.com/drive/folders/1p_QoSQ7gzs7hVEwfJyhGMT-ELXZcWoA_?usp=sharing

## Command to generate ABC lookup table (generate one ABC lookup table for each period/optimal allele combination)
Note: If you choose to generate your own lookup tables, you must generate 23 ABC lookup tables, one for each period/optimal allele combination below  
List of optimal alleles per period to simulate: per 2, opt 11-20; per 3, opt 5-13; per 4, opt 7-10  

Example command:
```
python ABC_lookup.py \
     --out-folder ./../sistr_resources_test/abc_lookup/ \
     --per 3 \
     --opt-allele 9 \
     --num-sims 5
```

Required parameters:  
* `--out-folder <string>` path to folder to which to output lookup table
* `--per <int>` period (repeat unit length) -> Note: SISTR currently only handles TRs with repeat unit lengths 2-4 bp  
* `--opt-allele <int>` optimal allele   

Optional parameters (default values are those used in Mitra, et. al. 2020):
* `--num_sims <int>` number of simulations  (default value: 10000)
* `--a_param <float>` a (shape) parameter of prior gamma distribution (default value: 0.0881)  
* `--b_param <float>` b (scale) parameter of prior gamma distribution (default value: 0.2541)
* `--file-name-addition <string>` customize ending of file name (default: no extra characters in file name) 

### Output file format:
The output file will be a tab delimited file with 2 columns: s (s value used), freqs (corresponding allele frequencies)  
Each row is a separate simulation.  
The file is named in the format <per>_<opt-allele>.txt  

## Command to generate LRT lookup table (generate one LRT lookup table for each period/optimal allele combination)
Note: If you want to generate your own lookup tables, you must generate 46 LRT lookup tables, two for each period/optimal allele combination below. One lookup table contains a wide range of s values from 0 to 1 inclusive, and one lookup table contains information from simulations using s = 0.
List of optimal alleles per period to simulate: per 2, opt 11-20; per 3, opt 5-13; per 4, opt 7-10  

Example command for lookup table with a variety of s values:
```
python LRT_lookup.py \
     --out-folder ./../sistr_resources_test/lrt_lookup/ \
     --per 3 \
     --opt-allele 9 \
     --s-vals 0,1e-06,1e-05 \
     --num-sims 2  
```

Example command for lookup table with only s = 0 (Note: In this case, --s-vals must be the string 0,0 and --file-name-addition must be the string _zero):
```
python LRT_lookup.py \
     --out-folder ./../sistr_resources_test/lrt_lookup/ \
     --per 3 \
     --opt-allele 9 \
     --s-vals 0,0 \
     --num-sims 2 \
     --file-name-addition _zero  
```

Required parameters:  
* `--out-folder <string>` path to folder to which to output lookup table
* `--per <int>` period (repeat unit length) -> Note: SISTR currently only handles TRs with repeat unit lengths 2-4 bp  
* `--opt-allele <int>` optimal allele   
* `s_vals <string>` list of s values, separated by commas, for which to generate lookup table 

Optional parameters (default values are those used in Mitra, et. al. 2020):
* `--num_sims <int>` number of simulations per s value (default: 2000)
* `--file-name-addition <string>` customize ending of file name (default: no extra characters in file name)
Note: `--file-name-addition _zero` is required when generating lookup table for s = 0

### Output files:
The output will be a tab delimited file with 2 columns: s (s value used), freqs (num-sims allele frequencies seqparated by semicolons)  
Each row is a separate s value.  
The file is named in the format <per>_<opt-allele><file-name-addition>_freqs.txt