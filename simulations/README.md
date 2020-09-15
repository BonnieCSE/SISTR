# Running simulations to generate lookup tables

**Note: Due to the compute time required to generate the lookup tables, we recommend downloading the tables used in Mitra et. al. 2020, which are available for community use and can be found [here](https://drive.google.com/drive/folders/1p_QoSQ7gzs7hVEwfJyhGMT-ELXZcWoA_?usp=sharing).** However, if you would like to create your own custom lookup tables, see instructions below.

## Command to generate an ABC lookup table  
Note: If you choose to generate your own lookup tables, you must generate 23 ABC lookup tables, one for each period/optimal allele combination below.    
List of optimal alleles per period to simulate: period 2, optimal alleles 11-20; period 3, optimal alleles 5-13; period 4, optimal alleles 7-10  

Example command:
```
python ABC_lookup.py \
     --out-folder ./../sistr_resources_test/abc_lookup/ \
     --per 3 \
     --opt-allele 9 \
     --num-sims 5
```

Required parameters:  
* __`--out-folder <string>`__ path to an existing folder to which to output lookup table
* __`--per <int>`__ period (repeat unit length) -> Note: SISTR currently only handles TRs with repeat unit lengths 2-4 bp  
* __`--opt-allele <int>`__ optimal allele   

Optional parameters (default values are those used in Mitra, et. al. 2020):
* __`--num_sims <int>`__ number of simulations  (default value: 10000)
* __`--a_param <float>`__ a (shape) parameter of prior gamma distribution (default value: 0.0881)  
* __`--b_param <float>`__ b (scale) parameter of prior gamma distribution (default value: 0.2541)
* __`--num_gens <int>`__ number of generations to run the simulations for (default: 55920)
* __`--file-name-custom <string>`__ customize ending of file name (default: no extra characters in file name) 

### Output file format:
The output file will be a tab delimited file with 2 columns: s (s value used), freqs (corresponding allele frequencies)  
Each row is a separate simulation.  
The file is named in the format \[per\]_\[opt-allele\]\[file-name-custom\].txt  

## Command to generate a LRT lookup table  
Note: If you choose to generate your own lookup tables, you must generate 46 LRT lookup tables, two for each period/optimal allele combination below. One lookup table contains a variety of s values from 0 to 1 inclusive, while the other table only contains information regarding s = 0.  
List of optimal alleles per period to simulate: period 2, optimal alleles 11-20; period 3, optimal alleles 5-13; period 4, optimal alleles 7-10  

Example command for lookup table with a variety of s values:
```
python LRT_lookup.py \
     --out-folder ./../sistr_resources_test/lrt_lookup/ \
     --per 3 \
     --opt-allele 9 \
     --s-vals 0,1e-06,1e-05,0.0001,0.001 \
     --num-sims 2  
```

Example command for lookup table with only s = 0 (Note: In this case, `--s-vals` must be the string `0,0` and `--file-name-addition` must be the string `_zero`):
```
python LRT_lookup.py \
     --out-folder ./../sistr_resources_test/lrt_lookup/ \
     --per 3 \
     --opt-allele 9 \
     --s-vals 0,0 \
     --num-sims 5 \
     --file-name-custom _zero  
```

Required parameters:  
* __`--out-folder <string>`__ path to an existing folder to which to output lookup table
* __`--per <int>`__ period (repeat unit length) -> Note: SISTR currently only handles TRs with repeat unit lengths 2-4 bp  
* __`--opt-allele <int>`__ optimal allele   
* __`--s_vals <string>`__ list of s values, separated by commas, for which to generate lookup table 

Optional parameters (default values are those used in Mitra, et. al. 2020):
* __`--num_sims <int>`__ number of simulations per s value (default: 2000)
* __`--num_gens <int>`__ number of generations to run the simulations for (default: 55920)
* __`--file-name-custom <string>`__ customize ending of file name (default: no extra characters in file name)  
Note: The flag __`--file-name-custom _zero`__ is required when generating a lookup table with only s = 0

### Output file format:
The output will be a tab delimited file with 2 columns: s (s value used), freqs (num-sims allele frequencies seqparated by semicolons)  
Each row is a separate s value.  
The file is named in the format \[per\]_\[opt-allele\]\[file-name-custom\]_freqs.txt