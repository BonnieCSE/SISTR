# This file contains helper functions for performing the likelihood ratio test (LRT) for SISTR
# to test whether a model with selection fits better than a model without selection.

########## Imports ##########

from ABC_functions import *
from scipy.stats.distributions import chi2

########## LRT Helper Functions ##########

# Survival function for mixture distribution
def SF(x):
    if x > 0:
        return 0
    if x <= 0:
        return 1
    
"""Return list of summary statistics for given s

Parameters
----------
file_prefix: string
suffix_list: list
    List of summary statistics to return
s: float
    S value to return summary statistics for
    
Returns
-------
list_of_lists: list
    List of list of summary statistics
"""
def GetAllLRTLists(file_prefix, suffix_list, s):
    
    list_of_lists = []
    for suffix in suffix_list:
        
        lrt_file = open(file_prefix + '_' + suffix + '.txt', 'r')
    
        header = lrt_file.readline().strip()
    
        for line in lrt_file:
            info = line.strip().split('\t')
            s_val = float(info[0])
            
            if s_val == s:
                summ_stats = info[1]
            
                if suffix == 'bins': 
                    summ_stats= [s for s in summ_stats.split(';')]
                elif suffix == 'common':
                    summ_stats= [int(s) for s in summ_stats.split(';')]
                else:
                    summ_stats= [float(s) for s in summ_stats.split(';')]
            
                list_of_lists.append(summ_stats)
                break
           
        lrt_file.close()
    
    return list_of_lists

# Out of all s values in list, get s value closest to given s
def getNearestS(s_ABC_round, s_list_available):
    min_dist = 100000000
    nearest_s = -2
    for elem in s_list_available:
        dist = abs(s_ABC_round - elem)
        if dist < min_dist:
            min_dist = dist
            nearest_s = elem
    return nearest_s

# Get LRT bin for given s
def get_LRT_bin(s):
    if s < 0.00001:
        return 0
    if s <= 0.0001:
        return round(s, 5)
    if s <= 0.001:
        return round(s, 4)
    if s <= 0.01:
        return round(s, 3)
    if s <= 0.1:
        return round(s,2)
    return round(s,1)

# Get list of frequencies for given s value
def GetLRTListFreq(lrtFile, s_ABC):
    lrt_file = open(lrtFile, 'r')
    header = lrt_file.readline().strip()
    lrt_list = []
    
    for line in lrt_file:
        info = line.strip().split('\t')
        s = float(info[0])
        if s == s_ABC:
            freqs = info[1]
            freqs_list = [freq for freq in freqs.split(';')]
            lrt_file.close()
            return freqs_list
       
    # Return empty list if given s value is not found in file
    lrt_file.close()
    freqs_list = []
    return freqs_list

# Get list of 200 frequencies for given s value
def GetLRTListFreq200(lrtFile, s_ABC, return_first_200=True):
    lrt_file = open(lrtFile, 'r')
    header = lrt_file.readline().strip()
    lrt_list = []
    
    for line in lrt_file:
        info = line.strip().split('\t')
        s = float(info[0])
        if s == s_ABC:
            freqs = info[1]
            freqs_list = [freq for freq in freqs.split(';')]
            lrt_file.close()
            if return_first_200==True:
                return freqs_list[0:200]
            else:
                return freqs_list[200:400]
       
    # Return empty list if given s value is not found in file
    lrt_file.close()
    freqs_list = []
    return freqs_list

# Get list of frequencies from certain row in file
def GetLRTListByRow(lrtFile, row_num):
    lrt_file = open(lrtFile, 'r')
    header = lrt_file.readline().strip()
    lrt_list = []
    i = 0
    for line in lrt_file:
        
        if row_num == i:
            info = line.strip().split('\t')
            s = float(info[0])
            freqs = info[1]
            freqs_list = [freq for freq in freqs.split(';')]
            lrt_file.close()
            return freqs_list
        else:
            i = i + 1
       
    # Return empty list if given s value is not found in file
    lrt_file.close()
    freqs_list = []
    return freqs_list

# Get likelihood that observed summary statistics match given value of s
def GetLikelihoodFromTable(lookup_table_het, lookup_table_common, lookup_table_bins, num_sims, \
                           obs_het, obs_common, obs_bins, constant_het, denom_het, constant_common, \
                           denom_common, eps_bins, use_het, use_common, use_bins):
    
    # Get error tolerances to use when comparing summary statistics
    EPSILON_het = GetEpsilonHet(obs_het, constant_het, denom_het)
    EPSILON_common = GetEpsilonCommon(obs_common, constant_common, denom_common)
    EPSILON_bins = eps_bins
    
    num_accepted = 0
    
    # List of summary statistics to use 
    # (indicated by 0, 1, and 2 for heterozygosity, number of common alleles, allele frequency bins respectively)
    stats_to_check = []
    
    if use_het == 'y':
        stats_to_check.append(0)
    if use_common == 'y':
        stats_to_check.append(1)
    if use_bins == 'y':
        stats_to_check.append(2)
        
    for i in range(0, len(lookup_table_het)):
        
        stats = [False, False, False]
        if abs(lookup_table_het[i] - obs_het) < EPSILON_het:
            stats[0] = True
        if abs(lookup_table_common[i] - obs_common) < EPSILON_common:
            stats[1] = True
        if GetVectorDistance(lookup_table_bins[i], obs_bins) < EPSILON_bins:
            stats[2] = True
            
        # Check whether to accept s
        append = True
        for elem in stats_to_check:
            if stats[elem] == False:
                append = False
        if append == True:
            num_accepted = num_accepted + 1
    
    # Get likelihood with pseudocount
    likelihood = (num_accepted + 1)/num_sims
    
    return likelihood

# Get p-value from LRT
def LikelihoodRatioTest(LRT_table_0_het, LRT_table_s_het, LRT_table_0_common, LRT_table_s_common, \
                        LRT_table_0_bins, LRT_table_s_bins, LRT_num_sims, obs_het, obs_common, \
                        obs_bins, constant_het, denom_het, constant_common, denom_common, eps_bins, \
                        use_het, use_common, use_bins):
    
    # Get likelihood s = 0 
    likelihood_0 = GetLikelihoodFromTable(LRT_table_0_het, LRT_table_0_common, LRT_table_0_bins, \
                                          LRT_num_sims, obs_het, obs_common, obs_bins, constant_het, \
                                          denom_het, constant_common, denom_common, eps_bins, use_het, \
                                          use_common, use_bins)
   
    # Get likelihood s = ABC_s
    likelihood_s_ABC = GetLikelihoodFromTable(LRT_table_s_het, LRT_table_s_common, LRT_table_s_bins, \
                                              LRT_num_sims, obs_het, obs_common, obs_bins, constant_het, \
                                              denom_het, constant_common, denom_common, eps_bins, use_het, \
                                              use_common, use_bins)
     
    # Calculate likelihood ratio
    LR = likelihood_0/likelihood_s_ABC
   
    # Calculate LogLR 
    LogLR = -2*np.log(LR)
  
    # LogLR ~ Mixture distribution (50% 0, 50% Chi-square (df=1))
    pval = 0.5*SF(LogLR) + 0.5*chi2.sf(LogLR, 1)
    return likelihood_0, likelihood_s_ABC, LR, LogLR, pval