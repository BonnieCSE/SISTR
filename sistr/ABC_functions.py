# This file contains helper functions for performing ABC (approximate Bayesian computation) in SISTR
# to obtain a posterior distribution of s (selection coefficient) for each locus.

########## Imports ##########

import math
import numpy as np

########## ABC Helper Functions ##########

### Get bins of allele frequencies summary statistic ###
def GetBins(allele_freqs, num_bins):
    
    bins = [0] * num_bins # List of binned allele frequencies
   
    middle_index = int(len(allele_freqs)/2)
    
    boundary_low = middle_index - int((num_bins - 1)/2) # Everything below boundary_low is combined together into a bin
    bins[0] = sum(allele_freqs[0:boundary_low + 1])
    
    boundary_high = middle_index + int((num_bins - 1)/2) # Everything above boundary_high is combined together into a bin
    bins[num_bins - 1] = sum(allele_freqs[boundary_high:len(allele_freqs)])
    
    bins_index = 1
    
    # Fill in rest of the bins between lower and upper boundary
    for i in range(boundary_low + 1, boundary_high):
        bins[bins_index] = allele_freqs[i]
        bins_index = bins_index + 1
                           
    return bins

### Get summary statistics ###
# Return heterozygosity, number of common alleles, bins of allele frequencies (if num_bins > 0)
def GetSummStats(freq_string, num_bins):
    allele_freqs = [float(freq) for freq in freq_string.split(',')]
    obs_het = 1-sum([item**2 for item in allele_freqs])
    obs_common = len([i for i in allele_freqs if i >= 0.05]) 
    if num_bins == 0:
        return obs_het, obs_common
    else:
        obs_bins = GetBins(allele_freqs, num_bins)
        return obs_het, obs_common, obs_bins

### Get average allele ###
def GetAvg(allele_freqs):

    allele_sizes = list(range(-1*int(len(allele_freqs)/2), int(len(allele_freqs)/2)+1))
    avg = np.dot(allele_freqs, allele_sizes)
    
    return avg

### Get variance ###
def GetVar(allele_freqs):

    allele_sizes = list(range(-1*int(len(allele_freqs)/2), int(len(allele_freqs)/2)+1))
    avg = np.dot(allele_freqs, allele_sizes)
    difference = np.zeros(len(allele_freqs))
    for i in range(0, len(allele_sizes)):
        difference[i] = (allele_sizes[i] - avg)**2

    return np.dot(allele_freqs, difference)

### Process allele frequencies ###
# Return optimal allele repeat units, allele frequencies
def Process_Freqs(freqs, per, end, start, return_freqs=True, motif=False):
    # Process freqs
    freqs_list = freqs.split(',')

    # Create dictionary of freqs: key is allele and value is freq
    freqs_dic = {}
    
    for pair in freqs_list:
        # Processing allele frequencies in non-motif (alternate) format
        if motif == False:
            allele = int(pair.split(':')[0])
            
        # Processing allele frequencies in motif format
        
        else:
            allele = pair.split(':')[0]
            allele = len(allele)
        freq = int(pair.split(':')[1])
        freqs_dic[int(allele/per)] = freq
        
    # Note: Actual population size is actually half the value of the variable pop_size 
    # because pop_size is calculated by adding up diploid allele freqs
    pop_size = 0
    for elem in freqs_dic:
        pop_size = pop_size + freqs_dic[elem]
   
    # Get optimal allele (how many ru away from reference allele) - allele with highest frequency
    opt_allele_rel = max(freqs_dic, key=freqs_dic.get)

    # Get info about reference allele (length in base pairs, repeat units)
    # Get optimal allele repeat unit length 
    if motif == False:
        ref_length_bp = end - start + 1
        ref_length_ru = int(ref_length_bp/per)
        opt_allele = ref_length_ru + opt_allele_rel
                         
    else:
        opt_allele = opt_allele_rel                
    
    if return_freqs == False:
        return opt_allele
    
    else:
        
        # Get allele frequencies
        freqs_dic_final = {}
        allele_list = []
        
        for allele in freqs_dic:
            new_allele = allele - opt_allele_rel
            freqs_dic_final[new_allele] = freqs_dic[allele]
            allele_list.append(new_allele)

        # Get highest absolute value in list
        max_allele = abs(max(allele_list, key=abs))

        # Put allele freqs in list
        allele_freqs = [0] * (2*max_allele+1) 
        for key in freqs_dic_final:
            allele_freq = freqs_dic_final[key]/pop_size
            allele_freqs[max_allele+key] = allele_freq

        return opt_allele, allele_freqs

### Get epsilon (error tolerance) for heterozygosity ###
def GetEpsilonHet(obs_het, constant_het, denom_het):
    epsilon = (obs_het + constant_het)/denom_het
    return epsilon

### Get epsilon (error tolerance) for number of common alleles ###
def GetEpsilonCommon(obs_common, constant_common, denom_common):
    epsilon = obs_common/denom_common + constant_common
    return math.floor(epsilon)

### Get distance between 2 vectors ###
def GetVectorDistance(vector1, vector2):
    distance = 0
    for i in range(0, len(vector1)):
        distance = distance + abs(vector1[i] - vector2[i])
    return distance

### Get list of s with corresponding summary statistics ###
def GetABCList(abcFile, num_bins):
    abc_file = open(abcFile, 'r')
    header = abc_file.readline().strip().split('\t')
    
    abc_list = []
    
    if num_bins == 0:
        for line in abc_file:
            info = line.strip().split('\t')
            s = float(info[0])
            het = float(info[1])
            common = int(info[2])
            stats_list = [s, het, common]
            abc_list.append(stats_list)
     
    # Get bins summary statistic
    else:
        # Get column number of freqs column in file
        freqs_column = 0
        for i in range(0, len(header)):
            if header[i] == 'freqs':
                freqs_column = i

        for line in abc_file:
            info = line.strip().split('\t')
            s = float(info[0])
            freq_string = info[freqs_column]
            allele_freqs = [float(freq) for freq in freq_string.split(',')]

            abc_het = 1-sum([item**2 for item in allele_freqs])
            abc_common = len([i for i in allele_freqs if i >= 0.05]) 
            abc_bins = GetBins(allele_freqs, num_bins)

            stats_list = [s, abc_het, abc_common, abc_bins]
            abc_list.append(stats_list)
        
    abc_file.close()
    return abc_list
        
# Get posterior estimate of s using ABC
def Get_S_ABC(abc_list, obs_het, obs_common, obs_bins, constant_het, 
              denom_het, constant_common, denom_common, eps_bins, use_het, 
              use_common, use_bins):
    
    s_accepted = []
    
    # Get error tolerances to use when comparing summary statistics
    EPSILON_het = GetEpsilonHet(obs_het, constant_het, denom_het)
    EPSILON_common = GetEpsilonCommon(obs_common, constant_common, denom_common)
    EPSILON_bins = eps_bins
  
    # List of summary statistics to use 
    # (indicated by 0, 1, and 2 for heterozygosity, number of common alleles, allele frequency bins respectively)
    stats_to_check = []
    
    if use_het == 'y':
        stats_to_check.append(0)
    if use_common == 'y':
        stats_to_check.append(1)
    if use_bins == 'y':
        stats_to_check.append(2)
    
    for combo in abc_list:
        stats = [False, False, False]
        if abs(obs_het - combo[1]) < EPSILON_het: 
            stats[0] = True
        
        if abs(obs_common - combo[2]) < EPSILON_common: 
            stats[1] = True
                
        if GetVectorDistance(obs_bins, combo[3]) < EPSILON_bins:
            stats[2] = True
            
        # Check whether to accept s
        append = True
        for elem in stats_to_check:
            if stats[elem] == False:
                append = False
        if append == True:
            s_accepted.append(combo[0])
         
    num_accepted = len(s_accepted)
    
    # Get posterior estimate of s along with 95% CI
    if num_accepted >= 10:
        median_s = np.median(s_accepted)
        lower_bound = np.percentile(s_accepted, 2.5) 
        upper_bound = np.percentile(s_accepted, 97.5) 
        return median_s, lower_bound, upper_bound, num_accepted, s_accepted
    
    else:
        return -1, -1, -1, 0, s_accepted