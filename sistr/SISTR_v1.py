# Script to run SISTR version 1: 
# A method to obtain a posterior estimate of s and corresponding p value for each STR given allele frequency data

# Imports 
import copy
import argparse
from scipy.stats import geom
from LRT_functions import *

### Main function ###
def main():
    
    # Load arguments from command line
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--eps-het-numerator", type=float, default=0.005)
    parser.add_argument("--eps-het-denominator", type=int, default=3)
    parser.add_argument("--eps-bins", type=float, default=0.3)
    parser.add_argument("--num-bins", type=int, default=5)
    parser.add_argument("--lrt-num-sims", type=int, default=2000)
    parser.add_argument("--abc-lookup-folder", default = 'sistr_resources/abc_lookup/')
    parser.add_argument("--lrt-lookup-folder", default = 'sistr_resources/lrt_lookup/')
    parser.add_argument("--in-file")
    parser.add_argument("--out-file")
    
    args = parser.parse_args()
    
    constant_common = 1
    denom_common = 1
    
    use_het = 'y'
    use_common = 'n'
    use_bins = 'y'
        
    # Open input and output files
    allele_freqs_file = open(args.in_file, 'r')
    header = allele_freqs_file.readline().strip().split('\t')
    results = open(args.out_file, "w")

    # Write results header
    results.write("chrom" + "\t" + "start" + "\t" + "end" + "\t" + "total" + "\t" + "period" + \
                  "\t" + "optimal_ru" + "\t" + "motif" + "\t" + "ABC_s_median" + "\t" + "ABC_s_95%_CI" + \
                  "\t" + "Percent_s_accepted" + "\t" + "Likelihood_0" + "\t" + "Likelihood_s" + "\t" + \
                  "LR" + "\t" + "LogLR" + "\t" + "LRT_p_value" + "\n")
     
    # Preprocess ABC lookup table
    # Get ABC tables
    ABC_tables = {}
    
    opt_allele_dic_w_per = {}
    opt_allele_dic_w_per[3] = [(3,5),(3,6),(3,7),(3,8),(3,9),(3,10),(3,11),(3,12),(3,13)]
    opt_allele_dic_w_per[2] = [(2,11),(2,12),(2,13),(2,14),(2,15),(2,16),(2,17),(2,18),(2,19),(2,20)]
    opt_allele_dic_w_per[4] = [(4,7),(4,8),(4,9),(4,10)]
    
    pers = [2,3,4]
    for per in pers:
        for per_opt in opt_allele_dic_w_per[per]:
            period = per_opt[0]
            opt_allele = per_opt[1]
            file = args.abc_lookup_folder + str(period) + '_' + str(opt_allele) + '.txt' 
            table = GetABCList(file, args.num_bins)

            ABC_tables[per_opt] = table
        
    # Perform ABC on each locus
    for line in allele_freqs_file:
        
        # Get information from line
        info = line.strip().split('\t')
        chrom = info[0]
        start = int(info[1])
        end = int(info[2])
        freqs = info[3]
        total = info[4]
        per = int(info[5])
        motif = info[6]

        # Get optimal allele and allele_freqs
        opt_allele, allele_freqs = Process_Freqs(freqs, per, end, start)
        freq_string = ','.join(str(round(item, 5)) for item in allele_freqs)
        
        # Add 0s to allele frequency list if number of alleles less than number of bins
        if len(allele_freqs) < args.num_bins:
            num_zeros_to_add = int((args.num_bins - len(allele_freqs))/2)
            for i in range(0, num_zeros_to_add):
                freq_string = '0.0,' + freq_string
                freq_string = freq_string + ',0.0'
                
        # Get summary stats
        obs_het, obs_common, obs_bins = GetSummStats(freq_string, args.num_bins)
        
        results.write(chrom + '\t' + str(start) + '\t' + str(end) + '\t' + total + '\t' + str(per) + \
                      '\t' + str(opt_allele) + '\t' + motif + '\t')
        
        # All optimal alleles > 13 have the same mutation rate as optimal allele 13
        if per == 3 and opt_allele > 13: 
            opt_allele = 13
        if per == 3 and opt_allele < 5:
            opt_allele = 5
            
        if per == 4 and opt_allele > 10:
            opt_allele = 10
        if per == 4 and opt_allele < 7:
            opt_allele = 7
        
        if per == 2 and opt_allele > 20:
            opt_allele = 20
        if per == 2 and opt_allele < 11:
            opt_allele = 11
            
        # Read abcFile line by line and place in lookup table in the form of a list
        abc_list = ABC_tables[(per, opt_allele)]
        
        # Perform ABC
        s_ABC, lower_bound, upper_bound, num_accepted, s_accepted = Get_S_ABC(abc_list, 
                                       obs_het, obs_common, obs_bins, args.eps_het_numerator, 
                                       args.eps_het_denominator, constant_common, denom_common, args.eps_bins, use_het, 
                                       use_common, use_bins)
        
        # Write N/A if <10 s accepted during ABC
        if s_ABC == -1:
            s_ABC = 'N/A'
            ABC_conf_int = 'N/A'
            num_accepted = 'N/A'
            likelihood_0 = 'N/A'
            likelihood_s_ABC = 'N/A'
            LR = 'N/A'
            LogLR = 'N/A'
            pval = 'N/A'
        
            results.write(str(s_ABC) + '\t' + ABC_conf_int + '\t' + num_accepted + '\t')
            results.write(str(likelihood_0) + '\t' + str(likelihood_s_ABC) + '\t' + \
                          str(LR) + '\t' + str(LogLR) + '\t' + str(pval) + '\n')
         
        # Perform LRT if can get posterior estimate of s
        else:
            # Round s values < 10^-5 to 0
            if s_ABC < 10**-5:
                s_ABC = 0
            else:
                s_ABC = round(s_ABC, 5)

            if lower_bound < 10**-5:
                lower_bound = 0
            else:
                lower_bound = round(lower_bound, 5)

            if upper_bound < 10**-5:
                upper_bound = 0
            else:
                upper_bound = round(upper_bound, 5)

            ABC_conf_int = '(' + str(lower_bound) + ' , ' + str(upper_bound) + ')'

            results.write(str(s_ABC) + '\t' + ABC_conf_int + '\t' + str(round(num_accepted/len(ABC_tables[(3,5)])*100,4)) + '\t')

            s_ABC_round = get_LRT_bin(s_ABC)

            # Perform LRT

            ### Get list of s values in LRT simulations file ###
            s_list_available = []
            lrtFile = args.lrt_lookup_folder + str(per) + '_' + str(opt_allele) + '_freqs.txt' 
            
            lrt_file = open(lrtFile, 'r')
            header = lrt_file.readline().strip()

            for line in lrt_file:
                info = line.strip().split('\t')
                s = float(info[0])
                if s not in s_list_available:
                    s_list_available.append(s)

            # Get nearest s in LRT file to perform LRT
            if s_ABC_round not in s_list_available:
                s_ABC_round = getNearestS(s_ABC_round, s_list_available)

            # Get LRT summary statistic tables for s = 0
            lrtFile_for_s_0 = args.lrt_lookup_folder + str(per) + '_' + str(opt_allele) + '_zero_freqs.txt' 
            freqs_list_raw_0 = GetLRTListByRow(lrtFile_for_s_0, 0)
            LRT_table_0_het = []
            LRT_table_0_common = []
            LRT_table_0_bins = []

            # Get summary statistics from allele frequencies
            for freq_string in freqs_list_raw_0:

                obs_het_0, obs_common_0, obs_bins_0 = GetSummStats(freq_string, args.num_bins)
                LRT_table_0_het.append(obs_het_0) 
                LRT_table_0_common.append(obs_common_0) 
                LRT_table_0_bins.append(obs_bins_0)

            # Get LRT summary statistic tables for s = s_ABC_round
            if s_ABC_round == 0:
                freqs_list_raw_s = GetLRTListByRow(lrtFile_for_s_0, 1)
            else:
                freqs_list_raw_s = GetLRTListFreq(lrtFile, s_ABC_round)
            
            LRT_table_s_het = []
            LRT_table_s_common = []
            LRT_table_s_bins = []

            # Get summary statistics from allele frequencies
            for freq_string in freqs_list_raw_s:

                obs_het_s, obs_common_s, obs_bins_s = GetSummStats(freq_string, args.num_bins)
                LRT_table_s_het.append(obs_het_s) 
                LRT_table_s_common.append(obs_common_s) 
                LRT_table_s_bins.append(obs_bins_s)

            # Perform LRT
            likelihood_0, likelihood_s_ABC, LR, LogLR, pval = LikelihoodRatioTest(LRT_table_0_het, \
                                    LRT_table_s_het, LRT_table_0_common, LRT_table_s_common, \
                                    LRT_table_0_bins, LRT_table_s_bins, args.lrt_num_sims, \
                                    obs_het, obs_common, obs_bins, args.eps_het_numerator, args.eps_het_denominator, \
                                    constant_common, denom_common, args.eps_bins, use_het, use_common, use_bins)

            results.write(str(round(likelihood_0, 7)) + '\t' + str(round(likelihood_s_ABC, 7)) + '\t' + \
                              str(round(LR, 7)) + '\t' + str(round(LogLR ,7)) + '\t' + str(round(pval, 7)) + '\n')

    allele_freqs_file.close()
            
    results.close()
    
if __name__ == '__main__':
    main()