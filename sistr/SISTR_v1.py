# Script to run SISTR: 
# Method to obtain posterior estimate of s and corresponding p value for each STR

# Imports 
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from matplotlib import pyplot as plt
from scipy.stats import geom
import copy
import sys
sys.path.append("/storage/BonnieH/selection_project/helper_functions")
from LRT_functions import *

### Main function ###
def main():
    
    # Load parameters
    constant_het = float(sys.argv[1])
    denom_het = int(sys.argv[2])
    constant_common = int(sys.argv[3])
    denom_common = int(sys.argv[4])
    eps_bins = float(sys.argv[5])
    inFile = sys.argv[6]
    use_het = sys.argv[7]
    use_common = sys.argv[8]
    use_bins = sys.argv[9]
    num_bins = int(sys.argv[10])
    abc_model = sys.argv[11]
    lrt_model = sys.argv[12]
    LRT_num_sims = int(sys.argv[13])
    
    # Name output files
    file_type = ''
    if inFile == '/storage/BonnieH/selection_project/ssc_files/0810/allele_freqs_filt_test.txt':
        file_type = '_test_0810'
    if inFile == '/storage/BonnieH/selection_project/ssc_files/0810/allele_freqs_filt.txt':
        file_type = '_all_per_0810'
        
    filename = str(constant_het) + "_" + str(denom_het) + "_" + str(eps_bins) + \
               "_" + use_het + use_common + use_bins + str(num_bins) + "_" + abc_model + file_type
    
    outFile = '/storage/BonnieH/selection_project/per_locus/SISTR_results_0825/' + filename + '.txt'
    figFile = '/storage/BonnieH/selection_project/per_locus/SISTR_results_0825/figs/' + filename + '.png'
    statsFile = '/storage/BonnieH/selection_project/per_locus/SISTR_results_0825/stats/' + filename + '.txt'
    
    # Open output files
    allele_freqs_file = open(inFile, 'r')
    stats_file = open(statsFile, 'w')
    results = open(outFile, "w")

    # Write results header
    results.write("chrom" + "\t" + "start" + "\t" + "end" + "\t" + "total" + "\t" + "period" + \
                  "\t" + "optimal_ru" + "\t" + "motif" + "\t" + "coding" + "\t"  + 'intron' + \
                  '\t' + 'UTR5' + '\t' + 'UTR3' + '\t' + 'promoter5kb' + '\t' + 'intergenic' + \
                  "\t" + "tss_gene" + "\t" + "het" + "\t" + "common" + "\t" + "bins" + "\t" + \
                  "ABC_s_median" + "\t" + "ABC_s_95%_CI" + "\t" + "Num_s_accepted_(max_10000)" + "\t" + \
                  "Likelihood_0" + "\t" + "Likelihood_s" + "\t" + "LR" + "\t" + "LogLR" + \
                  "\t" + "LRT_p_value" + "\n")

    total_lines_acc = 0 # Total lines with >10 s accepted
    total_lines = 0 
    lower_0 = 0 # Total lines where lower bound is 0
    s_acc = [] # List of percent of s accepted
     
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
            file = '/gymreklab-tscc/bonnieh/abc/results/' + abc_model + '/' + str(period) + '_' + str(opt_allele) + '.txt' 
            table = GetABCList(file, num_bins)

            ABC_tables[per_opt] = table
        
    # Perform ABC on each locus
    for line in allele_freqs_file:
        
        # Get information from line
        total_lines = total_lines + 1
        info = line.strip().split('\t')
        chrom = info[0]
        start = int(info[1])
        end = int(info[2])
        freqs = info[3]
        total = info[4]
        per = int(info[5])
        motif = info[6]
        coding = info[7]
        intron = info[8]
        UTR5 = info[9]
        UTR3 = info[10]
        promoter5kb = info[11]
        intergenic = info[12]
        gene = info[13]

        # Get optimal allele and allele_freqs
        opt_allele, allele_freqs = Process_Freqs(freqs, per, end, start)
        freq_string = ','.join(str(round(item, 5)) for item in allele_freqs)
        
        # Add 0s to allele frequency list if number of alleles less than number of bins
        if len(allele_freqs) < num_bins:
            num_zeros_to_add = int((num_bins - len(allele_freqs))/2)
            for i in range(0, num_zeros_to_add):
                freq_string = '0.0,' + freq_string
                freq_string = freq_string + ',0.0'
                
        # Get summary stats
        obs_het, obs_common, obs_bins = GetSummStats(freq_string, num_bins)
        
        results.write(chrom + '\t' + str(start) + '\t' + str(end) + '\t' + total + '\t' + str(per) + \
                      '\t' + str(opt_allele) + '\t' + motif + '\t' + coding + '\t'  + intron + \
                      '\t' + UTR5 + '\t' + UTR3 + '\t' + promoter5kb + '\t' + intergenic + \
                      '\t' + gene + '\t' + str(round(obs_het, 7)) + '\t' + str(obs_common) + '\t' + \
                      ','.join(str(round(item,4)) for item in obs_bins) + '\t')
        
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
                                       obs_het, obs_common, obs_bins, constant_het, 
                                       denom_het, constant_common, denom_common, eps_bins, use_het, 
                                       use_common, use_bins)
        
        # Write N/A if <10 s accepted during ABC
        if s_ABC == -1:
            s_ABC = 'N/A'
            ABC_conf_int = 'N/A'
            num_accepted = '<10'
            likelihood_0 = 'N/A'
            likelihood_s_ABC = 'N/A'
            LR = 'N/A'
            LogLR = 'N/A'
            pval = 'N/A'
        
            results.write(str(s_ABC) + '\t' + ABC_conf_int + '\t' + str(num_accepted) + '\t')
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

            if lower_bound == 0:
                lower_0 = lower_0 + 1

            if lower_bound != -1:
                total_lines_acc = total_lines_acc + 1

            if num_accepted >= 10:
                s_acc.append(num_accepted)

            ABC_conf_int = '(' + str(lower_bound) + ' , ' + str(upper_bound) + ')'

            results.write(str(s_ABC) + '\t' + ABC_conf_int + '\t' + str(num_accepted) + '\t')

            s_ABC_round = get_LRT_bin(s_ABC)

            # Perform LRT

            ### Get list of s values in LRT simulations file ###
            s_list_available = []
            lrtFile = '/gymreklab-tscc/bonnieh/lrt/results/' + lrt_model + '/' \
            + str(per) + '_' + str(opt_allele) + '_freqs.txt' 
            
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
            lrtFile_for_s_0 = '/gymreklab-tscc/bonnieh/lrt/results/' + lrt_model + '/' + str(per) + '_' + str(opt_allele) + '_15_freqs.txt' 
            freqs_list_raw_0 = GetLRTListByRow(lrtFile_for_s_0, 0)
            LRT_table_0_het = []
            LRT_table_0_common = []
            LRT_table_0_bins = []

            # Get summary statistics from allele frequencies
            for freq_string in freqs_list_raw_0:

                obs_het_0, obs_common_0, obs_bins_0 = GetSummStats(freq_string, num_bins)
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

                obs_het_s, obs_common_s, obs_bins_s = GetSummStats(freq_string, num_bins)
                LRT_table_s_het.append(obs_het_s) 
                LRT_table_s_common.append(obs_common_s) 
                LRT_table_s_bins.append(obs_bins_s)

            # Perform LRT
            likelihood_0, likelihood_s_ABC, LR, LogLR, pval = LikelihoodRatioTest(LRT_table_0_het, \
                                    LRT_table_s_het, LRT_table_0_common, LRT_table_s_common, \
                                    LRT_table_0_bins, LRT_table_s_bins, LRT_num_sims, \
                                    obs_het, obs_common, obs_bins, constant_het, denom_het, \
                                    constant_common, denom_common, eps_bins, use_het, use_common, use_bins)

            results.write(str(round(likelihood_0, 7)) + '\t' + str(round(likelihood_s_ABC, 7)) + '\t' + \
                              str(round(LR, 7)) + '\t' + str(round(LogLR ,7)) + '\t' + str(round(pval, 7)) + '\n')

    allele_freqs_file.close()
            
    results.close()
    
    # Output statistics
    stats_file.write("Total lines: " + str(total_lines) + "\n")
    percent_acc = total_lines_acc/total_lines
    stats_file.write("Total accepted lines: " + str(total_lines_acc) + " " + str(percent_acc) + "\n")
    percent_0 = lower_0/total_lines
    stats_file.write("Lower bound for s is zero: " + str(lower_0) + " " + str(percent_0) + "\n")
    
    for i in range(0, len(s_acc)):
        s_acc[i] = s_acc[i]/10000 * 100
    
    s_acc_med = np.median(s_acc)
    s_acc_mean = np.mean(s_acc)

    stats_file.write("s_acc median: " + str(s_acc_med) + "\n")
    stats_file.write("s_acc mean: " + str(s_acc_mean) + "\n")
 
    plt.figure(1)
    buckets = list(np.arange(0, 101, 1)) 
    
    plt.hist(s_acc, bins=buckets, weights=np.ones(len(s_acc)) / len(s_acc))

    plt.title('Percent s accepted during ABC %s model \n const_het %.5f, denom_het = %d, num_bins = %d \n Summ stats used: het = %s common alleles = %s bins = %s'%(abc_model, constant_het, denom_het, num_bins, use_het, use_common, use_bins))
    plt.xlabel("% s accepted")
    plt.ylabel("Frequency")
    plt.savefig(figFile, bbox_inches='tight')
    
if __name__ == '__main__':
    main()