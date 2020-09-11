# Script to generate LRT lookup table on TSCC

# Imports
import sys 
sys.path.append("/projects/ps-gymreklab/bonnieh/helper_functions")
from Simulation_functions import *
from ABC_functions import *

# Main function
def main():    
    
    # Load parameters from command line, see README.md for argument explanations
    per = int(sys.argv[1])
    opt_allele = int(sys.argv[2])
    numfile = int(sys.argv[3])
    s_vals = sys.argv[4]
    LRT_num_sims = int(sys.argv[5])
    use_var_gens = int(sys.argv[6])
    s_list = [float(s) for s in s_vals.split(',')]
    
    # Name and open output files
    outFolder1 = '20k_0825/'
    outFolder2 = '50k_0825/'
    outFolder3 = 'eurodem_0825/'
    outFile1 = '/projects/ps-gymreklab/bonnieh/lrt/results/' + outFolder1 + str(per) + '_' + str(opt_allele) + '_' + str(numfile) 
    outFile2 = '/projects/ps-gymreklab/bonnieh/lrt/results/' + outFolder2 + str(per) + '_' + str(opt_allele) + '_' + str(numfile) 
    outFile3 = '/projects/ps-gymreklab/bonnieh/lrt/results/' + outFolder3 + str(per) + '_' + str(opt_allele) + '_' + str(numfile) 
   
    outFile1a = outFile1 + '_het.txt'
    outFile1b = outFile1 + '_common.txt'
    outFile1c = outFile1 + '_bins.txt'
    outFile1d = outFile1 + '_freqs.txt'
   
    outFile2a = outFile2 + '_het.txt'
    outFile2b = outFile2 + '_common.txt'
    outFile2c = outFile2 + '_bins.txt'
    outFile2d = outFile2 + '_freqs.txt'
    
    outFile3a = outFile3 + '_het.txt'
    outFile3b = outFile3 + '_common.txt'
    outFile3c = outFile3 + '_bins.txt'
    outFile3d = outFile3 + '_freqs.txt'
    outFile3e = outFile3 + '_gens.txt'
    
    results1a = open(outFile1a, "w")
    results1b = open(outFile1b, "w")
    results1c = open(outFile1c, "w")
    results1d = open(outFile1d, "w")
    
    results2a = open(outFile2a, "w")
    results2b = open(outFile2b, "w")
    results2c = open(outFile2c, "w")
    results2d = open(outFile2d, "w")
    
    results3a = open(outFile3a, "w")
    results3b = open(outFile3b, "w")
    results3c = open(outFile3c, "w")
    results3d = open(outFile3d, "w")
    results3e = open(outFile3e, "w")
    
    results_list = []
    results_list.append(results1a)
    results_list.append(results1b)
    results_list.append(results1c)
    results_list.append(results1d)
    
    results_list.append(results2a)
    results_list.append(results2b)
    results_list.append(results2c)
    results_list.append(results2d)
    
    results_list.append(results3a)
    results_list.append(results3b)
    results_list.append(results3c)
    results_list.append(results3d)
    results_list.append(results3e)
    
    # Write results header
    results1a.write("s" + "\t" + "het" +"\n")
    results1b.write("s" + "\t" + "common" +"\n")
    results1c.write("s" + "\t" + "bins" +"\n")
    results1d.write("s" + "\t" + "freqs" +"\n")
    
    results2a.write("s" + "\t" + "het" +"\n")
    results2b.write("s" + "\t" + "common" +"\n")
    results2c.write("s" + "\t" + "bins" +"\n")
    results2d.write("s" + "\t" + "freqs" +"\n")
    
    results3a.write("s" + "\t" + "het" +"\n")
    results3b.write("s" + "\t" + "common" +"\n")
    results3c.write("s" + "\t" + "bins" +"\n")
    results3d.write("s" + "\t" + "freqs" +"\n")
    results3e.write("s" + "\t" + "num_gens" +"\n")
    
    # Mutation model parameters
    period_info = {}

    L1_log = 0.04
    L2_log = 0.15 
    L3_log = 0.33
    L4_log = 0.45 

    # List contents: mu, beta, p, l, optimal ru for the mu value
    period_info[1] = [10**-4.2, 0.5, 1, L1_log, 13]
    period_info[2] = [10**-5, 0.3, 0.6, L2_log, 6]
    period_info[3] = [10**-5.5, 0.3, 0.9, L3_log, 5] 
    period_info[4] = [10**-6, 0.3, 0.9, L4_log, 3]

    num_alleles = 25
    N_e = 7310
    max_iter = 55920
    end_samp_n = 6500
    
    # Get list of TMRCA values and put in TMRCA_list
    TMRCAFile = '/projects/ps-gymreklab/bonnieh/TMRCA.txt'
    TMRCA_file = open(TMRCAFile, 'r')
    TMRCA_list = []
    for line in TMRCA_file:
        TMRCA = line.strip()
        TMRCA = float(TMRCA)
        TMRCA = int(TMRCA)
        if TMRCA > 5920:
            TMRCA_list.append(TMRCA)
    
    # Get mutation model parameters based on STR class inputted on command line
    log_mu_prime = np.log10(period_info[per][0])+period_info[per][3]*(opt_allele - period_info[per][4])
    mu_prime = 10**log_mu_prime
    if mu_prime < 10**-8: mu_prime = 10**-8
    if mu_prime > 10**-3: mu_prime = 10**-3
            
    mu = mu_prime
    beta = period_info[per][1]
    p = period_info[per][2]
    L = period_info[per][3]
    
    # Index from which to obtain TMRCA values
    index = -1 + numfile*len(s_list)*LRT_num_sims
    
    # Run ABC_num_sims number of simulations for each s value
    for s in s_list:
        for result_file in results_list:
            result_file.write(str(s) + "\t")
        
        for i in range(0, LRT_num_sims):
            index = index + 1
            if use_var_gens == 0:
                max_iter = TMRCA_list[index]
                
            # Simulate allele frequencies
            allele_freqs_20k, allele_freqs_50k, allele_freqs_euro = Simulate(num_alleles, N_e, mu, beta, p, L, s, max_iter, end_samp_n)

            # Compute summary statistics of simulated allele frequencies
            het_20k = 1-sum([item**2 for item in allele_freqs_20k])
            het_50k = 1-sum([item**2 for item in allele_freqs_50k])
            het_euro = 1-sum([item**2 for item in allele_freqs_euro]) 

            common_20k = (allele_freqs_20k>=0.05).sum()
            common_50k = (allele_freqs_50k>=0.05).sum()
            common_euro = (allele_freqs_euro>=0.05).sum()

            bins_20k = GetBins(allele_freqs_20k, 5)
            bins_50k = GetBins(allele_freqs_50k, 5)
            bins_euro = GetBins(allele_freqs_euro, 5)
            
            # Write summary statistics and allele frequencies to file
            results1a.write(str(het_20k))
            results1b.write(str(common_20k))
            results1c.write(','.join(str(item) for item in bins_20k))
            results1d.write(','.join(str(item) for item in allele_freqs_20k))
            
            results2a.write(str(het_50k))
            results2b.write(str(common_50k))
            results2c.write(','.join(str(item) for item in bins_50k))
            results2d.write(','.join(str(item) for item in allele_freqs_50k))
            
            results3a.write(str(het_euro))
            results3b.write(str(common_euro))
            results3c.write(','.join(str(item) for item in bins_euro))
            results3d.write(','.join(str(item) for item in allele_freqs_euro))
            results3e.write(str(max_iter))
            
            if i == LRT_num_sims - 1:
                for result_file in results_list:
                    result_file.write('\n')
            else:
                for result_file in results_list:
                    result_file.write(';')
                
    for result_file in results_list:    
        result_file.close()
    
if __name__ == '__main__':
    main()