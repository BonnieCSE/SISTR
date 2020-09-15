# Script to generate ABC lookup table

# Imports
import argparse
from Simulation_functions import *

# Main function
def main():
    
    # Load arguments from command line
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--out-folder")
    parser.add_argument("--per", type=int)
    parser.add_argument("--opt-allele", type=int)
    parser.add_argument("--num-sims", type=int, default=10000)
    parser.add_argument("--a-param", type=float, default=0.0881)
    parser.add_argument("--b-param", type=float, default=0.2541)
    parser.add_argument("--file-name-custom", default='')
    parser.add_argument("--use-var-gens", default = 'n')
   
    args = parser.parse_args()
    
    # Name and open output files
    
    outFile = args.out_folder + str(args.per) + '_' + str(args.opt_allele) + str(args.file_name_custom) + '.txt'
    
    results = open(outFile, "w")
    
    # Write results header
    results.write("s" + "\t" + "freqs" + "\n")

    # Mutation model parameters
    period_info = {}
    
    L2_log = 0.15 
    L3_log = 0.33
    L4_log = 0.45 

    # List contents: mu, beta, p, l, optimal ru for the mu value
    period_info[2] = [10**-5, 0.3, 0.6, L2_log, 6]
    period_info[3] = [10**-5.5, 0.3, 0.9, L3_log, 5] 
    period_info[4] = [10**-6, 0.3, 0.9, L4_log, 3]

    number_alleles = 25
    n_effec = 7310
    max_iter = 55920
    end_samp_n = 6500
    
    # Get list of TMRCA values and put in TMRCA_list
    if args.use_var_gens == 'y':
        TMRCAFile = './../lookup_tables/TMRCA.txt'
        TMRCA_file = open(TMRCAFile, 'r')
        TMRCA_list = []
        for line in TMRCA_file:
            TMRCA = line.strip()
            TMRCA = float(TMRCA)
            TMRCA = int(TMRCA)
            if TMRCA > 5920:
                TMRCA_list.append(TMRCA)
    
    # Get mutation model parameters based on STR class inputted on command line
    log_mu_prime = np.log10(period_info[args.per][0])+period_info[args.per][3]*(args.opt_allele - period_info[args.per][4])
    mu_prime = 10**log_mu_prime
    if mu_prime < 10**-8: mu_prime = 10**-8 
    if mu_prime > 10**-3: mu_prime = 10**-3

    mu = mu_prime
    beta = period_info[args.per][1]
    p_param = period_info[args.per][2]
    L = period_info[args.per][3]
    
    # Index from which to obtain TMRCA values
    index = -1
    
    # Run ABC_num_sims number of simulations
    for i in range(0, args.num_sims):

        # Draw s from prior
        # Assume only one possible prior value for other parameters (mu, beta, p, l)
        s = np.random.gamma(args.a_param, args.b_param)
       
        if s > 1:
            s = 1
            
        index = index + 1
        if args.use_var_gens == 'y':
            max_iter = TMRCA_list[index]
            
        # Simulate allele frequencies
        allele_freqs_20k, allele_freqs_50k, allele_freqs_euro = Simulate(number_alleles, n_effec, mu, beta, p_param, L, s, max_iter, end_samp_n)

        # Write allele frequencies to file
        results.write(str(s) + "\t" + ','.join(str(item) for item in allele_freqs_euro) + "\n")
        
    results.close()

if __name__ == '__main__':
    main()