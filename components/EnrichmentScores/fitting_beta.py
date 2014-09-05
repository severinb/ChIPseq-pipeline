import re, os
import numpy as np
from enrichment_score import calculate_enrichment_scores


def sum_of_posteriors_foreground_regions(fname):    
    posteriors = {}
    total_posterior = 0.
    with open(fname) as file_handler:
        for line in file_handler:            
            row = line.split()
            posterior = float(row[2])
            if re.search('_reg\d+', row[-1]):
                posteriors.setdefault(row[-1].strip(), 0.0)
                posteriors[row[-1].strip()] += posterior
            total_posterior += posterior
    return posteriors, total_posterior



def likelihood_derivative_beta(normalized_sitecount, beta, N_over_L, B):
    return beta + N_over_L - B*np.power(np.sum(1.0 / (beta + normalized_sitecount)) , -1 )


def fit_beta(siteFile, interm_dir, wmFile, number_of_windows):
    binding_regions_sitecount, total_sitecount = sum_of_posteriors_foreground_regions(siteFile)
    total_length = np.sum([l for l in number_of_windows.values()])
    M = len(number_of_windows.keys())
    beta_min, beta_max = 1.0e-12, 1.0e+0
    step_size = 1.0e-3
    
    # motifName = os.path.basename(wmFile)    
    # res_file = open(os.path.join(interm_dir, '%s.beta' % motifName), 'w')    
    beta_lower_bound, beta_upper_bound = beta_min, beta_min+step_size
    try:
        N_over_L = total_sitecount / float(total_length)
    except ZeroDivisionError:
        print wmFile
        exit()
    B = len(binding_regions_sitecount.keys())
    tmp = [] # it's defined to speed up the analysis
    # instead of calculating \sum_{r \in B} \frac{l_r}{beta*l_r + n_r}
    # I calculate: \sum_{r \in B} \frac{1}{beta + n_r/l_r}
    # where n_r/l_r is like site per nucleotide for a sequence
    for region, sitecount in binding_regions_sitecount.items():
        tmp.append(sitecount / float(number_of_windows[region]))
    normalized_sitecount = np.array(tmp)
    # # updating the lower and upper bounds for the value of Beta
    # # It's where the two bounds (upper and lower) have opposite signs
    # likelihood_beta_lower = likelihood_derivative_beta(normalized_sitecount, beta_lower_bound, N_over_L, B, number_of_windows)
    # # res_file.write('%f\t%0.12f\n' % (beta_lower_bound, likelihood_beta_lower))
    # if likelihood_beta_lower > .0:
    #     while True:            
    #         likelihood_beta_upper = likelihood_derivative_beta(normalized_sitecount, beta_upper_bound, N_over_L, B)
    #         if likelihood_beta_upper < 0.:
    #             beta_lower_bound = beta_upper_bound - step_size
    #             break
    #         beta_upper_bound += step_size
    #         # res_file.write('%f\t%0.12f\n' % (beta_upper_bound, likelihood_beta_upper))            
    #         if beta_upper_bound > beta_max:  # prevents to go very very large numbers for Beta
    #             return beta_max
    # else:  # once the derivative of likelihood for the lowest beta is already negative
    #     return beta_lower_bound

    # bisection method for refining the search area
    while (beta_upper_bound - beta_lower_bound) > 1.0e-12:
        mid_beta = (beta_upper_bound + beta_lower_bound) / 2.0        
        # print beta_lower_bound, mid_beta, beta_upper_bound
        mid_likelihood = likelihood_derivative_beta(normalized_sitecount, mid_beta, N_over_L, B)
        # print beta_lower_bound, mid_beta, beta_upper_bound, mid_likelihood
        if mid_likelihood > .0:
            beta_lower_bound = mid_beta
        else:
            beta_upper_bound = mid_beta

    # with open('tmp_enrichment', 'w') as outf:
    #     for beta in np.linspace(beta_min, 3.8201662787927147e-05, 200):        
    #         scores =calculate_enrichment_scores(siteFile, beta, number_of_windows, '')
    #         outf.write('%f\t%f\n' % (beta, scores['mean']))
    #         print beta, '\t', scores
    # print mid_beta
    # res_file.write('beta:\t%0.12f\n' % mid_beta)
    # res_file.close()
    return mid_beta
