import re, os
import numpy as np


def sum_of_posteriors_foreground_regions(fname):    
    posteriors = {}
    number_of_windows = {}
    number_of_regions = {}
    total_posterior = 0.
    total_number_of_windows = 0.
    with open(fname) as file_handler:
        for line in file_handler:            
            row = line.split()
            posterior = float(row[2])
            if re.search('_reg\d+', row[-1]):
                posteriors.setdefault(row[-1].strip(), 0.0)
                posteriors[row[-1].strip()] += posterior
                number_of_windows.setdefault(row[-1].strip(), 0)
                number_of_windows[row[-1].strip()] += 1
            number_of_regions.setdefault(row[-1].strip(), 0)
            total_posterior += posterior
            total_number_of_windows += 1.
    return posteriors, total_posterior, number_of_windows, total_number_of_windows, len(number_of_regions.keys())


def likelihood_derivative_beta(normalized_sitecount, beta, N_over_L, B):
    return beta + N_over_L - B*np.power(np.sum(1.0 / (beta + normalized_sitecount)) , -1 )


def fit_beta(siteFile, interm_dir, wmFile):
    binding_regions_sitecount, total_sitecount, number_of_windows, total_length, M = sum_of_posteriors_foreground_regions(siteFile)
    beta_min, beta_max = 1.0e-12, 1.0e+3
    step_size = 5.0
    motifName = os.path.basename(wmFile)    
    # res_file = open(os.path.join(interm_dir, '%s.beta' % motifName), 'w')    
    beta_lower_bound, beta_upper_bound = beta_min, beta_min+step_size
    N_over_L = total_sitecount / float(total_length)
    B = len(binding_regions_sitecount.keys())
    tmp = [] # it's defined to speed up the analysis
    # instead of calculating \sum_{r \in B} \frac{l_r}{beta*l_r + n_r}
    # I calculate: \sum_{r \in B} \frac{1}{beta + n_r/l_r}
    # where n_r/l_r is like site per nucleotide for a sequence
    for region, sitecount in binding_regions_sitecount.items():
        tmp.append(sitecount / float(number_of_windows[region]))
    normalized_sitecount = np.array(tmp)
    # updating the lower and upper bounds for the value of Beta
    # It's where the two bounds (upper and lower) have opposite signs
    likelihood_beta_lower = likelihood_derivative_beta(normalized_sitecount, beta_lower_bound, N_over_L, B)
    # res_file.write('%f\t%0.12f\n' % (beta_lower_bound, likelihood_beta_lower))
    if likelihood_beta_lower > .0:
        while True:            
            likelihood_beta_upper = likelihood_derivative_beta(normalized_sitecount, beta_upper_bound, N_over_L, B)
            if likelihood_beta_upper < 0.:
                beta_lower_bound = beta_upper_bound - step_size
                break
            beta_upper_bound += step_size
            # res_file.write('%f\t%0.12f\n' % (beta_upper_bound, likelihood_beta_upper))            
            if beta_upper_bound > beta_max:  # prevents to go very very large numbers for Beta
                return beta_max
    else:  # once the derivative of likelihood for the lowest beta is already negative
        return beta_lower_bound

    # bisection method for refining the search area
    while (beta_upper_bound - beta_lower_bound) > 1.0e-9:
        mid_beta = (beta_upper_bound + beta_lower_bound) / 2.0
        # print beta_lower_bound, mid_beta, beta_upper_bound
        mid_likelihood = likelihood_derivative_beta(normalized_sitecount, mid_beta, N_over_L, B)
        # print beta_lower_bound, mid_beta, beta_upper_bound, mid_likelihood        
        if mid_likelihood > .0:
            beta_lower_bound = mid_beta
        else:
            beta_upper_bound = mid_beta
    # print mid_beta
    # res_file.write('beta:\t%0.12f\n' % mid_beta)
    # res_file.close()
    return mid_beta

# def fit_beta(siteFile, interm_dir, wmFile):
#     binding_regions_sitecount, total_sitecount, number_of_windows, total_length, M = sum_of_posteriors_foreground_regions(siteFile)         
#     beta_min, beta_max = 1.0e-12, 1.0e+3
#     step_size = 5.0
#     beta_lower_bound, beta_upper_bound = beta_min, beta_min+step_size
#     N_over_L = total_sitecount / float(total_length)
#     B = len(binding_regions_sitecount.keys())
#     tmp = [] # it's defined to speed up the analysis
#     # instead of calculating \sum_{r \in B} \frac{l_r}{beta*l_r + n_r}
#     # I calculate: \sum_{r \in B} \frac{1}{beta + n_r/l_r}
#     # where n_r/l_r is like site per nucleotide for a sequence
#     for region, sitecount in binding_regions_sitecount.items():
#         tmp.append(sitecount / float(number_of_windows[region]))
#     normalized_sitecount = np.array(tmp)
#     # updating the lower and upper bounds for the value of Beta
#     # It's where the two bounds (upper and lower) have opposite signs
#     likelihood_beta_lower = likelihood_derivative_beta(normalized_sitecount, beta_lower_bound, N_over_L, B)
#     if likelihood_beta_lower > .0:
#         while True:
#             likelihood_beta_upper = likelihood_derivative_beta(normalized_sitecount, beta_upper_bound, N_over_L, B)
#             if likelihood_beta_upper < 0.:
#                 beta_lower_bound = beta_upper_bound - step_size
#                 break
#             beta_upper_bound += step_size
#             if beta_upper_bound > beta_max:  # prevents to go very very large numbers for Beta
#                 return beta_max
#     else:  # once the derivative of likelihood for the lowest beta is already negative
#         return beta_lower_bound

#     print beta_lower_bound, beta_upper_bound
#     # bisection method for refining the search area
#     while (beta_upper_bound - beta_lower_bound) > 1.0e-9:
#         mid_beta = (beta_upper_bound + beta_lower_bound) / 2.0
#         # print beta_lower_bound, mid_beta, beta_upper_bound
#         mid_likelihood = likelihood_derivative_beta(normalized_sitecount, mid_beta, N_over_L, B)
#         print beta_lower_bound, mid_beta, beta_upper_bound, mid_likelihood        
#         if mid_likelihood > .0:
#             beta_lower_bound = mid_beta
#         else:
#             beta_upper_bound = mid_beta
#     print mid_beta
#     motifName = os.path.basename(wmFile)
#     # res_file = open(os.path.join(interm_dir, '%s.beta' % motifName), 'w')
#     # res_file.write('beta:\t%0.12f\n' % mid_beta)
#     # res_file.close()
#     return mid_beta
