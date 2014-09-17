import re, os
import numpy as np

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


def calculate_enrichment_scores(siteFile, beta, length, res_filename):    
    sites, N = sum_of_posteriors_foreground_regions(siteFile)
    L = np.sum([l for l in length.values()])
    averageLength = np.mean([l for l in length.values()])
    M = len(length.keys())
    fg_regions = M / 11  # 10 times more than fg sequences there's shuffled ones
    denumerator = np.log( N + L*beta ) - np.log(M)
    enrichmentScores = np.repeat(np.log(beta*averageLength) - denumerator, fg_regions)
    try:
        index = 0
        # with open(res_filename, 'w') as outfile:        
        for region, sitecount in sites.items():
                # outfile.write('\t'.join([
                #         region,
                #         '%0.10f\n' % ( np.log( sitecount + length[region]*beta )  - denumerator)]))
            enrichmentScores[index] = (np.log( sitecount + length[region]*beta ) - denumerator)
            index += 1
    except:
        return {}
    return {'mean':np.exp(np.mean(enrichmentScores)),
            'std':np.exp(np.std(enrichmentScores)),
            'LL_ratio':np.sum(enrichmentScores)} 
