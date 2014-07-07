#!/usr/bin/env python

import sys, os, re
from pylab import *
from string import *

def get_peak_posteriors(TFBS_file, peaks_file):
    """This funtion returns a dictionary with the peak (as from peakrefinement) as key and the sum of the posteriors of all TFBS inside the peak as value. The input is a <sites> file from motevo.
    """

    posteriors_dict = {}

    # initialize posterior_dict by the peak ids of all peaks where the WM was ran on. If I don't do this, no one will know whether a peak has posterior 0.
    # peaks file: >>hg19_chr12_123237285_123237405_reg1000001.p1_25.235_+

    for line in open(peaks_file):
        if line.startswith('>>'):
            peakid = line.lstrip('>>').strip()
            posteriors_dict[peakid] = 0.0

    # TFBS file: 1-20 + 0.000000 ZNF143.p2 hg19_chr12_123237285_123237405_reg1000001.p1_25.235_+
    f = open(TFBS_file, 'r')

    emptyFile = True
    while True:
        line = f.readline()
        if line:
            emptyFile = False
            TFBS = line.split()
            post = float(TFBS[2])
            peak = TFBS[4].strip()
            try:
                posteriors_dict[peak] += post
            except KeyError:
                posteriors_dict[peak] = post
        else:
            break
    f.close()

    if emptyFile:
        posteriors_dict['somekey'] = 0.0


    return posteriors_dict


def TP_P_FN_N(trueSites_dict, bgSites_dict):
    """This funtion calculates true positives TN, false positives P, true negatives TN and false negatives N for different cut offs for posteriors.
    """

    max_post = max(trueSites_dict.values() + bgSites_dict.values()) #*0.6
    #max_post = 1
    senslist = []
    ppvlist = []
    step = 0.0001
  
    sens = lambda TP, FN: (float(TP)/(TP + FN))
    ppv = lambda  TP, FP: (float(TP)/(TP + FP))

    for co in arange(0, max_post, step):
        TP = len([x for x in trueSites_dict.values() if float(x) >= float(co)])
        FN = len(trueSites_dict.values()) - TP 
        FP = len([x for x in bgSites_dict.values() if float(x) >= float(co)]) 
        TN = len(bgSites_dict.values()) - FP 

        #add pseudo counts proportionally to the probabilty of each category (10 times more false peaks!)
        #better without pseudo counts, otherwise there can be strange behaviour at the left top corner of the curve...
        # TP += 0.1
        # FN += 1.0
        # FP += 0.1
        # TN += 1.0

        senslist.append(sens(TP,FN))
        ppvlist.append(ppv(TP,FP))

        #print '\t'.join(['co %s :' , 'tp %s', 'p %s', 'fn %s', 'n %s' ,'tp+fn %s', 'p+n %s', 'sens %s','ppv %s']) %(co, TP, FP, FN, TN, TP+FN, FP+TN, sens(TP,FN), ppv(TP,FP))

    return senslist, ppvlist

def plot_ROC(nonre_senslist, nonre_ppvlist, re_senslist, re_ppvlist, outpdf):
    """This function plots sensitivity on x-axis and ppv on y-axis for both refined and non refined predictions.
    """
    """
    plt.plot(nonre_senslist, nonre_ppvlist, 'g', label='non refined Sites')
    plt.plot(re_senslist, re_ppvlist, 'k', label='refined Sites')
    plt.xlabel('sensitivity')
    plt.ylabel('ppv')
    plt.legend(loc='lower left')
    plt.savefig(outpdf)
    plt.close()
    """

    #print nonre_senslist, nonre_ppvlist
    figure()
    plot(nonre_senslist, nonre_ppvlist, 'g', label='non refined WM')
    plot(re_senslist, re_ppvlist, 'k', label='refined WM')
    xlabel('sensitivity')
    ylabel('ppv')
    ylim([0,1.1])
    legend(loc='lower left')
    savefig(outpdf)
    
    
def main():
    """Function that plots a ROC curve (sensitivity, positive predictive value) for TFBS predictions from motevo with refined and nonrefined WM.
       Input are four site files from motevo. Each one on true peaks and background peaks for each WM.
    """

    if len(sys.argv) != 6:
        print '\nUsage: ./ROC.py nonref_sites_true(on odd) ref_sites_true(on odd) nonref_sites_bg ref_sites_bg outdir\n'
        sys.exit(0)

    non_sites_true = sys.argv[1]
    ref_sites_true = sys.argv[2]
    non_sites_bg = sys.argv[3]
    ref_sites_bg = sys.argv[4]
    outdir = sys.argv[5]

    nonref_trueSitesPost_dict = get_peak_posteriors(non_sites_true)
    nonref_bgSitesPost_dict = get_peak_posteriors(non_sites_bg)
    ref_trueSitesPost_dict = get_peak_posteriors(ref_sites_true)
    ref_bgSitesPost_dict = get_peak_posteriors(ref_sites_bg)
    
    nonref_senslist,nonref_ppvlist = TP_P_FN_N(nonref_trueSitesPost_dict, nonref_bgSitesPost_dict)
    ref_senslist, ref_ppvlist = TP_P_FN_N(ref_trueSitesPost_dict, ref_bgSitesPost_dict)

    plot_ROC(nonref_senslist, nonref_ppvlist, ref_senslist, ref_ppvlist, outdir)
    print 'Non Refined'
    print nonref_senslist
    print nonref_ppvlist
    print 'Refined'
    print ref_senslist
    print ref_ppvlist

    """
    Maybe just take bg and true sites from one wm as input and return senslist and ppvlist so that several wm can be tested. (maybe good for testing different backgrounds in same plot or so)
    """


if __name__ == '__main__':
    main()
