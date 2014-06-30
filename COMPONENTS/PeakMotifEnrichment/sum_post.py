#!/usr/bin/env python

import sys, os
from pylab import *

def main():

    sites = sys.argv[1]

    # 25-36 + 0.012856 /import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/PipeLineSource/TESTRUN/NRF1_Z3/OUTPUT/NRF1_FgBg-phylogibbs2/Logo2 hg19_chr6_apd_hap1_137830_138079_MACS_peak_13639_52.02
    # 140-151 + 0.056816 /import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/PipeLineSource/TESTRUN/NRF1_Z3/OUTPUT/NRF1_FgBg-phylogibbs2/Logo2 hg19_chr6_apd_hap1_137830_138079_MACS_peak_13639_52.02
    # 155-166 - 0.043819 /import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/PipeLineSource/TESTRUN/NRF1_Z3/OUTPUT/NRF1_FgBg-phylogibbs2/Logo2 hg19_chr6_apd_hap1_137830_138079_MACS_peak_13639_52.02


    d = {}
    for line in open(sites):
        t = line.strip().split()
        try:
            d[t[4]].append(float(t[2]))
        except KeyError:
            d[t[4]] = []
            d[t[4]].append(float(t[2]))            

    p = []
    s = []

    for k in d:
        t = k.split('_')
        sco = float(t[-2])
        p.append(sum(d[k]))
        s.append(sco)

    plot(s, p, '.')
    savefig('macs_sco_post.pdf')

    parr = array(p)
    sarr = array(s)

    print len(parr[where(parr >= 0.2)])
    print len(parr)


    # compute fraction of explained variance of each.                                                                                                                                                                                        
    covmat = cov(sarr, parr)
    fov1 = (1 + sqrt( (covmat[0][1]**2)/(covmat[0][0] * covmat[1][1]) ))/2
    print 'FOV: ', fov1


if __name__ == '__main__':
    main()
