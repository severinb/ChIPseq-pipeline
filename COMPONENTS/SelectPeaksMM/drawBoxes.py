#! /usr/bin/env python

import sys
from pylab import *

if __name__ == '__main__':

    infile = '/import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/COMPONENTS/COMPONENTS.v5/SelectPeaksMM/STATS'

    meanFraglens = []
    boxes = []
    TFs = []

    for line in open(infile):
        if line.startswith('#'):
            continue
        else:
            t = line.strip().split()

            TFs.append(t[0])
            meanFraglens.append(int( mean( map(float, t[4:-1]) )))

            a = loadtxt(t[-1], usecols=[1,2])

            s = (a.T[1] - a.T[0])/2.

            boxes.append(s)


    print TFs
    print meanFraglens

    totmean = mean(map(mean, boxes))

    boxplot(boxes, positions=meanFraglens)
    xticks(rotation=30)
    xlabel('mean Fraglen')
    ylabel('sigmas')
    #xticks(meanFraglens, TFs)
    plot([min(meanFraglens), max(meanFraglens)], [totmean, totmean], 'r', label='total sigma mean')
    legend()
    savefig('boxplotsigma.pdf')

