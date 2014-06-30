#! /usr/bin/env python

import sys
from pylab import *

if __name__ == '__main__':

    infile = '/import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/COMPONENTS/COMPONENTS.v5/SelectPeaksMM/STATS1gauss'

    meanFraglens = []
    minFraglens = []
    maxFraglens = []
    boxes = []
    TFs = []

    for line in open(infile):
        if line.startswith('#'):
            continue
        else:
            t = line.strip().split()

            TFs.append(t[0])
            meanFraglens.append(float( mean( map(float, t[4:-1]) )))
            minFraglens.append(float( min( map(float, t[4:-1]))) )
            maxFraglens.append(float( max( map(float, t[4:-1]))) )


            a = loadtxt(t[-1], usecols=[1,2])

            s = (a.T[1] - a.T[0])/2.0

            boxes.append(s)

    print TFs
    print meanFraglens
    print minFraglens
    print maxFraglens

    totmean = mean(map(mean, boxes))

    fraglens = array(maxFraglens) -10

    print fraglens

    boxplot(boxes, positions=fraglens)
    xticks(rotation=30)
    xlabel('max Fraglen')
    ylabel('sigmas')
    #xticks(meanFraglens, TFs)
    plot([min(fraglens), max(fraglens)], [totmean, totmean], 'r', label='total sigma mean')
    #plot([min(fraglens), max(fraglens)], [100, 100])
    #plot([min(fraglens), max(fraglens)], [30, 30])
    plot(fraglens, (1.0/(2*sqrt(2*log(2))))*(array(fraglens)), label='1/(2*sqrt(2*log(2)))*fraglen')
    plot(fraglens, 0.3*array(fraglens)+40, label='0.3*fraglen+60')
    plot(fraglens, 0.3*array(fraglens)+20, label='0.3*fraglen')
    legend()
    savefig('boxplotsigma1g.pdf')

