#!/usr/bin/env python

import os, re
from string import *
import component_skeleton.main
from pylab import *


def execute(cf):


    outplot = cf.get_output("correlation_plot")

    peakstats_string = cf.get_parameter("peakstats_string", "string")

    statsfiles = peakstats_string.split()

    motif_num = len(statsfiles)

    # #chrom  start   end     peakID  height  quality summed_posterior
    # chr11   45355348        45355469        reg1006710.p1   3.616   0.033969388741  0.86503
    # chrX    153829376       153829535       reg1003967.p4   3.486   0.0470921572293 0.038494
    # chr6    14984749        14984852        reg1008286.p2   3.336   0.0389569986723 0.587459
    # chr6    10441258        10441370        reg1000562.p1   7.151   0.0184399263633 0.014908


    # a list of dictionaries for each given motif. key: peakID, value: peak posterior
    WM_dicts = []

    for tf in statsfiles:
        WM1_d = {}
        for line in open(tf):
            if line.startswith('#'):
                continue
            t = line.strip().split()

            WM1_d[t[3]] = float(t[6])

        WM_dicts.append(WM1_d)

    peakids = [k for d in WM_dicts for k in d.keys()]

    unique_peakids = unique(peakids)

    # make a matrix of all posteriors for all motifs. First list in lists contains posteriors of WM_1 and so on. I make sure here that each row of the matrix is one peak.
    # So rows are peaks and columns are motifs.

    mat = []

    for wm in WM_dicts:

        mat.append([])

        for k in unique_peakids:
            try:
                post = wm[k]
            except KeyError:
                psot = 0.0

            mat[-1].append(post)


    print cov(array(mat))
    print corrcoef(array(mat))

    figure()

    cmat = corrcoef(array(mat))
    try:
        m_shape = cmat.shape[0]
        mask = tri(m_shape, k=-1).T
        C = ma.array(cmat, mask=mask)
        imshow(C, interpolation="nearest")
        colorbar()
        xticks(range(motif_num),['WM_%i' %i for i in range(1,motif_num+1)], rotation=30)
        yticks(range(motif_num),['WM_%i' %i for i in range(1,motif_num+1)])
        title('Correlation of Motif Predictions')

    except AttributeError: # exception for one dimensional C, because we only have one matrix
        pass

    savefig(outplot)


    return 0


component_skeleton.main.main(execute)


