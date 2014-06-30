#!/usr/bin/env python

import sys, os
import component_skeleton.main

def execute(cf):

    regions = cf.get_input("regions") #alignments or sequences (containing peak ids)
    peakstats = cf.get_input("peakstats") #peakstats file of GetTrueRegion. In this file all the stats of a peak are stored: heights, rmsd and posts
    post_co = cf.get_parameter("posterior_cutoff", "float") #everything smaller is taken
    outfile = cf.get_output("outfile")

    #chrom  start   end     peakID  height  quality summed_posterior        summed_weighted_posterior
    peakids = []
    i = 1
    for line in open(peakstats):
        if i == 1:
            i += 1
            continue
        t = line.split()
        if float(t[6]) <= post_co:
            peakids.append(t[3])

    #>>hg19_chr8_37595265_37595360_reg1000661.p2_12.7821180556_+
    o = open(outfile, 'w')
    writeit = False
    for line in open(regions):
        if line.startswith('>>'):
            t = line.split('_')
            peakid = t[-3]
            if peakid in peakids:
                writeit = True
                o.write(line)
            else:
                writeit = False
        else:
            if writeit:
                o.write(line)
            else:
                continue

    return 0

component_skeleton.main.main(execute)
