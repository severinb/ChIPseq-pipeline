#!/usr/bin/env python

import os, re
from string import *
import component_skeleton.main
from pylab import *

def findFraglen(readfiles):

    fraglens = []
    rlist = readfiles.split()
    for rf in rlist:
        ##get Fragment length from the log of FragmentLength component
        fraglenRoot = os.path.join(os.path.split(rf)[0], '..', 'intermediate')
        interdir = os.listdir(fraglenRoot)
        for f in interdir:
            m = re.match('\S*.res$', os.path.join(fraglenRoot, f))
            if m:
                fraglenF = m.group()
                break
        fraglen = int(float(open(fraglenF).read().strip().split()[-1]))
        fraglens.append(fraglen)
        print '\nFound: %s with fragment length: %s\n' %(rf,fraglen)
            
    meanFraglen = mean(fraglens)
    maxFraglen= max(fraglens)
    print '\nFragment Length: %s' %fraglens

    return meanFraglen


def execute(cf):

    peaks = cf.get_input("peaks")

    # peaks with annotations
    # #chrom  start   stop    peakID  Z-score strand  nearest_upstream_promoter_members       distance        2nearest_upstream_promoter_members      distance        3nearest_upstream_promoter_members      distance        nearest_downstream_promoter_members     distance        2nearest_downstream_promoter_members    distance        3nearest_downstream_promoter_members    genes[member|symbol|gid|name]
    # chr11   6640614 6640717 reg1000000.p2   12.526  +
    # chr17   80246010        80246121        reg1000003.p1   11.668  +
    # chr7    915954  916079  reg1000001.p2   11.597  +
    # chr19   3990331 3990439 reg1000006.p2   11.345  +

    outfile = cf.get_output("peaks_with_sites")

    knowntfbs_string = cf.get_parameter("knownTFBSs_string", "string")
    denovotfbs_string = cf.get_parameter("denovoTFBSs_string", "string")
    knownpeakposts_string = cf.get_parameter("known_peakposts_string", "string")
    denovopeakposts_string = cf.get_parameter("denovo_peakposts_string", "string")
    FGfiles_string = cf.get_parameter("FGfiles_string", "string")
    minposterior = cf.get_parameter("minposterior", "float")


    # *tfbs_string and *peakposts_strin:
    # WM_path filepath; WM_path filepath; ...
    known_names = []
    denovo_names = []
    knowntfbsfiles = []
    denovotfbsfiles = []
    knownpeakfiles = []
    denovopeakfiles = []

    for pair in knowntfbs_string.split(';'):
        t = pair.split()
        known_names.append(os.path.split(t[0])[1])
        knowntfbsfiles.append(t[1])

    for pair in knownpeakposts_string.split(';'):
        t = pair.split()
        knownpeakfiles.append(t[1])

    for pair in denovotfbs_string.split(';'):
        t = pair.split()
        denovo_names.append(os.path.split(t[0])[1])
        denovotfbsfiles.append(t[1])

    for pair in denovopeakposts_string.split(';'):
        t = pair.split()
        denovopeakfiles.append(t[1])



    # TFBS file
    # chr15   72410267        72410277        reg1002725.p1   320.179 0.012779        2.0
    # chr15   72410416        72410426        reg1002725.p1   171.179 0.030356        4.5
    # chr15   72410436        72410446        reg1002725.p1   151.179 0.215111        5.5
    # chr15   72410437        72410447        reg1002725.p1   150.179 0.020911        5.5


    # peak posts file:

    fraglen = findFraglen(FGfiles_string)

    # this list contains a dictionary for each given motif. key: peakID, value: list of sites
    knownWM_dicts = []

    for tf in knowntfbsfiles:
        WM1_d = {}
        for line in open(tf):
            if line.startswith('#'):
                continue
            t = line.strip().split()
            dist = int(float(t[4]))

            if float(t[5]) >= minposterior:
                if abs(dist) <= fraglen:
                    try:
                        WM1_d[t[3]].append([dist, float(t[5])])
                    except KeyError:
                        WM1_d[t[3]] = [[dist, float(t[5])]]

        knownWM_dicts.append(WM1_d)


    denovoWM_dicts = []

    for tf in denovotfbsfiles:
        WM1_d = {}
        for line in open(tf):
            if line.startswith('#'):
                continue
            t = line.strip().split()
            dist = int(float(t[4]))

            if float(t[5]) >= minposterior:
                if abs(dist) <= fraglen:
                    try:
                        WM1_d[t[3]].append([dist, float(t[5])])
                    except KeyError:
                        WM1_d[t[3]] = [[dist, float(t[5])]]

        denovoWM_dicts.append(WM1_d)

    # list containing dictionaries. key: peakID, value: summed_posterior
    wmnames = []

    knownWM_dicts_peaks = []

    for i, tf in enumerate(knownpeakfiles):
        WM1_d = {}
        for line in open(tf):
            if line.startswith('#'):
                continue
            t = line.strip().split()
            post = t[6]

            WM1_d[t[3]] = post

        knownWM_dicts_peaks.append(WM1_d)
        wmnames.append('%s_site_number' %(known_names[i]))



    denovoWM_dicts_peaks = []

    for i, tf in enumerate(denovopeakfiles):
        WM1_d = {}
        for line in open(tf):
            if line.startswith('#'):
                continue
            t = line.strip().split()
            post = t[6]

            WM1_d[t[3]] = post

        denovoWM_dicts_peaks.append(WM1_d)
        wmnames.append('%s_site_number' %(denovo_names[i]))


    o = open(outfile, 'w')

    # chr11   6640616 6640715 reg1000000.p2   12.456  +       WM_1:-4:0.676   WM_1:6:0.741    WM_2:-69:0.490  WM_2:-102:0.864 WM_2:111:0.380

    # o.write("#chromosome\tstart\tstop\tpeakID\tZscore\tstrand\tnearest_upstream_promoter_members\tdistance\t2nearest_upstream_promoter_members\tdistance\t3nearest_upstream_promoter_members\tdistance\tnearest_downstream_promoter_members\tdistance\t2nearest_downstream_promoter_members\tdistance\t3nearest_downstream_promoter_members\tgenes[member|symbol|gid|name]\t")

    # o.write('\t'.join(wmnames))

    # o.write("\tbinding_sites (WM_name:relative distance to peak center:posterior binding probability)\n")


    for line in open(peaks):
        if line.startswith('#'): #write header line
            o.write(line.strip())
            o.write('\t')
            o.write('\t'.join(wmnames))
            o.write("\tbinding_sites (WM_name:relative distance to peak center:posterior binding probability)\n")
            continue

        t = line.strip().split()

        peakID = t[3]

        posts = []

        for j, wmdict in enumerate(knownWM_dicts_peaks):
            try:
                post = wmdict[peakID]
                posts.append(post)
            except KeyError:
                posts.append('0')

        for j, wmdict in enumerate(denovoWM_dicts_peaks):
            try:
                post = wmdict[peakID]
                posts.append(post)
            except KeyError:
                posts.append('0')

        sites = []

        for j, wmdict in enumerate(knownWM_dicts):
            try:
                sitelist = wmdict[peakID]
                for site in sitelist:
                    sites.append('%s:%i:%.3f' %(known_names[j], site[0], site[1]))
            except KeyError:
                pass

        for j, wmdict in enumerate(denovoWM_dicts):
            try:
                sitelist = wmdict[peakID]
                for site in sitelist:
                    sites.append('%s:%i:%.3f' %(denovo_names[j], site[0], site[1]))
            except KeyError:
                pass

        sorted_sites = sorted(sites, key = lambda k: abs(int(k.split(':')[1])))

        o.write(line.strip() + '\t' + '\t'.join(posts) + '\t' + '\t'.join(sorted_sites) + '\n')

    o.close()

 
    return 0


component_skeleton.main.main(execute)


