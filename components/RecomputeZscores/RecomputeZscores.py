#!/usr/bin/env python

import sys, os, re
from pylab import *
import datetime
import component_skeleton.main
from string import *

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

    peakfile = cf.get_input('peakfile') #from SelectPeaksMM
    peakstats =  cf.get_input('peakstats') #from SelectPeaksMM
    readcount_file = cf.get_input('binned_reads') #from BinReads
    fit_log = cf.get_input('fit_log')
    zco_file = cf.get_input('Z_cutoff_file')

    outfile = cf.get_output('allpeaks') #height gets replaced by Z-value. Height can still be found in peakstats file of SelectPeaksMM
    tops = cf.get_output('outfile')
    outstats = cf.get_output('peakstats') #same file as peakstats infile but with additional column containing Z-scores
    log_file = cf.get_output('log_file')
    scatterplot = cf.get_output('height_z_scatter')
    zhist = cf.get_output('z_hist')
    revcumplot = cf.get_output('z_revcum')

    winsize = cf.get_parameter('FGwinsize', 'int')
    Z_co = cf.get_parameter('Z_cutoff', 'float')
    FG_read_string = cf.get_parameter('read_files', 'string') #for fragment length
    toppeaks = cf.get_parameter('topPeaks', 'int')

    if Z_co <= 0:
        Z_co = float(open(zco_file).readline().strip())

    #Also take total counts for foreground and background from PeakCall_newNoiseModel. Now I compute it myself, although it's not exactly the same numbers

    # peakfile:
    # chr19   39894602        39894725        reg1000001.p8   1715.069        +
    # chr3    47517431        47517540        reg1000011.p1   1271.075        +
    # chr15   63796704        63796799        reg1000019.p4   1178.874        +

    # peakstats:
    # reg1050315.p1   1.110   0.0329590802201 80.578  41.250 3.234
    # reg1050315.p2   4.925   0.0958940756527 441.062 41.250 3.234
    # reg1021921.p1   2.715   0.0501350096756 277.059 105.453 2.345

    # readcount_file
    # #chr    start   stop    middle  fg_0    fg_1    bg_0    bg_1
    # chr1    9000    9500    9250    0       0       2.682076        2.454219
    # chr1    9250    9750    9500    0       0       4.281291        4.251841
    # chr1    9500    10000   9750    0       0       4.340815        4.655877

    # fit_log:
    # Fitted parameters:
    #         -sigma: 0.230846850688983
    #         -rho: 0.0163927715194372
    #         -mu: -0.106772262679966

    T1 = datetime.datetime.now()

    fraglen = findFraglen(FG_read_string)

    # read mu and sigma from noise distribution
    for line in open(fit_log):
        t = line.strip().split()
        try:
            if t[0] == '-sigma:':
                noise_sigma = float(t[1])
            elif t[0] == '-mu:':
                noise_mu = float(t[1])
        except IndexError:
            continue

    if not noise_sigma or not noise_mu:
        print 'Error: Could not find sigma and rho of noise distribution\n'
        sys.exit(1)

    # make a dictionary with peakID as key and height, rmsd, mu, sigma and uniform height as value, from peakstats.
    stats_dict = {}
    for l in open(peakstats):
        t = l.strip().split()
        stats_dict[t[0]] = [float(t[1]), float(t[2]), float(t[3]), float(t[4]), float(t[5])]


    # load background counts into memory
    f = open(readcount_file)
    header = f.readline()
    f.close()

    fgnum = 0
    bgnum = 0
    for i in header.strip().split():
        if i.startswith('fg'):
            fgnum += 1
        elif i.startswith('bg'):
            bgnum += 1

    chroms = loadtxt(readcount_file, skiprows=1, usecols=[0], dtype=str)

    # get total counts
    N = len(chroms) #number of windows

    idx = [3+i+1 for i in arange(fgnum)]
    a = loadtxt(readcount_file, skiprows=1, usecols=idx)
    totfg = sum(a)

    if not bgnum == 0:
        idx = [3] + [3+fgnum+i+1 for i in arange(bgnum)] 
        coords_counts = loadtxt(readcount_file, skiprows=1, usecols=idx)
    else:
        coords_counts = array( zip( loadtxt(readcount_file, skiprows=1, usecols=3), ones(len(chroms)) * totfg/N))
        print 'No background sample given. Uniform background of %s assumed (%s/%s).' %(totfg/N, totfg, N)

    #make a dictionary where each chromosome is a key, to make it faster afterwards:
    chrom_dict = {}
    for c in unique(chroms):
        chrom_dict[c] = coords_counts[where(chroms == c)]

    totbg = sum(coords_counts[:, 1:]) #middle coordinates in first column

    pseudocount = 0.5 #defined like this in PeakCall_newNoiseModel

    totdenom = totfg + totbg + N*pseudocount
    logcomp = log((totbg + N*pseudocount) / (totfg + N*pseudocount))


    T2 = datetime.datetime.now()

    # now compute Z-score for every peak:
    o1 = open(outstats, 'w')

    allstats = os.path.join(os.path.split(outstats)[0], 'allpeaks_allstats')
    o2 = open(allstats, 'w')
    o2.write('#chrom\tstart\tstop\tpeakID\tZscore\tstrand\theight\trmsd\tmu\tsigma\tbackheight\n')

    tmpfile = os.path.join(os.path.split(outstats)[0], 'tmpfile')
    o = open(tmpfile, 'w')
    goodpeaks = 0
    allpeaks = 0

    for l in open(peakfile):
        t = l.strip().split()
        chrom = t[0]
        start = int(t[1])
        stop = int(t[2])
        peakID = t[3]
        strand = t[5]

        height, rmsd, mu, sigma, backheight = stats_dict[peakID]

        # salvage read_counts
        readcount_peak = height * sqrt(2*3.141592653589793) * sigma / fraglen
        readcount_uniform = backheight * winsize / fraglen
        readcount_fg = readcount_peak + readcount_uniform

        candidates = chrom_dict[chrom]
        nearest_index = argmin(abs(candidates.T[0] - ((start+stop)/2)))
        nearest_win = candidates[nearest_index]

        # Now i sum up the background counts, BUT I took the average of coverage for the coverage profiles. So the peak height is average of replicates. This is of course wrong. 
        # So check how this looks like and then take also sum for the coverage profiles!
        readcount_fg *= fgnum #reverse mean coverage (that is height)
        readcount_bg = sum(nearest_win[1:])

        fg_count_binned = totfg * ((readcount_fg + readcount_bg + pseudocount) / totdenom)
        bg_count_binned = totbg * ((readcount_fg + readcount_bg + pseudocount) / totdenom)

        denom = sqrt(2*(noise_sigma**2) + 1.0/fg_count_binned + 1.0/bg_count_binned) #sqrt because we want standard deviation not variance

        numerator = log((readcount_fg + pseudocount) / (readcount_bg + pseudocount)) + logcomp - noise_mu

        Z = numerator/denom
        o2.write('%s\t%s\t%s\t%s\t%.3f\t%s\t%s\t%s\t%s\t%s\t%s\n' %(chrom, start, stop, peakID, Z, strand, height, rmsd, mu, sigma, backheight))
        allpeaks += 1

        if Z >= Z_co:
            o.write('%s\t%s\t%s\t%s\t%.3f\t%s\n' %(chrom, start, stop, peakID, Z, strand))
            o1.write('%s\t%s\t%s\t%s\t%s\t%s\t%.3f\n' %(peakID, height, rmsd, mu, sigma, backheight, Z))
            goodpeaks += 1

    o1.close()
    o.close()
    o2.close()


    # sort tmpfile and write toppeaks to extra file
    os.system('sort -k5gnr %s > %s' %(tmpfile, outfile))
    os.system('rm %s' %tmpfile)

    o = open(tops, 'w')
    i = 1
    for l in open(outfile):
        if i <= toppeaks:
            t = l.strip().split()
            lastz = float(t[4])
            o.write(l)
            i += 1
        else:
            break

    o.close()

    # plot scatter of height versus Z score
    a = loadtxt(allstats, usecols=[6,4])
    figure()
    plot(a[:,0], a[:,1], '.', rasterized=True)
    xlabel('peak height')
    ylabel('peak Z-score')
    savefig(scatterplot)
    close()


    h = hist(a[:,1],300, histtype='step')
    close()
    hx = [mean([h[1][i], h[1][i+1]]) for i in arange(len(h[0]))]
    figure()
    plot(hx, log(h[0]))
    plot([lastz, lastz],[0, log(max(h[0]))], label='lowest top Z-score: %s' %lastz)
    plot([Z_co, Z_co],[0, log(max(h[0]))], label='Z-score cut-off: %s' %Z_co)
    title('Histogram of Peak Z-scores (log-counts)')
    xlabel('peak Z-score')
    ylabel('log-number of peaks')
    legend()
    savefig(zhist)
    close()

    figure()
    plot(sorted(a.T[1], reverse=1), log(arange(1, len(a.T[1])+1, 1)))
    plot([lastz, lastz],[0, log(len(a.T[0])-1)], label='lowest top Z-score: %s' %lastz)
    plot([Z_co, Z_co],[0, log(len(a.T[0])-1)], label='Z-score cut-off: %s' %Z_co)
    xlabel('peak Z-score')
    ylabel('log-number of peaks with at least Z-score')
    title('Reverse Cumulative Distribution of Peak Z-scores (log-counts)')
    legend()
    savefig(revcumplot)
    close()

    T3 = datetime.datetime.now()

    # write log
    l = open(log_file, 'w')
    timetext = '\n'.join(['Running Time:',
                          '\t-Read in data: %s' %str(T2-T1),
                          '\t-Recalculating Z-score: %s' %(T3-T2),
                          '\t-Overall: %s' %(T3-T1)])
    text = '\n'.join(['Total number of peaks: %i' %allpeaks,
                      'Number of peaks above Z-score of %s: %i\n' %(Z_co, goodpeaks)])

    l.write(text)
    #l.write(timetext)
    l.close()


    return 0


component_skeleton.main.main(execute)
