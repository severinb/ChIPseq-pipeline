#!/usr/bin/env python

import sys
import os
from string import *
from pylab import *
import component_skeleton.main
import datetime

def normalizeHeights(infile, covdir, outfile, peakstatsfile, peakstatsoutfile):

    #make dictionary with peakstats
    statsDict = {}
    for line in open(peakstatsfile):
        t = line.strip().split()
        statsDict[t[0]] = map(float, t[1:])

    o = open(outfile, 'w')

    totnum = 0 #number of input peaks
    notfound = 0 #counts how many peaks don't have any background coverage

    for line in open(infile):

        # chr9    121791404       121791495       reg1000058.p4   93.617  +
        # chr17   29251036        29251147        reg1000126.p2   93.586  +
        # chr16   48788898        48788992        reg1000093.p2   93.419  +
        # chr2    52784183        52784281        reg1000145.p1   93.295  +

        totnum += 1

        t = line.split()

        ID = '_'.join(t[:-3])
        chrom = '_'.join(t[:-5])
        start = int(t[-5])
        end = int(t[-4])

        height = float(t[-2])
        h1 = height

        #peak has to be inside the coordinates of a background coverage file
        found = False
        m = int((float(start) + float(end))/2.0) #middle of peak
        for covfile in os.listdir(covdir):

            tCov = covfile.strip().split('_')
            chromCov = '_'.join(tCov[:-2])
            startCov = int(tCov[-2])
            endCov = int(tCov[-1])

            ml = int(m - startCov)
            mr = int(endCov - m)

            if chromCov == chrom:
                if abs(ml-mr) < 5: #actually it should be exactly the same, but maybe integer function does some rounding and it gets different...
                    found = covfile
                    break
                else:
                    continue
            else:
                continue


        if found:
            regcov = os.path.join(covdir, found)

            #add some sort of pseudo count of 1
            a = loadtxt(regcov, usecols=[5]) + 1

        else:
            a = [1,1]
            notfound += 1
            print '%s not found' %ID

        height /= mean(a)
        t[-2] = str(height)

        peakID = t[-3]  #e.g. reg10000232.p1
        statsDict[peakID][0] /= mean(a) #also normalize height in peakstats file

        o.write('\t'.join(t) + '\n')
        
    o.close()

    po = open(peakstatsoutfile, 'w')
    for key in statsDict:
        po.write('\t'.join([key] + map(str, statsDict[key])) + '\n')
    po.close()
    
    return totnum, notfound


def execute(cf):

    in_file = cf.get_input("in_file") #allpeaks (or peakscores)
    covdir = cf.get_input("coverage_dir")
    peakstats = cf.get_input("peakstats_infile")
    peakscores = cf.get_output("outfile")
    allpeaks = cf.get_output("allpeaks")
    log_file = cf.get_output("log_file")
    peakstatsout = cf.get_output("peakstats")

    height_revcum = cf.get_output("height_revcum")
    height_sigma_scatter = cf.get_output("height_sigma_scatter")
    height_hist = cf.get_output("height_hist")
    height_before_after_scatter = cf.get_output("height_before_after_scatter")
    #height_rmsd_scatter = cf.get_output("height_rmsd_scatter")

    topPeaks = cf.get_parameter("topPeaks", "int")

    T1 = datetime.datetime.now()

    normedfile = os.path.join(os.path.split(allpeaks)[0], 'normedfile')
    totnum, notfound = normalizeHeights(in_file, covdir, normedfile, peakstats, peakstatsout)

    #sort file, as it may has different order after normalizing
    os.system("sort -k5gnr %s > %s" %(normedfile, allpeaks))

    #write toppeaks
    p = open(peakscores, 'w')
    for i, line in enumerate(open(allpeaks)):
        if i >= topPeaks:
            t = line.strip().split()
            lastheight = float(t[-2])
            break
        else:
            p.write(line)

    p.close()

    #plot stuff
    a = loadtxt(allpeaks, usecols=[1,2,4])
    sigmas = (a.T[1]-a.T[0])/2.
    heights = a.T[2]+0.01 #add pseudocount
    #quals = a.T[3] (maybe get quals from peakstats file of selectPeaksMM, but I don't know whether this is interesting)

    figure()
    plot(log10(heights), arange(1,len(heights)+1,1), '.')
    lh = log10(lastheight)
    plot([lh, lh], [0, len(heights)], label='lowest given out peak: %s' %lastheight)
    xlabel('log10(heights)')
    ylabel('number of peaks with up to peak height')
    legend()
    savefig(height_revcum)
    close()

    figure()
    hist(heights, 50)
    title('histogram of peak heights')
    savefig(height_hist)
    close()

    figure()
    plot(sigmas, heights, '.')
    xlabel('peak sigma')
    ylabel('peak height')
    savefig(height_sigma_scatter)
    close()


    #produce peakheight before versus after normalization scatter
    heightbeforeDict = {}
    for line in open(peakstats):
        t = line.strip().split()
        heightbeforeDict[t[0]] = float(t[1])

    before = []
    after = []

    for line in open(peakstatsout):
        t = line.strip().split()
        ID = t[0]
        hafter = float(t[1])
        after.append(hafter)
        before.append(heightbeforeDict[ID])

    figure()
    plot(array(before), array(after), '.')
    xlabel('peak height before normalization')
    ylabel('peak height after normalization')
    savefig(height_before_after_scatter)
    close()


    # figure()
    # plot(heights, quals, '.')
    # xlabel('peak height')
    # ylabel('peak quality')
    # savefig(height_rmsd_scatter)
    # close()


    T2 = datetime.datetime.now()

    text = '\n'.join(['Running time: %s' %str(T2-T1),
                      'Number of input peaks: %i' %totnum,
                      'Number of peaks that had no background coverage: %i' %notfound])
    lf = open(log_file, 'w')
    lf.write(text)
    lf.close()


    return 0


component_skeleton.main.main(execute)
