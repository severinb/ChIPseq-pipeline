#!/usr/bin/env python
import component_skeleton.main
import subprocess
import os
from string import *
import datetime
from pylab import *

def get_zco(x, y, outplot):

    figure()

    #to make a stable analysis, sample the data linearly spaced in x
    size = 100
    xn = linspace(min(x),max(x),size)
    yn = []

    for xi in xn:
        yidx = argmin(abs(x-xi))
        yn.append(y[yidx])

    x = array(xn)
    y = array(yn)

    plot(x,y)
    
    # 1. find linear fit at the beginning of the curve. Take first 100 data points
    r1 = size/20
    X = vstack([x[:r1], np.ones(len(x[:r1]))]).T
    params = lstsq(X, y[:r1])

    fx1 = x*params[0][0]+params[0][1]

    # more save to just take a line with slope 0
    fx1 = x*0.0+max(y)
    plot(x,fx1)

    # 2. compute deviation from linear fit and actual y to find the point where actual y goes down
    dy = fx1 - y
    idx1 = where(dy>1)[0][0]
    plot([x[idx1], x[idx1]], [0,max(y)])

    # 3. fit a second line from this index on. Try to find the steepest line in some region
    slopes = []
    inters = []
    minslope = 10
    for i in arange(2,size/5, 1):
        r2 = i

        X = vstack([x[idx1:idx1+r2], np.ones(len(x[idx1:idx1+r2]))]).T
        params = lstsq(X, y[idx1:idx1+r2])

        fx2 = x*params[0][0]+params[0][1]

        slopes.append(params[0][0])
        inters.append(params[0][1])


    slope = min(slopes)
    inte = inters[argmin(slopes)]
    upperb = arange(2,size/5, 1)[argmin(slopes)]

    fx2 = x*slope+inte

    plot([x[idx1+upperb], x[idx1+upperb]], [0,max(y)])
    plot(x[idx1:idx1+upperb],fx2[idx1:idx1+upperb])

    # 3. again compute deviation
    dy2 = fx2 - y
    idx2 = where(dy2<=-0.5)[0][0]

    co = x[idx2]

    plot([co,co],[0,12], label='Z cut-off %i' %int(co))
    xlabel('predicted WM-score')
    ylabel('Windows with up to predicted WM-score')

    legend()
    savefig(outplot)


    return int(co)


def merge_peaks(infile, outfile, cutoff):
    """This function filters peaks by a cut off and merges overlapping peaks. 
       Stores highest WM-score for each window.
       chrom start end strand(+)	 
    """
    
    f = open(infile, 'r')
    peaks = {}
    i = 1

    #skip first header row
    f.readline()

    #initializing by finding the first read above cutoff
    init = True
    while init:
        line = f.readline().split()
        zval = float(line[-1])
        if zval >= cutoff:
            peaks[i] = [line[0], int(line[1]), int(line[2]), zval]  #[chrom, start, end, zval]
            init = False
            break
        else:
            continue
    
    eof = False #tag for end of file
    while not eof:
        line = f.readline().split()
        if not line:
            eof = True
            break
        else:
            if float(line[-1]) >= cutoff:
                if line[0] == peaks[i][0]:                   #check if chromosome is the same as for the preceding one
                    if int(line[1]) <= peaks[i][2]:          #test if read is overlapping the preceding one
                        peaks[i][2] = int(line[2])           #replace end by end of the new read
                        if float(line[-1]) >= peaks[i][3]:   #test if Wm-score is higher as preceding one. Store WM-score of the window with highest WM-score.
                            peaks[i][3] = float(line[-1])
                    else:
                        i += 1
                        peaks[i] = [line[0], int(line[1]), int(line[2]), float(line[-1])]
                else:
                    i += 1
                    peaks[i] = [line[0], int(line[1]), int(line[2]), float(line[-1])]
            else:
                continue

 
    #file with merged peaks, sorted by descending WM-scores
    #format: chr7    127471196  127472363  pred1000001  WM-score  +
    o = open(outfile, 'w')
    peaklist = sorted(peaks.values(), key=lambda k: k[-1], reverse=True)

    i = 1000000
    for line in peaklist:
        o.write('\t'.join([str(line[0]), str(line[1]), str(line[2]), 'pred'+str(i), str(line[3]), '+\n']))
        i += 1
    o.close()

    return len(peaklist)


def execute(cf):
    """
    This component is a modified version of PeakMerger
    """

    infile = cf.get_input("infile")

    outfile = cf.get_output("outfile")
    intermdir = cf.get_output("intermediate")
    plotfile = cf.get_output("plotfile")
    logfile = cf.get_output("log_file")

    perlPATH = cf.get_parameter("perlPATH", "string")
    cutoff = cf.get_parameter("wmscore_cutoff", "int") #manually given WM-score cutoff, used when the value is non negative
 
    T1 = datetime.datetime.now()

    os.mkdir(intermdir)

    #BG cut-off estimation if not defined. For this run Phil's BG estimator.:
    if cutoff < 0:
        
        ##chr    start   stop    middle    fg_0
        #chr1    8750    9250    9000      4.2086659

        a = loadtxt(infile, usecols=[4], skiprows=1)

        x = sorted(a)
        y = log(arange(1, len(x)+1, 1))[::-1]

        cutoff = get_zco(x,y, plotfile)

        h = hist(x, sqrt(len(x)))
        close()
        b = [mean([h[1][i], h[1][i+1]]) for i in arange(len(h[0]))]
        figure()
        plot(b, h[0])

        xlabel('WM-scores')
        yscale('log')
        savefig(os.path.join(os.path.split(plotfile)[0], 'WMscore_hist.pdf'))
        close()


    else:
        #produce dummy plot
        figure()
        savefig(plotfile)
        close()


    T2 = datetime.datetime.now()

    print cutoff

    #Now filter infile with WM-score cut-off and merge overlapping windows:
    #outfile is a bed like file with predicted peak ID in 4th column and WM-score and strand in 5th and 6th.

    peaknum = merge_peaks(infile, outfile, cutoff)

    T3 = datetime.datetime.now()


    logtext = '\n'.join(['Used WM-score cut-off: %i' %cutoff,
                         '%i merged windows that made cut-off.' %(peaknum),
                         'Running time:',
                         '\tCut-off estimation: %s' %str(T2-T1),
                         '\tApplying cut-off and merge windows: %s' %str(T3-T2),
                         '\tOverall: %s' %str(T3-T1)
                         ])
    lf = open(logfile, 'w')
    lf.write(logtext)
    lf.close


    return 0

component_skeleton.main.main(execute)
