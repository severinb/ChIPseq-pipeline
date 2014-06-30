#!/usr/bin/env python

import sys
import os
from string import *
import component_skeleton.main
import datetime
from pylab import *


def get_zco(x, y, outroot):

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
    idx1 = where(dy>2)[0][0]
    plot([x[idx1], x[idx1]], [0,max(y)])

    # 3. fit a second line from this index on. Try to find the steepest line in some region
    slopes = []
    inters = []
    minslope = 10
    for i in arange(size/50,size/5, 1):
        r2 = i

        X = vstack([x[idx1:idx1+r2], np.ones(len(x[idx1:idx1+r2]))]).T
        params = lstsq(X, y[idx1:idx1+r2])

        fx2 = x*params[0][0]+params[0][1]

        slopes.append(params[0][0])
        inters.append(params[0][1])


    slope = min(slopes)
    inte = inters[argmin(slopes)]
    upperb = arange(2,size/5, 1)[argmin(slopes)]
    print x[idx1+upperb]

    fx2 = x*slope+inte

    plot([x[idx1+upperb], x[idx1+upperb]], [0,max(y)])
    plot(x[idx1:idx1+upperb],fx2[idx1:idx1+upperb])

    # 3. again compute deviation
    dy2 = fx2 - y
    try: # it can happen that there is no enrichment, and then no deviation will be found that is lower than -.5. So just use fixed 3.0 cut-off
        idx2 = where(dy2[where(x >= x[idx1+upperb])] <=- 0.5)[0][0]
        co = x[idx2 + idx1 + upperb] # only allow the cut-off to be higher than upper bound
    except IndexError:
        co = 3.0
        print 'Warning: Could not find a cut-off automatically. Use fixed cut-off of 3.0'


    plot([co,co],[0,12], label='Z cut-off')

    legend()
    savefig('%s.pdf' %(os.path.join(outroot, 'Z-co-fit')))


    return round(co,2)


def merge_peaks(infile, outfile, zfile, z_cutoff, topPeaks, idfile):
    """This function filters peaks by a cut off and merges overlapping peaks. 
       Stores highest z value for each peak.
       Returns list of top 1000 peaks sorted by chromosome.
       chrom start end strand(+)	 
    """
    
    f = open(infile, 'r')
    peaks = {}
    i = 1

    #initializing by finding the first read above z_cutoff
    init = True
    while init:
        line = f.readline().split()
        zval = float(line[-2])
        if zval >= z_cutoff:
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
            if float(line[-2]) >= z_cutoff:
                if line[0] == peaks[i][0]:                   #check if chromosome is the same as for the preceding one
                    if int(line[1]) <= peaks[i][2]:          #test if read is overlapping the preceding one
                        peaks[i][2] = int(line[2])           #replace end by end of the new read
                        if float(line[-2]) >= peaks[i][3]:   #test if z value is higher as preceding one. Store z value of the window with highest z value.
                            peaks[i][3] = float(line[-2])
                    else:
                        i += 1
                        peaks[i] = [line[0], int(line[1]), int(line[2]), float(line[-2])]
                else:
                    i += 1
                    peaks[i] = [line[0], int(line[1]), int(line[2]), float(line[-2])]
            else:
                continue

 
    #file with merged peaks, sorted by descending z-values
    #format: chr7    127471196  127472363  reg1000001  Z  +
    oso = open(zfile, 'w')
    idf = open(idfile, 'w')
    peaklist = sorted(peaks.values(), key=lambda k: k[-1], reverse=True)
    i = 1000000
    for line in peaklist:
        oso.write('\t'.join([str(line[0]), str(line[1]), str(line[2]), 'reg'+str(i), str(line[3]), '+\n']))
        idf.write('%s_%s_%s\t%s\n' %(str(line[0]), str(line[1]), str(line[2]), 'reg'+str(i)))
        i += 1
    oso.close()
    idf.close()

    #choose 500 top peaks and sort them by chromosome and position
    #add '+' for strandfor each peak and print 'chrom start end +'
    i = 0
    o = open(outfile, 'w')
    for line in open(zfile):
        if i >= topPeaks:
            break
        else:
            o.write(line)
            i += 1
    o.close()


def execute(cf):
    in_dir = cf.get_input("in_dir")
    outfile = cf.get_output("out_file")
    allpeaks = cf.get_output("allpeaks")
    outplot = cf.get_output("revcum")
    logfile = cf.get_output("PeakMerger_log")
    IDfile = cf.get_output("IDfile")
    zco_file = cf.get_output("Z_cutoff_file")
    z_co = cf.get_parameter("z_cutoff", "float")
    topPeaks = cf.get_parameter("topPeaks", "int")

    T1 = datetime.datetime.now()

    infile = os.path.join(in_dir, 'outzvals')

    # load Z-values to plot and get Z-cutoff
    a = loadtxt(infile, usecols=[-2])
    x = array(sorted(a))
    y = log(arange(1, len(x)+1, 1))[::-1]

    if z_co < 0:
        z_co1 = get_zco(x, y, os.path.split(outfile)[0])

        # in case fitting goes wrong, constrain it:
        z_co = max(min(z_co1, 6.),2.)

    # plot reverse cumulative
    figure()
    plot(x,y)
    plot([z_co, z_co], [min(y), max(y)], label='chosen cut off %s' %z_co)
    ylim(0,int(max(y)+1))
    xlabel("Window Z-score")
    ylabel("log(Number of Windows)")
    legend()
    title('Reverse Cumulative Distribution of Window Z-scores (log-counts)')
    savefig(outplot)
    close()


    merge_peaks(infile, outfile, allpeaks, z_co, topPeaks, IDfile)

    T2 = datetime.datetime.now()
    l = open(logfile, 'a')
    l.write('Z value cut off: %s \n' %z_co)
    #l.write("Running time: " + str(T2-T1) + '\n')
    l.close()
 
    o = open(zco_file, 'w')
    o.write('%s\n' %z_co)
    o.close()

    return 0


component_skeleton.main.main(execute)
