#!/usr/bin/env python
import component_skeleton.main
import subprocess
import os
from string import *
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

    plot([co,co],[0,12], label='Z cut-off')

    legend()
    savefig('%s.pdf' %(os.path.join(outroot, 'Z-co-fit')))


    return round(co,2)


def execute(cf):

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

        print datetime.datetime.now(), 'starting estimation'

        a = loadtxt(infile, usecols=[4], skiprows=1)

        x = sorted(a)
        y = log(arange(1, len(x)+1, 1))[::-1]

        cutoff = get_zco(x,y, intermdir)

    else:
        #produce dummy plot
        figure()
        savefig(plotfile)
        close()

    T2 = datetime.datetime.now()

    print cutoff

    #Now filter infile with WM-score cut-off:
    #outfile is a bed like file with middle coordinate in 4th column and WM-score and strand in 5th and 6th.
    f = open(infile)
    o = open(outfile, 'w')

    header = f.readline()
    o.write('#chr\tstart\tstop\tmiddle\tWM-score\tstrand\n')

    fout = 0 #counts number of above cut-off windows

    ei = 0 #count for all windows
    for line in f:
        ei += 1
        t = line.strip().split()
        if float(t[4]) >= cutoff:
            o.write(line.strip()+'\t+\n')
            fout += 1
        else:
            continue

    f.close()
    o.close()

    T3 = datetime.datetime.now()

    logtext = '\n'.join(['Used WM-score cut-off: %i' %cutoff,
                         '%i of %i input windows made cut-off (%s percent).' %(fout, ei, ((float(fout)/ei)*100)),
                         'Running time:',
                         '\tCut-off estimation: %s' %str(T2-T1),
                         '\tApplying cut-off: %s' %str(T3-T2),
                         '\tOverall: %s' %str(T3-T1)
                         ])
    lf = open(logfile, 'w')
    lf.write(logtext)
    lf.close


    return 0

component_skeleton.main.main(execute)
