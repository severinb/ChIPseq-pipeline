#!/usr/bin/env python

import sys
import os
from string import *
import numpy
import math
import datetime
from pylab import *


def get_zco(x, y, outroot):

    figure()

    #to make a stable analysis, sample the data linearly spaced in x
    size = 1000
    xn = linspace(min(x), min(2000, max(x)), size)
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

    #plot(x, fx2)
    plot([x[idx1+upperb], x[idx1+upperb]], [0,max(y)])
    plot(x[idx1:idx1+upperb],fx2[idx1:idx1+upperb])

    savefig('%s.pdf' %(os.path.join(outroot, 'co-fit'))) # intermediate save of plot

    # 3. again compute deviation
    dy2 = fx2 - y
    idx2 = where(dy2[where(x >= x[idx1+upperb])] <= -0.5)[0][0]

    co = x[idx2 + idx1 + upperb] #cut-off has to be higher than upperb... makes it more stable

    plot([co,co],[0,12], label='cut-off: %.3f' %co)

    xlim([0, min(max(co*3,1000), max(x))])
    legend()
    savefig('%s.pdf' %(os.path.join(outroot, 'co-fit')))

    return int(co)


def main():

    infile = sys.argv[1]
    interm = sys.argv[2]
    outplot = sys.argv[3]
    cofile = sys.argv[4]

    # load Z-values to plot and get Z-cutoff
    a = loadtxt(infile, usecols=[-1])
    x = array(sorted(a))
    y = log(arange(1, len(x)+1, 1))[::-1]

    co = get_zco(x, y, interm)

    # plot reverse cumulative
    figure()
    plot(x,y)
    plot([co, co], [min(y), max(y)], label='chosen cut-off %s' %co)
    xlim([0,min(max(3*co,1000), max(x))])
    ylim(0,int(max(y)+1))
    xlabel("Read Count")
    ylabel("log-Number of Background Windows")
    legend()
    title('Reverse Cumulative Distribution of Background Read Counts')
    savefig(outplot)
    close()


    o = open(cofile, 'w')
    o.write(str(co))
    o.close()


    return 0


if __name__ == "__main__":
    main()
