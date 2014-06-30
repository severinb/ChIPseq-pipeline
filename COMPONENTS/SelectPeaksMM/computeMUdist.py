#!/usr/bin/env python

import sys
from pylab import *

if __name__ == '__main__':

    mufile = sys.argv[1]
    outplot = sys.argv[2]

    a = loadtxt(mufile)

    dt = []
    for i in arange(len(a)):
        l = a[i]

        d = 0
        for j in arange(len(l)):
            for k in arange(j):
                d += sqrt((l[j] - l[k])**2)

        dt.append(d)


    dt = sorted(dt)[:-4]
    #d = sum([ [sqrt((a.T[i] - a.T[j])**2) for j in arange(i)] for i in arange(4)], axis=0)

    print dt

    hist(dt)

    savefig(outplot+'.pdf')
