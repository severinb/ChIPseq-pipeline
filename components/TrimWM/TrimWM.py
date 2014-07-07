#!/usr/bin/env python
import component_skeleton.main
import os, re
from string import *
import datetime, time
from pylab import *


def trimWM(wmfile, co, outwm):
    """
    This function takes a WM file and trims the edges until one position has an information content of above co.
    """    

    f = open(wmfile, 'r')

    ls = f.readlines()
    lines = [i for i in ls if i.strip()] #sort out empty lines (i.e. '\n' lines)

    header = lines[:3]
    header[1] = 'NA\t%s\n' %os.path.split(outwm)[1]

    footer = lines[-1]

    cols = lines[3:-1]

    start = 0
    stop = len(cols) -1

    def informCont(a,c,g,t):
        pc = 0.001
        tot = a + c + g + t + 4*pc

        #return 2 - sum([((-i/tot)+pc) * log2((i/tot)+pc) for i in [a,c,g,t]])
        return 2 - sum([(-(i+pc)/tot) * log2((i+pc)/tot) for i in [a,c,g,t]])
        

    for i in arange(len(cols)):

        t = cols[i].split()
        A = float(t[1])
        C = float(t[2])
        G = float(t[3])
        T = float(t[4])

        if informCont(A,C,G,T) < co:
            start = i
            continue
        else:
            start = i
            break


    for i in arange(len(cols))[::-1]:

        t = cols[i].split()
        try:
            A = float(t[1])
        except IndexError:
            print t
            print cols
            print i
            print wmfile
        C = float(t[2])
        G = float(t[3])
        T = float(t[4])

        if informCont(A,C,G,T) < co:
            stop = i + 1
            continue
        else:
            stop = i + 1
            break


    o = open(outwm, 'w')
    for i in header:
        o.write(i)

    for i,j in enumerate(cols[start:stop]):
        t = j.split()
        o.write('\t'.join([str(i+1).zfill(2)] + t[1:]) + '\n')

    o.write(footer)
    o.close()
    

    #return how much WM was trimmed from left and from right
    #if WM was trimmed everywhere, i.e. the WM had no informative columns, then return 0,0 to not give out the WM.
    if start >= stop:
        return 0, 0
    return start, (len(cols)-stop)


def execute(cf):
    """
    This component trims all given WMs from left and right edges until information content of a column gets above cut-off.
    """

    ##Ports and parameters
    WM1 = cf.get_input("WM1")
    WM2 = cf.get_input("WM2")
    WM3 = cf.get_input("WM3")
    WM4 = cf.get_input("WM4")
    WM5 = cf.get_input("WM5")
    WM6 = cf.get_input("WM6")
    WM7 = cf.get_input("WM7")
    WM8 = cf.get_input("WM8")
    WM9 = cf.get_input("WM9")
    WM10 = cf.get_input("WM10")
    WM11 = cf.get_input("WM11")
    WM12 = cf.get_input("WM12")
    WM13 = cf.get_input("WM13")
    WM14 = cf.get_input("WM14")
    WM15 = cf.get_input("WM15")
    WM16 = cf.get_input("WM16")
    WM17 = cf.get_input("WM17")
    WM18 = cf.get_input("WM18")
    WM19 = cf.get_input("WM19")
    WM20 = cf.get_input("WM20")
    WMdir = cf.get_input("WMdir")

    outdir = cf.get_output("outdir")
    log_file = cf.get_output("log_file")

    co = cf.get_parameter("information_cutoff", "float")

    os.mkdir(outdir)
    
    WMs = [i for i in[WM1, WM2, WM3, WM4, WM5, WM6, WM7, WM8, WM9, WM10, WM11, WM12, WM13, WM14, WM15, WM16, WM17, WM18, WM19, WM20] if i]

    if WMdir:
        WMs += [os.path.join(WMdir, wm) for wm in  os.listdir(WMdir)]

    WMtrims = []

    for i, WM in enumerate(WMs):
        wm = os.path.split(WM)[1]
        outwm = os.path.join(outdir, '%s.trimmed' %wm)
        left, right = trimWM(WM, co, outwm)
        if left == 0 and right == 0:
            os.system('rm %s' %outwm)
        else:
            WMtrims.append((WM, os.path.split(outwm)[1], left, right))

    l = open(log_file, 'w')
    l.write('Weight matrix trimming done with column information content cut-off of %s:\nTrimmed WMs are stored here: %s\n' %(co,outdir))

    for wm in WMtrims:
        l.write('\t%s was trimmed to %s %i columns from left and %i columns from right\n' %(wm[0], wm[1], wm[2], wm[3]))
        
    l.close()


    return 0


component_skeleton.main.main(execute)
                                                                 
