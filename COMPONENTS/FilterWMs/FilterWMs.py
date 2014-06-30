#!/usr/bin/env python

import component_skeleton.main
import sys, os
from pylab import *


def filterQuarterWMs(wms):
    """
    Sometimes WMs just contain 1 everywhere (I produce them in RunMotevo).
    This function filters those out.
    """

    goodwms = []
    badwms = []

    for wm in wms:
        name, mat, matrev = wm2mat(open(wm))
        if mat.shape[0] * mat.shape[1] == len(where(mat == 0.25)[0]):
            badwms.append(wm)
        else:
            goodwms.append(wm)


    return goodwms, badwms


def wm2mat(wmFile,ps=0.5):
    """ Read swiss regulon style WM file """
    M = False
    ok = False
    lines = wmFile.read().splitlines()
    for line in lines:
        if not (line.startswith("//") or line==""):
            if ok:
                elm = reshape(array(map(float,line.split()[1:5])),(1,4))
                if M is False:
                    M = elm
                else:
                    M = vstack((M,elm))
            else:
                if line[0:2]=='NA':
                    na,name = line.split()
                elif line[0:2]=='PO' or line[0:2]=='P0': # WM start after PO/0 line
                    ok = True
                else:
                    pass
    # add pseudo count
    M += ps
    # convert counts to frequencies
    M = M/M.sum(axis=1)[:,newaxis]
    # also build the reverse complement matrix
    Mrev = vstack((M[::-1,3],M[::-1,2],M[::-1,1],M[::-1,0])).T

    return(name,M,Mrev)


def execute(cf):
    """
    This component reduces candidate WMs.
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
    outdir = cf.get_output("WMdir")
    log_file = cf.get_output("log_file")

    os.mkdir(outdir)

    wms = [i for i in[WM1, WM2, WM3, WM4, WM5, WM6, WM7, WM8, WM9, WM10, WM11, WM12, WM13, WM14, WM15, WM16, WM17, WM18, WM19, WM20] if i]

    goodwms, badwms = filterQuarterWMs(wms)


    for i, wm in enumerate(goodwms):
        os.system('cp %s %s/WM_%i' %(wm, outdir, i+1))


    print badwms

    l = open(log_file, 'w')
    l.write('%s good WMs, %s bad WMs' %(len(goodwms), len(badwms)))
    l.close()


    return 0

component_skeleton.main.main(execute)
