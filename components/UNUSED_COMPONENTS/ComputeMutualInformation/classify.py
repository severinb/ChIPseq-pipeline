#!/usr/bin/env python

import os
from pylab import *

def GetMutInf(data, z_co, wm_co):
    """
    This function fits cut-offs on Z-score and WM-score by maximizing mutual information I.
    I(x,y) = sum_xy (P(x,y) * log(P(x,y)/(P(x)*P(y)))). 
    Binary Z and WM scores 11, 10, 01, 00. But cut-offs to make it binary are unknown. Maximize mutual information by fitting these cut-offs. 

    I(z,wm) = f++*log(f++/(f+.f.+)) + f+-*log(f+-/(f+.f.-)) + f-+*log(f-+/(f-.f.+)) + f--*log(f--/(f-.f.-))

    data is a matrix with a row for each region (500 bp windows) and the first column containing Z-scores, the second column containing WM-scores
    """

    def categorizeData(mat, z_co, wm_co):

        Zpidx = set(where(mat.T[0] >= z_co)[0])
        Znidx = set(where(mat.T[0] < z_co)[0])
        WMpidx = set(where(mat.T[1] >= wm_co)[0])
        WMnidx = set(where(mat.T[1] < wm_co)[0])

        Fpp = 1.0 + len(Zpidx - WMnidx)
        Fpn = 1.0 + len(Zpidx - WMpidx)
        Fnp = 1.0 + len(WMpidx - Zpidx)
        Fnn = 1.0 + len(set(arange(len(mat.T[0]))) - Zpidx - WMpidx)


        return array([[Fpp, Fpn], [Fnp, Fnn]])


    Fmat = categorizeData(data, z_co, wm_co) #get counts
    Fmat /= sum(Fmat) #get frequencies

    mi = 0.0
    for zrow in [0,1]:
        for wmcol in [0,1]:
            mi += Fmat[zrow, wmcol] * log2( Fmat[zrow, wmcol] / ( sum(Fmat[zrow, :]) * sum(Fmat[:, wmcol]) ) )

    # normalize mi by computing entropies
    P_z_p = sum(Fmat[0, :])
    P_z_n = sum(Fmat[1, :])

    P_wm_p = sum(Fmat[:, 0])
    P_wm_n = sum(Fmat[:, 1])

    H_z = -(P_z_p*log2(P_z_p) + P_z_n*log2(P_z_n))
    H_wm = -(P_wm_p*log2(P_wm_p) + P_wm_n*log2(P_wm_n))

    nmi = mi/min([H_z, H_wm])

    return mi, nmi


def main():
    """
    This component calls the program that computes WM scores over whole genome
    """

    ##Ports and parameters
    datafile = sys.argv[1]
    infileroot = sys.argv[2]
    outfileroot = sys.argv[3]

    jobid = int(os.environ['SGE_TASK_ID'])

    a = loadtxt(datafile, usecols=[1,2])

    o = open(outfileroot + str(jobid), 'w')

    for line in open(infileroot + str(jobid)):
        t = line.strip().split()
        zco = float(t[0])
        wmco = float(t[1])

        mutinf, norm_mutinf = GetMutInf(a, zco, wmco)

        o.write('%s\t%s\t%s\t%s\n' %(zco, wmco, mutinf, norm_mutinf))

    o.close()


if __name__ == '__main__':
    main()
