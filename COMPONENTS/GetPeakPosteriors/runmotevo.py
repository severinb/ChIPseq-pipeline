#!/usr/bin/env python

import subprocess
import os, re
from string import *
from pylab import *
import sys
import pickle

    
def runMotevo(motevo_path, seqs, params, WM):
    """
    runs Motevo
    """
    
    pwd = os.getcwd()
    os.chdir(os.path.split(params)[0])

    print '\nrun Motevo'
    proc = subprocess.Popen(motevo_path + ' %s %s %s' %(seqs, params, WM),
                            stdout=subprocess.PIPE,
                            stderr= subprocess.PIPE,
                            shell=True
                            )

    stdout_value, stderr_value = proc.communicate()

    os.chdir(pwd)
    
    if proc.poll() > 0:
        print '\tstderr:', repr(stderr_value.rstrip())
        return -1
    else:
        return 0


def getDicts(sites, statsfile, regcov_dir):
    """
    sitesDict: here IDs are just IDs of the region, thus no .p1 or .rep1 etc...
    13-22 - 0.573747 BestWM hg19_chr7_915000_916000_reg1000004_9.279_+
    This function returns a dictionary: reg1000004: [(13-22, -, 0.573747), (13-22, - ,0.573747) ...]

    IDstats: store all stats for one region under region ID key
    statsfile: reg1000087.p1   4.393   1.77485636798e-05       409.430 34.822 (height, rmsd, mu, sigma)
    reg1000087: [[4.393, 1.77485636798e-05, 409.430, 34.822], ]

    IDcoords:
    regcovfile: chr10   121356000       121357250       reg1000053      1       0.5
    chr10_121356000_121357250: reg1000053
    """


    sitesDict = {}
    for line in open(sites): #e.g. 1500-1511 + 0.001162 BestWM hg19_chr13_99962000_99964000_reg1002740_14.9861052478_+
        t = line.strip().split()
        regID = t[4].split('_')[-3]
        try:
            sitesDict[regID].append( (t[0], t[1], float(t[2])) )
        except KeyError:
            sitesDict[regID] = []
            sitesDict[regID].append( (t[0], t[1], float(t[2])) )


    IDstats = {} #the peakstats file from SelectPeaks or from NormalizeHieghts contains stats for ALL fitted regions of PeakMerger outfile
    for line in open(statsfile):
        t = line.strip().split()
        try:
            IDstats[t[0].split('.')[0]].append(map(float, t[1:]) + [t[0]])
        except KeyError:
            IDstats[t[0].split('.')[0]] = []
            IDstats[t[0].split('.')[0]].append(map(float, t[1:]) + [t[0]])


    IDcoords = {}
    for fi in os.listdir(regcov_dir):
        f = open(os.path.join(regcov_dir, fi))
        t = f.readline().split()
        f.close()
        IDcoords['_'.join(t[:-3])] = t[-3]


    return sitesDict, IDstats, IDcoords



def main():
    """
    This script runs motevo to predict sites. It also creates dictionaries.
    """

    motevo_path = sys.argv[1]
    seqs = sys.argv[2]
    parameter_file =  sys.argv[3]
    WM =  sys.argv[4]

    pickled_sitesdict = sys.argv[5]
    pickled_idstats = sys.argv[6]
    pickled_idcoords = sys.argv[7]

    sitesfile = sys.argv[8]
    statsfile = sys.argv[9]
    regcov_dir = sys.argv[10]


    runMotevo(motevo_path, seqs, parameter_file, WM)

    #get true sites and a post dict that contains non-refined coordinates as keys and lists of all posteriors as values
    sitesDict, IDstats, IDcoords = getDicts(sitesfile, statsfile, regcov_dir)


    pickle.dump(sitesDict, open(pickled_sitesdict, 'w'))
    pickle.dump(IDstats, open(pickled_idstats, 'w'))
    pickle.dump(IDcoords, open(pickled_idcoords, 'w'))


if __name__ == '__main__':
    main()
