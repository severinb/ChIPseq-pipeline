#!/usr/bin/env python
import component_skeleton.main
import os, re
from operator import itemgetter
from itertools import groupby
from subprocess import *
from string import *

def concatRepeats2Bed(repeatPath, repeatfile):
    """
    To filter out repeats I use bedtools. Therefore the repeats of all the chromosomes have to be in one bed file in format.
    repeat file format: 778     167     7       20      chr2L   1       154     -23011390       +       HETRP_DM        Satellite       Satellite       1519    1669    -203    1
    """

    o = open(repeatfile, 'w')
    repdir = os.listdir(repeatPath)
    for chrom in repdir:
        if chrom.startswith('chr') and chrom.endswith('.fa.out'):
            for line in open(os.path.join(repeatPath,chrom)):
                if line.strip():
                    try:  #do this because there can be header lines that start with strings. And there can be empty lines
                        t = line.strip().split()
                        i = int(t[0])
                        o.write('%s\t%s\t%s\t%s\n' %(t[4], t[5], t[6], t[8]))
                    except ValueError:
                        continue
    o.close()


def execute(cf):
    peaks = cf.get_input("peaks") #result file from RefinePeaks component. But any peaks in this format can be given
    interm = cf.get_output("intermediate")
    noRepPeaks = cf.get_output("out_file")
    log_file = cf.get_output("log_file")
    repeatPath = cf.get_parameter("repeatPath", "string")
    BedToolsPath = cf.get_parameter("BedToolsPath", "string")
    
    os.mkdir(interm)
    ##concatenate all chromosome repeat files into one.
    repeatfile = os.path.join(interm, 'repeats')
    concatRepeats2Bed(repeatPath, repeatfile)
            
    ##run intersectBed
    intersectfile = os.path.join(interm, 'intersects')
    #proc = subprocess.Popen('%s intersect -a %s -b %s -wb -wa > %s' %(os.path.join(BedToolsPath, 'bedtools'), peaks, repeatfile, intersectfile),
    cmd = ' '.join([
        BedToolsPath,
        'subtract',        
        '-a %s' % peaks,
        '-b %s' % repeatfile ])
    p1 = Popen(cmd, stdout=PIPE, shell=True)
    p2 = Popen("awk '{if(($3-$2)>=50) print $0}'", stdout=PIPE, stdin=p1.stdout, shell=True) # testing weather the sequences are at least 50 bp long
    with open(noRepPeaks, 'w') as file:
        for i in xrange(1000):   # selecting top 1000 sequences, this can be changed or even provided as an input paramter to the component 
            file.write(p2.stdout.readline())

    # cleaning up the intermediate directory
    cmd = 'gzip %s' % repeatfile
    os.system(cmd)
    
    return 0


component_skeleton.main.main(execute)
                                                                 
