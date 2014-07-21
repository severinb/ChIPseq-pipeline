#!/usr/bin/env python
import component_skeleton.main
import os, re
from operator import itemgetter
from itertools import groupby
import subprocess
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
        '-b %s' % repeatfile,
        '|',
        'awk \'{if(($3-$2)>=50) print $0}\'',
        '|',
        'head -1000',
        '>',
        noRepPeaks,
        ])
    proc = subprocess.Popen(cmd,
                            stdout=subprocess.PIPE,
                            stderr= subprocess.PIPE,
                            shell=True
                            )

    # proc = subprocess.Popen('%s intersect -a %s -b %s -wb -wa > %s' %(BedToolsPath, peaks, repeatfile, intersectfile),
    #                         stdout=subprocess.PIPE,
    #                         stderr= subprocess.PIPE,
    #                         shell=True
    #                         )
    stdout_value, stderr_value = proc.communicate()
    print stdout_value
    print stderr_value

    if proc.poll() > 0:
        print '\tstderr:', repr(stderr_value.rstrip())
        return -1
    else:
        print '\nRepeat Region Intersection Done.\n'


    ##filter out repeats
    # PeakDict = {}

    ###initalize dictionary with all refined top scoring peaks
    # for line in open(peaks):
    #     tmp = line.strip().split()
    #     peak = '-'.join(tmp[:5])
    #     PeakDict[peak] = []

    ##format of intersect file: chr5   32444873    32445002  id    +     51.7691803279  0.000348359165405    chr5    32444926  32444968    GC_rich   21     +
    ##chr10   121356171       121356353       reg1000053.p1   6.747   +       chr10   121356151       121356176       +
    ##find peak overlapping repeats and put them into dictionary
    # for line in open(intersectfile):
    #     tmp = line.strip().split()
    #     peak = '-'.join(tmp[:5])
    #     PeakDict[peak].append(tmp[6:9])

    ##tooshort = 0   #counts number of peaks that were shorter than 20 bp
    # norepeats = 0  #counts number of peaks that don't contain repeats
    # givenout = 0   #counts number of given out peaks (long enough, ones with and without repeats)
    # onlyrepeat = 0 #counts how many peaks lie completely in a repeat

    # o = open(noRepPeaks, 'w')
    ##sort out repeats and then group continuos peaks that are longer than 50
    # for peak in PeakDict:
    #     p = peak.split('-')
    #     if not PeakDict[peak]:  #if there are no repeats in this peak, just print the peak, if it is longer than 20 (window length of PhyloGibbs)
    #         norepeats += 1
    #         o.write('\t'.join([p[0],p[1],p[2], p[3], p[4], '+\n']))
    #         givenout += 1
            ##restrict length at RunPhyloGibbs step
            ## if int(p[2])-int(p[1]) > 20:
            ##     o.write('\t'.join([p[0],p[1],p[2], p[3], p[4], '+\n']))
            ##     givenout += 1
            ## else:
            ##     tooshort += 1
            ## continue

        # else:                   #else sort out repeats
        #     prange = set(range(int(p[1]),int(p[2])+1))
        #     for r in PeakDict[peak]:
        #         rrange = set(range(int(r[1]), int(r[2])+1))
                ##make a list with bp positions that arent repeats
            #     prange -= rrange

            # newranges = list(prange)
            # newranges.sort()
            # if len(newranges) < 1:
            #     onlyrepeat += 1
            ##group continuous bp positions into regions
            # ij = 1
            # for k, g in groupby(enumerate(newranges), lambda (i,x):i-x):
            #     a = map(itemgetter(1), g)
            #     o.write('\t'.join([p[0],str(a[0]),str(a[-1]), p[3] + '.rep' + str(ij), p[4],  '+\n']))
            #     givenout += 1
            #     ij += 1

                ##restrict length at RunPhyloGibbs step
                ## if len (a) > 20:
                ##     o.write('\t'.join([p[0],str(a[0]),str(a[-1]), p[3] + '.rep' + str(ij), p[4],  '+\n']))
                ##     givenout += 1
                ##     ij += 1
                ## else:
                ##     tooshort += 1


    # o.close()
                                                                                                                                                
    ##os.system('rm %s' %repeatfile)

    # l = open(log_file, 'w')

    # text = '\n'.join(['-Number of input peaks: %i' %len(PeakDict.keys()),
    #                   '-Number of given out peaks: %i' %givenout,
    #                   #'-Number of too short peaks (shorter than 20bp): %i' %tooshort,
    #                   '-Number of non repeat containing peaks: %i' %norepeats,
    #                   '-Number of peaks that lay completely in a repeat: %i' %onlyrepeat])

    # l.write(text)
    # l.close()

    return 0


component_skeleton.main.main(execute)
                                                                 
