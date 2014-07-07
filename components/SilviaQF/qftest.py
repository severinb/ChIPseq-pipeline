#!/usr/bin/env python
import os
import sys
import math

if __name__ == '__main__':

    bseq = 'AAACCTTTGTTTCCCATTCTGTTACCTNCTNCT'#TNT'
    bq = 'IIHHIIHIIIIIIHHIIHIGIIIEIII#DI#II'#I#I'
    min_quality = 20
    max_N = 2
    min_length = 25
    q_sum = 0
    curr_l = 0
    offset = 33
    report = "--"

    for i in range(len(bq)):
        q_sum += 10 * math.log10(10 ** ((ord(bq[i]) - offset) / 10.0) + 1)
        
        q_mean = q_sum / float(i+1)
        print q_mean
        if q_mean >= min_quality:
            curr_l = i+1
    if bseq[:curr_l].count("N") <= max_N:
        if curr_l >= min_length:
            report = "accepted"
        else:
            report = "not"
    else:
        report = "passes not"


    print 'report: ' +report
    print 'q_mean: ' +str(q_mean)
    print 'curr_l: ' +str(curr_l)
    print 'q_sum: ' +str(q_sum)
