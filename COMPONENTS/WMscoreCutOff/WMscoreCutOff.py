#!/usr/bin/env python
import component_skeleton.main
import subprocess
import os
from string import *
import datetime
from pylab import *

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
        
        ##chr    start   stop    middle    bg_0
        #chr1    8750    9250    9000      4.2086659

        print datetime.datetime.now(), 'starting estimation'

        command = ' '.join([perlPATH,
                            'easy_bgcutoff_smoothingS.pl',
                            infile,
                            intermdir,
                            plotfile])

        proc = subprocess.Popen (command,
                                 stdout=subprocess.PIPE,
                                 stderr= subprocess.PIPE,
                                 shell=True
                                 )
        stdout_value, stderr_value = proc.communicate()
        print stdout_value
        print stderr_value

        if proc.poll() > 0:
            print '\tstderr:', repr(stderr_value1.rstrip())
            return -1

        cofile = os.path.join(intermdir, 'stats_bg_cutoff')
        cutoff = float(open(cofile).readline().strip())

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
