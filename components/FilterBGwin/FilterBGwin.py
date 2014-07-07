#!/usr/bin/env python
import component_skeleton.main
import subprocess
import os
from string import *
import datetime
from pylab import *

def execute(cf):

    infile = cf.get_input("infile")

    outfile = cf.get_output("outfile") #instead of outdir with output file with name all_counts
    intermdir = cf.get_output("intermediate")
    plotfile = cf.get_output("plotfile")
    logfile = cf.get_output("log_file")

    perlPATH = cf.get_parameter("perlPATH", "string")
    bg_cutoff = cf.get_parameter("bg_cutoff", "int") #manually given bg cutoff, used when the value is non negative
 
    T1 = datetime.datetime.now()

    #os.mkdir(outdir)
    os.mkdir(intermdir)
    #outfile = os.path.join(outdir, 'all_counts')

    #check how many bg replicates there are
    f = open(infile)
    header = f.readline()
    hlist = header.strip().split()
    BGnum = 0
    for i in hlist:
        if i.startswith('bg'):
            BGnum += 1

    print 'BGnum: ', BGnum

    #BG cut-off estimation if not defined. For this run Phil's BG estimator.:
    if bg_cutoff < 0:

        #because I don't know how to program in perl, I first create a file that looks similar to saeed's readcount output (with summed counts in last column).
        tmpfile = os.path.join(intermdir, 'tmp')
        t = open(tmpfile, 'w')
        
        ##chr    start   stop    middle  fg_0    fg_1    fg_2    bg_0
        #chr1    8750    9250    9000    0       0       0       4.2086659
        for line in f:
            l = line.strip().split()
            BGsum = 0
            for i in range(BGnum):
                BGsum += float(l[-BGnum+i])
            t.write('\t'.join(l[:4] + [str(BGsum)])+'\n')

        t.close()

        cofile = os.path.join(intermdir, 'cutoff_file')

        print datetime.datetime.now(), 'starting estimation'
        # command = ' '.join([perlPATH,
        #                     'easy_bgcutoff_smoothingS.pl',
        #                     tmpfile,
        #                     intermdir,
        #                     plotfile])

        command = ' '.join(['./newfilter.py',
                            tmpfile,
                            intermdir,
                            plotfile,
                            cofile])

        proc = subprocess.Popen (command,
                                 stdout=subprocess.PIPE,
                                 stderr= subprocess.PIPE,
                                 shell=True
                                 )
        stdout_value, stderr_value = proc.communicate()
        print stdout_value
        print stderr_value

        if proc.poll() > 0:
            print '\tstderr:', repr(stderr_value.rstrip())
            return -1

        os.system('rm %s' %tmpfile)

        bg_cutoff = float(open(cofile).read())

    else:
        #produce dummy plot
        figure()
        savefig(plotfile)
        close()

    T2 = datetime.datetime.now()

    print bg_cutoff


    #Now filter infile with BG cut-off:
    f = open(infile)
    o = open(outfile, 'w')

    header = f.readline()
    o.write(header)

    fout = 0 #counts number of filtered out windows

    ei = 0 #count for all windows
    for line in f:
        ei += 1
        t = line.strip().split()
        BGsum = 0
        for i in range(BGnum):
            BGsum += float(t[-BGnum+i])
        if BGsum < bg_cutoff:
            o.write(line)
        else:
            fout += 1
            continue

    f.close()
    o.close()

    T3 = datetime.datetime.now()

    logtext = '\n'.join(['Used BG cut-off: %i' %bg_cutoff,
                         '%i of %i input windows were filtered out (%s percent).' %(fout, ei, ((float(fout)/ei)*100))
                         #'Running time:',
                         #'\tBG cut-off estimation: %s' %str(T2-T1),
                         #'\tFiltering by BG cut-off: %s' %str(T3-T2),
                         #'\tOverall: %s' %str(T3-T1)
                         ])
    lf = open(logfile, 'w')
    lf.write(logtext)
    lf.close


    return 0

component_skeleton.main.main(execute)
