#!/usr/bin/env python
import component_skeleton.main
import subprocess
import re, os
from datetime import datetime
from string import *


def execute(cf):
    infile = cf.get_input("in_file")
    logfile = cf.get_output("importSample_log")
    desc = cf.get_parameter("desc", "string")
    FMIpath = cf.get_parameter("FMIpath", "string")
    FMIid = cf.get_parameter("FMIid", "string")
    perlPATH = cf.get_parameter("perlPATH", "string")
    FMI_output_dir = cf.get_parameter("FMI_output_dir", "string")

    T1 = datetime.now()

    #to be sure to use the right perl
    os.environ['PATH'] =  perlPATH + ':' + os.environ['PATH']

    print 'Importing %s\n' %FMIid

    #The FMIid can be different from the infile name (with stripped unique.fa), if there is a rename tag given.
    #-v option disables the filtering (this we have already done in Silvia's quality filter)
    proc = subprocess.Popen(perlPATH + ' ' + FMIpath+"/soft/importSample.pl %s %s \'%s\' %s -f fasta -v" % (infile, FMIid, desc, FMI_output_dir),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=True
                            )

    stdout_value, stderr_value = proc.communicate()
    print stdout_value
    print stderr_value

    if proc.poll() != 0 and proc.poll() != 1:
        print '\tstderr:', repr(stderr_value.rstrip())
        return -1
    
    T2 = datetime.now()
    time = 'Running time for sample import: ' + str(T2-T1) + '\n'
    lf = open(logfile, 'a')
    lf.write(time)
    lf.close	

    return 0


component_skeleton.main.main(execute)
