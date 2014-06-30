#!/usr/bin/env python

import component_skeleton.main
from string import *
import sys, os
import subprocess
from datetime import datetime


def run(cmd):
    """
    given a command, it runs it
    """
  
    proc = subprocess.Popen (cmd,                      
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                shell=True
                         )
    stdout_value, stderr_value = proc.communicate()
    print stdout_value
    print stderr_value

    if proc.poll() > 0:
        sys.stderr.write ( "\nError\n" )
        print '\tstderr:', repr(stderr_value.rstrip())
        return False
    else:
        return True


def execute(cf):

    infile = cf.get_input("read_file") #optional input file that contains reads (I use this for example to bin the genome wide predicted WM scores...)

    outfile = cf.get_output("out_file")
    logfile = cf.get_output("BinReads_log")

    FGfiles_string = cf.get_parameter("FGfiles_string", "string")
    BGfiles_string = cf.get_parameter("BGfiles_string", "string")
    FGwinsize = cf.get_parameter("FGwinsize", "int")
    BGwinsize = cf.get_parameter("BGwinsize", "int")
    stepsize = cf.get_parameter("stepsize", "int")
    chrominfo = cf.get_parameter("chrominfo", "string")
    perlPATH = cf.get_parameter("perlPATH", "string")
    FMIpath = cf.get_parameter("FMIpath", "string")

    T1 = datetime.now()

    FGfiles = FGfiles_string.split()
    FGargs = []
    if len(FGfiles) != 0:
        for f in FGfiles:
            FGargs.append('-f ' + f)

    BGfiles = BGfiles_string.split()
    BGargs = []
    if len(BGfiles) != 0:
        for b in BGfiles:
            BGargs.append('-b ' + b)

    if infile:
        FGargs.append('-f ' + infile)

    cmd = ' '.join(['%s %s/soft.bc2/binReads-bedops.pl' %(perlPATH, FMIpath)] \
                       + FGargs + BGargs + \
                       ['-l %s' %chrominfo,
                        '--fgwin=%i' %FGwinsize,
                        '--bgwin=%i' %BGwinsize,
                        '--shift=%i' %stepsize,
                        '--outfile=%s' %outfile,
                        #'--sorted',
                        '--fmi-repo=%s' %FMIpath
                        ])

    print cmd

    pwd = os.getcwd()
    os.chdir(os.path.split(outfile)[0])

    run(cmd)

    os.chdir(pwd)


    T2 = datetime.now()
    time = 'Running time for read counter: ' + str(T2-T1) + '\n'
    lf = open(logfile, 'w')
    lf.write(time)
    lf.close

    return 0

component_skeleton.main.main(execute)
    
   
