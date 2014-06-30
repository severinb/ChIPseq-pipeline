#!/usr/bin/env python
import component_skeleton.main
import subprocess
import re, os, sys
from string import *
from datetime import datetime
import shutil

def execute(cf):
    #in_dir = cf.get_input("in_dir")
    out_dir = cf.get_output("out_dir")
    error_profiles_plot = cf.get_output("error_profiles_plot")
    fraction_aligned_plot = cf.get_output("fraction_aligned_plot")
    repeatedness_plot = cf.get_output("repeatedness_plot")
    logfile = cf.get_output("createErrorPlots_log")
    FMIid = cf.get_parameter("FMIid", "string")
    annoType = cf.get_parameter("annoType","string")
    FMIpath = cf.get_parameter("FMIpath", "string")
    perlPATH = cf.get_parameter("perlPATH", "string")

    T1 = datetime.now()

    os.mkdir(out_dir)
 
    #os.environ['PATH'] = '/import/bc2/soft/app/perl/5.10.1/Linux/bin:' + os.environ['PATH']
    #to be sure to use the right perl
    os.environ['PATH'] =  perlPATH + ':' + os.environ['PATH']

    ##change present working directory to out_dir. Otherwise .Rout files get stored in component's directory
    pwd = os.getcwd()
    os.chdir(out_dir)

    outfile = os.path.join(out_dir, FMIid)

    command = ' '.join([perlPATH,
                        FMIpath+'/soft/extractData.pl',
                        FMIid,
                        annoType,
                        'genome -f | ',
                        FMIpath+'/soft/createErrorPlots.pl',
                        '-o %s -' %outfile
                       ])

    proc = subprocess.Popen (command,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True
                             )

    stdout_value, stderr_value = proc.communicate()
    print 'STDOUT: ', stdout_value
    print 'STDERR: ', stderr_value
    if proc.poll() > 0:
        print '\tstderr:', repr(stderr_value.rstrip())
        return -1	



    #/import/bc2/soft/bin/R CMD BATCH
    proc2 = subprocess.Popen ('R CMD BATCH ' + outfile + '.R',
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              shell=True
                              )

    stdout_value2, stderr_value2 = proc2.communicate()
    print 'R STDOUT: ', stdout_value2
    print 'R STDERR: ', stderr_value2
    if proc2.poll() > 0:
        print '\tstderr:', repr(stderr_value2.rstrip())
        return -1	

    shutil.copyfile(outfile + '-error_profiles.pdf', error_profiles_plot)
    shutil.copyfile(outfile + '-fraction_aligned.pdf', fraction_aligned_plot)
    shutil.copyfile(outfile + '-repeatedness.pdf', repeatedness_plot)

    os.chdir(pwd)

    T2 = datetime.now()
    time = 'Running time for creating error plots: ' + str(T2-T1) + '\n'
    lf = open(logfile, 'a')
    lf.write(time)
    lf.close	

    return 0

component_skeleton.main.main(execute)
