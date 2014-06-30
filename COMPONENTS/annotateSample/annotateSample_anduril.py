#!/usr/bin/env python
import component_skeleton.main
import subprocess
import re, os, sys
from datetime import datetime

def execute(cf):
    FMIid = cf.get_parameter("FMIid", "string")
    annoType = cf.get_parameter("annoType", "string")
    FMIpath = cf.get_parameter("FMIpath", "string")
    logfile = cf.get_output("annotateSample_log")
    perlPATH = cf.get_parameter("perlPATH", "string")

    T1 = datetime.now()

    #I need this because annotateSample.pl calls other perl programs (e.g. aln2.3 or so))
    os.environ['PATH'] = os.path.split(perlPATH)[0] + ':' + os.environ['PATH']

    FMIsoft = os.path.join(FMIpath, 'soft')

    print "submitting %s \n" %(FMIid)
    proc = subprocess.Popen (perlPATH + " " + FMIsoft + "/annotateSample.pl %s %s" %(FMIid, annoType),
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True
                             )
    stdout_value, stderr_value = proc.communicate()
    print stdout_value
    print stderr_value

    if proc.poll() > 0:
        print '\tstderr:', repr(stderr_value.rstrip())
        return -1

    T2 = datetime.now()
    l = open(logfile, 'w')
    log = 'Running time for sample annotation (mapping): ' + str(T2-T1)
    l.write(log)
    l.close()

    return 0


component_skeleton.main.main(execute)
